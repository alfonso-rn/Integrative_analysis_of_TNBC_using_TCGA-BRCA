# 6. Identification of TFBS affected by mutations in DEGs ----

# Intronic SNPs identified in differentially expressed genes (DEGs) from TNBC
# samples were evaluated for their potential impact on transcription factor
# binding sites (TFBS). SNPs were prioritized from the processed MAF dataset and
# analysed with motifbreakR using HOCOMOCO position weight matrices to predict
# motif gains, strengthening, losses, or weakening events.
#
# Strong predicted TFBS effects were then overlapped with prioritized
# self-targeting cis-regulatory modules (CRMs) to identify regulatory regions
# where intronic mutations may alter transcription factor binding. Summary

library(dplyr)
library(motifbreakR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

maf_intronic_TNBC <- readRDS("data/processed/maf_intronic_TNBC.rds")
selfTargetingCRMs <- readRDS("data/processed/selfTargetingCRMs.rds")

# 6.2. Prediction of SNPs effects on TFBS motifs ----

# Priorization of SNPs
variants_snv <- maf_intronic_TNBC %>%
  filter(Variant_Type == "SNP") %>%
  select(
    Tumor_Sample_Barcode,
    gene = Hugo_Symbol,
    chrom = Chromosome,
    pos = Start_Position,
    end = End_Position,
    ref = Reference_Allele,
    alt = Tumor_Seq_Allele2,
    variant_type = Variant_Type,
    expression_status
  )

cat("SNPs mutations for TFBS analysis:", nrow(variants_snv))
cat("Mutations out of TFBS analysis:", paste(
    names(table(maf_intronic_TNBC$Variant_Type[maf_intronic_TNBC$Variant_Type != "SNP"])),
    table(maf_intronic_TNBC$Variant_Type[maf_intronic_TNBC$Variant_Type != "SNP"]),
    sep = " = ",
    collapse = ", "
  )
)

snv_bed <- variants_snv %>%
  transmute(
    chrom = chrom,
    chromStart = pos - 1L,   # 0-based BED coordinates
    chromEnd = pos,
    name = paste(chrom, pos, ref, alt, sep = ":"),
    score = 0,
    strand = "."
  )

write.table(snv_bed, "data/processed/snv.bed", sep = "\t",
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

snps_mb <- snps.from.file("data/processed/snv.bed",
  search.genome = BSgenome.Hsapiens.UCSC.hg38,
  format = "bed"
)

data(hocomoco) # HOCOMOCO database

mb_results <- motifbreakR(
  snpList = snps_mb,
  pwmList = hocomoco,
  threshold = 0.85,
  method = "ic",
  show.neutral = FALSE,
  BPPARAM = BiocParallel::SerialParam()
)

mb_results_df <- as.data.frame(
  mb_results,
  row.names = NULL
) %>%
  as_tibble() %>%
  mutate(SNP_id = names(mb_results), .before = 1)

mb_results_df$TFBS_effect <- NA_character_

for (i in seq_len(nrow(mb_results_df))) {
  
  if (mb_results_df$alleleDiff[i] > 0) {
    mb_results_df$TFBS_effect[i] <- "gain_or_strengthening"
    
    } else if (mb_results_df$alleleDiff[i] < 0) {
      mb_results_df$TFBS_effect[i] <- "loss_or_weakening"
      
      } else { 
        mb_results_df$TFBS_effect[i] <- "no_change"
      }
  }

cat("Number SNPs with at least one predicted TFBS effect:", n_distinct(mb_results_df$SNP_id))
cat("Total TFBS predictions of motifbreakR:", nrow(mb_results_df))

# Summary of predicted TFBS effects
tfbs_effect_summary <- mb_results_df %>%
  count(TFBS_effect, effect) %>%
  mutate(percentage = round(100 * n / sum(n), 2))

write.csv(tfbs_effect_summary, "results/tables/tfbs_effect_summary.csv", row.names = FALSE)

# Summary of TFs with strong predicted effects
top_tf_strong <- mb_results_df %>%
  filter(effect == "strong") %>%
  group_by(geneSymbol) %>%
  summarise(
    strong_predictions = n(),
    affected_SNPs = n_distinct(SNP_id),
    gains = sum(TFBS_effect == "gain_or_strengthening", na.rm = TRUE),
    losses = sum(TFBS_effect == "loss_or_weakening", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(strong_predictions))

write.csv(top_tf_strong, "results/tables/top_tf_strong.csv", row.names = FALSE)

cat("Distinct TFs with strong predicted effects:", nrow(top_tf_strong))

# 6.2. Overlap between affected TFBS and prioritized CRMs ----

# Keep strong motifbreakR effects
mb_strong <- mb_results_df[mb_results_df$effect == "strong", ]

# Convert TFBS & CRM coordinates  to GRanges
tfbs_gr <- makeGRangesFromDataFrame(
  mb_strong,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

crm_gr <- makeGRangesFromDataFrame(
  selfTargetingCRMs,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

# Harmonize chromosome names between TFBS and CRM objects
seqlevels(tfbs_gr) <- sub("^chr", "", seqlevels(tfbs_gr))
seqlevels(crm_gr) <- sub("^chr-?", "", seqlevels(crm_gr))

# Find overlaps between affected TFBS and CRMs
tfbs_crm_hits <- findOverlaps(
  tfbs_gr,
  crm_gr,
  ignore.strand = TRUE
)

cat("Strong TFBS predictions overlapping prioritized CRMs:", length(tfbs_crm_hits))

# Create overlap table
gr_to_table <- function(gr, prefix) {
  as.data.frame(gr) %>%
    dplyr::rename(
      "{prefix}_chr" := seqnames,
      "{prefix}_start" := start,
      "{prefix}_end" := end,
      "{prefix}_width" := width,
      "{prefix}_strand" := strand
    )
}

tfbs_table <- gr_to_table(tfbs_gr[queryHits(tfbs_crm_hits)], "tfbs")
crm_table <- gr_to_table(crm_gr[subjectHits(tfbs_crm_hits)], "crm")

tfbs_crm_overlap <- bind_cols(tfbs_table, crm_table)

tfbs_crm_overlap$motifPos <- sapply(tfbs_crm_overlap$motifPos,
  function(x) paste(x, collapse = ", "))

write.csv(tfbs_crm_overlap, "results/tables/tfbs_crm_overlap.csv", row.names = FALSE)

cat("Strong TFBS predictions overlapping prioritized CRMs:", length(tfbs_crm_hits))
cat("Distinct affected TFs overlapping CRMs:", n_distinct(tfbs_crm_overlap$geneSymbol))
cat("Distinct CRMs containing affected TFBS:", n_distinct(tfbs_crm_overlap$crm_name))

# 6.4. Summary table ----
crm_tfbs_overlap_summary <- tfbs_crm_overlap %>%
  group_by(crm_name) %>%
  summarise(
    crm_chr = crm_chr[1],
    crm_start = crm_start[1],
    crm_end = crm_end[1],
    TFBS_CRM_overlaps = n(),
    affected_SNPs = n_distinct(SNP_id),
    affected_TFs = n_distinct(geneSymbol),
    gains_or_strengthening = sum(TFBS_effect == "gain_or_strengthening"),
    losses_or_weakening = sum(TFBS_effect == "loss_or_weakening"),
    TFs = paste(geneSymbol, collapse = ", "),
    .groups = "drop"
  ) %>%
  left_join(selfTargetingCRMs %>% select(crm_name, common_genes), by = "crm_name"
  ) %>%
  arrange(desc(TFBS_CRM_overlaps), desc(affected_TFs))

write.csv(crm_tfbs_overlap_summary,"results/tables/crm_tfbs_overlap_summary.csv", row.names = FALSE)
