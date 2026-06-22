# 7. Functional enrichment analysis of CRM-associated genes ----
#
# Functional enrichment analysis was performed on differentially expressed genes
# (DEGs) in TNBC using g:Profiler. Upregulated and downregulated genes were analysed
# separately to identify enriched Gene Ontology terms and biological pathways.
#
# Enriched terms containing genes associated with prioritized intragenic CRMs 
# were then identified to explore functional categories potentially linked to 
# regulatory regions affected by intronic mutations in TNBC.
#
# Enrichment sources:
# GO:BP = Biological Process / GO:MF = Molecular Function / GO:CC = Cellular Component
# REAC = Reactome pathways
# KEGG = Kyoto Encyclopedia of Genes and Genomes pathways
# WP = WikiPathways

library(dplyr)
library(gprofiler2)
library(ggplot2)

dataDEA <- readRDS("data/processed/dataDEA.rds")
top20_mutated_genes <- readRDS("data/processed/top20_mutated_genes.rds")
intragenicCRMs <- readRDS("data/processed/intragenicCRMs.rds")

# 7.1. Functional enrichment analysis ----

genes_up <- dataDEA %>%
  dplyr::filter(expression_status == "Up regulated in TNBC") %>%
  pull(gene_name)

genes_down <- dataDEA %>%
  dplyr::filter(expression_status == "Down regulated in TNBC") %>%
  pull(gene_name)

enrichment <- gost(
  query = list(Upregulated_TNBC = genes_up, Downregulated_TNBC = genes_down),
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG
              ", "WP"),
  correction_method = "fdr",
  evcodes = TRUE
)

saveRDS(enrichment, "data/processed/enrichment.rds")

# Results summary tables
enrichment_DEGs <- enrichment$result

enrichment_DEGs_terms <- enrichment_DEGs %>% count(query, source, name = "significant_terms")

write.csv(enrichment_DEGs_terms, "results/tables/enrichment_DEGs_terms.csv", row.names = FALSE)

enrichment_top_abundant <- enrichment_DEGs %>%
  arrange(desc(intersection_size)) %>%
  slice_head(n = 10) %>%
  mutate(parents = as.character(parents))

write.csv(enrichment_top_abundant, "results/tables/enrichment_top_abundant.csv", row.names = FALSE)

# Manhattan-like plot of enriched terms associated with top mutated DEGs
top_genes <- unique(unlist(top20_mutated_genes))

top_terms <- enrichment_DEGs %>%
  dplyr::filter(grepl(paste(top_genes, collapse = "|"), intersection)) %>%
  pull(term_id) %>%
  unique()

top_terms_significant <- enrichment_DEGs %>%
  dplyr::filter(term_id %in% top_terms) %>%
  arrange(p_value) %>%
  slice_head(n = 10)

DEGs_enrichmentplot1 <- gostplot(enrichment, capped = TRUE, interactive = FALSE)

DEGs_enrichmentplot2 <- publish_gostplot(DEGs_enrichmentplot1, 
                                  highlight_terms = top_terms_significant, 
                                  width = NA, height = NA, filename = NULL )

ggsave(filename = "results/figures/topmutatedDEGs_enrichmentplot.png", plot = DEGs_enrichmentplot2,
       width = 10, height = 7, dpi = 300, bg = "white")

# 7.2. Identification of enriched terms containing intragenic CRM genes ----

genesCRM <- intragenicCRMs$common_genes

enrichment_DEGs$CRM_genes_in_term <- NA_character_

for (i in seq_len(nrow(enrichment_DEGs))) {
  
  genes_intersection <- unlist(strsplit(enrichment_DEGs$intersection[i], ","))
  genes_intersection <- trimws(genes_intersection)
  
  CRM_genes <- intersect(genes_intersection, genesCRM)
  
  if (length(CRM_genes) > 0) {
    enrichment_DEGs$CRM_genes_in_term[i] <- paste(CRM_genes, collapse = "; ")
  }
}

enrichment_intragenicCRMs <- enrichment_DEGs[!is.na(enrichment_DEGs$CRM_genes_in_term), ]
enrichment_intragenicCRMs$parents <- as.character(enrichment_intragenicCRMs$parents)

write.csv(enrichment_intragenicCRMs, "results/tables/enrichment_intragenicCRMs.csv", row.names = FALSE)

# 7.3. Top 5 enriched terms associated with intragenic CRM genes ----

enrichment_intragenicCRMs_top <- enrichment_intragenicCRMs %>%
  mutate(
    CRM_gene_count = lengths(strsplit(CRM_genes_in_term, "; ")),
    minus_log10_p = -log10(p_value),
    term_label = paste(term_name, source, CRM_genes_in_term, sep = " | ")
  ) %>%
  arrange(query, source, p_value) %>%
  group_by(query, source) %>%
  slice_head(n = 5) %>%
  ungroup()

enrichment_intragenicCRMs_dotplot <- ggplot(
  enrichment_intragenicCRMs_top,
  aes(
    x = minus_log10_p,
    y = reorder(term_label, minus_log10_p),
    size = CRM_gene_count
  )
) +
  geom_point() +
  facet_wrap(~ query, scales = "free_y") +
  labs(
    x = "-log10(FDR-adjusted p-value)",
    y = NULL,
    size = "CRM genes\nin term"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold")
  )

ggsave("results/figures/enrichment_intragenicCRMs_dotplot.png", enrichment_intragenicCRMs_dotplot,
       width = 13, height = 9, dpi = 300, bg = "white")
