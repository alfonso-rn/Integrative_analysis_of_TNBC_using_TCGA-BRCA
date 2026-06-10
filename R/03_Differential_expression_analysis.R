# 4. Differential expression analysis: TNBC vs non-TNBC ----
#
# Differential expression analysis is performed between clinically defined TNBC
# and non-TNBC primary tumor samples.
#
# DEGs are identified using FDR < 0.01 and |log2FC| >= 1, corresponding to
# statistically significant genes with at least a two-fold expression change.
#
# Since TNBC is provided as the second condition, positive logFC values are
# interpreted as higher expression in TNBC compared with non-TNBC.
#
# The resulting DEGs are summarized by regulation direction and gene type.
# In addition, an unfiltered DEA is performed to retain all tested features
# for visualization in the volcano plot.

library(TCGAbiolinks)
library(ggplot2)

dataFilt <- readRDS("data/processed/dataFilt.rds")
samplesTNBC <- readRDS("data/processed/samplesTNBC.rds")
samplesNonTNBC <- readRDS("data/processed/samplesNonTNBC.rds")

# 4.1. Differential expression analysis (DEA) ----

dataDEA <- TCGAanalyze_DEA(
  mat1 = dataFilt[, samplesNonTNBC],
  mat2 = dataFilt[, samplesTNBC],
  Cond1type = "Non_TNBC",
  Cond2type = "TNBC",
  fdr.cut = 0.01, 
  logFC.cut = 1,
  method = "glmLRT"
)

 cat("Total DEGs:", nrow(dataDEA), "\n",
    "Upregulated DEGs in TNBC:", sum(dataDEA$logFC > 0), "\n",
    "Downregulated DEGs in TNBC:", sum(dataDEA$logFC < 0))

# 4.2. DEGs table with expression values ----

dataDEGsLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEA,
  typeCond1 = "Non_TNBC",
  typeCond2 = "TNBC",
  TableCond1 = dataFilt[, samplesNonTNBC],
  TableCond2 = dataFilt[, samplesTNBC]
)

# Add expression TNBC status column
dataDEA$expression_status <- ifelse(
  dataDEA$logFC > 0,
  "Up regulated in TNBC",
  "Down regulated in TNBC"
)

dataDEGsLevel$expression_status <- ifelse(
  dataDEGsLevel$logFC > 0,
  "Up regulated in TNBC",
  "Down regulated in TNBC"
)

saveRDS(dataDEA, "data/processed/dataDEA.rds")
saveRDS(dataDEGsLevel, "data/processed/dataDEGsLevel.rds")

# 4.3. Volcano plot representation ----

globalDEA <- TCGAanalyze_DEA(
  mat1 = dataFilt[, samplesNonTNBC],
  mat2 = dataFilt[, samplesTNBC],
  Cond1type = "Non_TNBC",
  Cond2type = "TNBC",
  method = "glmLRT"
)

TCGAVisualize_volcano(
  x = globalDEA$logFC,
  y = globalDEA$FDR,
  filename = "results/figures/volcano_TNBC_vs_NonTNBC.png",
  x.cut = 1,
  y.cut = 0.01,
  names = globalDEA$gene_name,
  show.names = "significant",
  color = c("grey", "red", "blue"),
  names.size = 2,
  xlab = "Gene expression fold change (Log2)",
  legend = "State",
  title = "Differential expression analysis: TNBC vs Non-TNBC",
  width = 10
)

# 4.4. Horizontal bar plot representation ----

barplot_gene_type <- ggplot(
  data = dataDEA,
  mapping = aes(
    x = factor(
      ifelse(logFC > 0, "Upregulated in TNBC", "Downregulated in TNBC"),
      levels = c("Downregulated in TNBC", "Upregulated in TNBC")
    ),
    fill = gene_type
  )
) +
  geom_bar(position = "fill") +
  coord_flip() +
  scale_y_continuous(
    labels = function(x) paste0(round(x * 100), "%")
  ) +
  labs(
    title = "Distribution of gene types among DEGs",
    subtitle = "TNBC vs non-TNBC",
    x = NULL,
    y = "Percentage of DEGs",
    fill = "Gene type"
  ) +
  theme_classic()

ggsave(filename = "results/figures/barplot_DEGs_type.png",
       plot = barplot_gene_type, width = 10, height = 6, dpi = 300 )

# 4.5. Gene type distribution summary table ----

gene_type_counts <- table(dataDEA$gene_type, dataDEA$expression_status)

gene_type_percent <- round(prop.table(gene_type_counts, margin = 2) * 100, 2)

gene_type_percent <- data.frame(
  gene_type = rownames(gene_type_percent),
  Down_regulated_TNBC = gene_type_percent[, "Down regulated in TNBC"],
  Up_regulated_TNBC = gene_type_percent[, "Up regulated in TNBC"],
  row.names = NULL
)

write.csv(gene_type_percent, "results/tables/DEGs_type_distribution.csv", 
          row.names = FALSE)
