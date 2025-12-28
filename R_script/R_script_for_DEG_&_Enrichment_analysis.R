### ===============================
### Package setup
### ===============================

# Increase timeout for large Bioconductor packages
options(timeout = 600)

# Check / install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
} else {
  message("BiocManager is already installed")
}

# Check / install DESeq2
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
} else {
  message("DESeq2 is already installed")
}
library(DESeq2)

# Check / install clusterProfiler
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
} else {
  message("clusterProfiler is already installed")
}
library(clusterProfiler)

# Check / install mouse annotation database
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Mm.eg.db", ask = FALSE)
} else {
  message("org.Mm.eg.db is already installed")
}
library(org.Mm.eg.db)

library(ggplot2)

### ===============================
### Load counts and metadata
### ===============================

counts <- read.table(
  "counts_matrix.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

coldata <- read.table(
  "metadata.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# Ensure sample alignment
counts <- counts[, rownames(coldata)]

### ===============================
### DESeq2 analysis
### ===============================

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "lactation", "virgin"))
res <- res[order(res$padj), ]

# Save full results
write.csv(as.data.frame(res), "deseq2_results.csv")

# Separate NA and non-NA genes
res_empty <- res[is.na(res$padj), ]
res_clean <- res[!is.na(res$padj), ]

write.csv(as.data.frame(res_empty), "deseq2_empty_results.csv")
write.csv(as.data.frame(res_clean), "deseq2_clean_results.csv")

# Out of the 18,418 genes retained after filtering from the original 27,129 genes,
# five genes were assigned NA adjusted p-values by DESeq2 due to sparse or unstable read counts.

res["215866", ]
res["18948", ]
res["19259", ]
res["212753", ]
res["382231", ]

counts(dds)["215866", ]
counts(dds)["18948", ]
counts(dds)["19259", ]
counts(dds)["212753", ]
counts(dds)["382231", ]

dim(res_clean)

### ===============================
### Volcano plot
### ===============================

res_clean$significance <- "Not Significant"
res_clean$significance[res_clean$padj < 0.05 & res_clean$log2FoldChange > 1] <- "Up"
res_clean$significance[res_clean$padj < 0.05 & res_clean$log2FoldChange < -1] <- "Down"

ggplot(res_clean, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot: Lactation vs Virgin",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_minimal()

### ===============================
### Pathway enrichment (DESeq2 genes)
### ===============================

res <- read.csv("deseq2_clean_results.csv", row.names = 1)

res$significance <- "Not Significant"
res$significance[res$padj < 0.05 & res$log2FoldChange > 1] <- "Up"
res$significance[res$padj < 0.05 & res$log2FoldChange < -1] <- "Down"

up_genes <- rownames(res)[res$padj < 0.05 & res$log2FoldChange > 1]
down_genes <- rownames(res)[res$padj < 0.05 & res$log2FoldChange < -1]

ego_up <- enrichGO(
  gene          = up_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dotplot(ego_up, showCategory = 15) +
  ggtitle("GO Biological Processes (Upregulated in Lactation)")

ekegg_up <- enrichKEGG(
  gene          = up_genes,
  organism      = "mmu",
  pAdjustMethod = "BH"
)

dotplot(ekegg_up, showCategory = 15) +
  ggtitle("KEGG Pathways (Upregulated in Lactation)")

### ===============================
### Export DESeq2 size factors
### ===============================

sf <- sizeFactors(dds)

sf_df <- data.frame(
  sample = names(sf),
  size_factor = as.numeric(sf)
)

write.table(
  sf_df,
  file = "size_factors.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

### ===============================
### Pathway enrichment (Generative DE genes)
### ===============================

# From here, pathway enrichment analysis of up- & down-regulated genes
# derived from the generative NB-VAE model begins.

gen_up <- read.csv("gen_up_entrez.csv", header = FALSE)[, 1]
gen_down <- read.csv("gen_down_entrez.csv", header = FALSE)[, 1]

# A relaxed cutoff is used for generative DE to account for shrinkage and uncertainty
ego_gen_up <- enrichGO(
  gene          = gen_up,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.2
)

ekegg_gen_up <- enrichKEGG(
  gene          = gen_up,
  organism      = "mmu",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.2
)

dotplot(ego_gen_up, showCategory = 15)
dotplot(ekegg_gen_up, showCategory = 15)