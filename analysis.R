# ==========================================
# Comparative Transcriptomic Analysis
# of Breast Cancer Using R
#
# Author: Diana Goldy
# Dataset: GSE3193
# Platform: GPL180
# ==========================================

# 1. Load required packages

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)

# 2. Download breast cancer dataset from GEO

gse <- getGEO("GSE3193", GSEMatrix = TRUE)

# Check available platforms
length(gse)

# 3. Extract expression data and sample metadata

eset <- gse[[1]]

pd <- pData(eset)

expr_matrix <- exprs(eset)

# Check dimensions
dim(expr_matrix)

# 4. Classify normal and tumour samples

normal_samples <- c(
  "GSM1887",
  "GSM1902",
  "GSM1906",
  "GSM73720"
)

group <- ifelse(
  colnames(expr_matrix) %in% normal_samples,
  "Normal",
  "Tumor"
)

group <- factor(group, levels = c("Normal", "Tumor"))

# Check sample distribution
table(group)

# 5. Filter probes based on data availability

normal_present <- rowSums(
  !is.na(expr_matrix[, group == "Normal"])
)

tumor_present <- rowSums(
  !is.na(expr_matrix[, group == "Tumor"])
)

keep <- normal_present >= 3 & tumor_present >= 20

expr_filtered <- expr_matrix[keep, , drop = FALSE]

# Check filtered expression matrix
dim(expr_filtered)

# 6. Retrieve platform annotation

gpl <- getGEO("GPL180")

annotation_data <- Table(gpl)[, c(1, 4, 5)]

colnames(annotation_data) <- c(
  "ID",
  "Gene_Symbol",
  "Gene_Name"
)

# Match annotation with filtered expression data

annotation_filtered <- annotation_data[
  match(rownames(expr_filtered), as.character(annotation_data$ID)),
]

# Check missing gene symbols
sum(is.na(annotation_filtered$Gene_Symbol))

# 7. Select samples with CR1 reference channel

keep_cr1 <- pd$source_name_ch1 == "CR1 reference"

expr_cr1 <- expr_filtered[, keep_cr1, drop = FALSE]

group_cr1 <- droplevels(group[keep_cr1])

# Check sample distribution
table(group_cr1)

# Check expression matrix dimensions
dim(expr_cr1)

# 8. Create design matrix

design_cr1 <- model.matrix(~0 + group_cr1)

colnames(design_cr1) <- levels(group_cr1)

# View design matrix
design_cr1

# 9. Define tumour versus normal comparison

contrast_cr1 <- makeContrasts(
  Tumor - Normal,
  levels = design_cr1
)

# 10. Perform differential expression analysis

fit <- lmFit(expr_cr1, design_cr1)

fit2 <- contrasts.fit(fit, contrast_cr1)

fit2 <- eBayes(fit2)

# 11. Extract differential expression results

results <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

# Count significant probes
sum(results$adj.P.Val < 0.05, na.rm = TRUE)

# 12. Identify significant probes

sig_results <- results[
  results$adj.P.Val < 0.05,
]

# Check number of significant probes
dim(sig_results)

# 13. Annotate significant probes with gene information

sig_results$ID <- rownames(sig_results)

sig_results_annotated <- merge(
  sig_results,
  annotation_filtered,
  by = "ID",
  all.x = TRUE
)

# View first five annotated results
head(
  sig_results_annotated[
    , c("ID", "logFC", "Gene_Symbol")
  ],
  5
)

# 14. Export significant gene results

write.csv(
  sig_results_annotated,
  "significant_genes_annotated.csv",
  row.names = FALSE
)

# 15. Generate volcano plot

results$Significance <- ifelse(
  results$adj.P.Val < 0.05,
  "Significant",
  "Not Significant"
)

volcano_plot <- ggplot(
  results,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Significance
  )
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c(
      "Significant" = "red",
      "Not Significant" = "grey"
    )
  ) +
  labs(
    title = "Volcano Plot: Breast Tumour vs Normal",
    x = "Log Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()

# Display volcano plot
print(volcano_plot)

# Save volcano plot
ggsave(
  "volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6
)

# 16. Generate heatmap of top 30 significant probes

# Sort significant probes by adjusted p-value
top_results <- sig_results[
  order(sig_results$adj.P.Val),
]

top_genes <- rownames(top_results)[1:30]

# Extract expression values
heatmap_data <- expr_cr1[
  top_genes,
  ,
  drop = FALSE
]

# Replace missing values with the row mean
for (i in seq_len(nrow(heatmap_data))) {
  heatmap_data[i, is.na(heatmap_data[i, ])] <-
    mean(heatmap_data[i, ], na.rm = TRUE)
}

# Save heatmap
png(
  "heatmap_top30.png",
  width = 1200,
  height = 1000,
  res = 150
)

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 30 Differentially Expressed Probes"
)

dev.off()

# 17. Prepare data for PCA analysis

# Select probes with the highest variability
gene_variance <- apply(
  expr_cr1,
  1,
  var,
  na.rm = TRUE
)

top_variable <- order(
  gene_variance,
  decreasing = TRUE
)[1:500]

pca_matrix <- expr_cr1[
  top_variable,
  ,
  drop = FALSE
]

# Replace missing values with row means
for (i in seq_len(nrow(pca_matrix))) {
  pca_matrix[i, is.na(pca_matrix[i, ])] <-
    mean(pca_matrix[i, ], na.rm = TRUE)
}

# 18. Perform PCA

pca_data <- prcomp(
  t(pca_matrix),
  scale. = TRUE
)

pca_df <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  Group = group_cr1
)

# Generate PCA plot
pca_plot <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Group)
) +
  geom_point(size = 4) +
  labs(
    title = "PCA: Breast Tumour vs Normal",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal()

print(pca_plot)

# Save PCA plot
ggsave(
  "PCA_plot.png",
  plot = pca_plot,
  width = 8,
  height = 6
)

# 19. Record R session information

sessionInfo()