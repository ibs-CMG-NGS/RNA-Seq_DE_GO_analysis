# 파일 경로: src/analysis/02a_generate_qc_plots.R
# 목적: DESeq2 Global QC 플롯 생성 (모든 샘플 대상)
# 사용법: Rscript 02a_generate_qc_plots.R --config config.yml --output_dir output/qc_plots

# --- 1. Setup: Load libraries and parse arguments ---
suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(optparse)
  library(ggplot2)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
})

# Argument parsing
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yml", 
              help = "Path to the config YAML file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", 
              help = "Directory to save QC plots", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Output directory must be supplied.", call. = FALSE)
}

# Load config
config <- yaml.load_file(opt$config)

# Create output directory if it doesn't exist
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== DESeq2 Global QC Plots Generation ===\n")
cat(paste0("Config file: ", opt$config, "\n"))
cat(paste0("Output directory: ", opt$output_dir, "\n\n"))

# --- 2. Load data and create DESeq2 object ---
cat("Loading data and creating DESeq2 object...\n")

source(here("src", "utils", "load_data.R"))
dds <- create_de_object(config_path = opt$config)

# Run DESeq2 analysis if not already done
if (!"dispersion" %in% names(mcols(dds))) {
  cat("Running DESeq2 analysis...\n")
  dds <- DESeq(dds)
}

# Get QC plot configuration
qc_config <- config$qc_plots
if (is.null(qc_config)) {
  # Default values if qc_plots section doesn't exist
  qc_config <- list(
    sample_distance_method = "euclidean",
    top_genes_n = 30,
    heatmap_cluster_rows = TRUE,
    heatmap_cluster_cols = TRUE
  )
}

# Get group variable for coloring
intgroup <- config$de_analysis$group_variable
meta_data <- as.data.frame(colData(dds))

# Get advanced options
adv_opts <- config$de_analysis$advanced_options
blind_mode <- if(!is.null(adv_opts$vst_blind)) adv_opts$vst_blind else FALSE

# --- 3. Variance Stabilizing Transformation ---
cat("\nApplying variance stabilizing transformation (VST)...\n")
vsd <- vst(dds, blind = blind_mode)

# --- 4. Plot 1: Sample Distance Heatmap ---
cat("Generating sample distance heatmap...\n")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL

# Choose color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Save heatmap without annotation to avoid compatibility issues
output_file <- file.path(opt$output_dir, "sample_distance_heatmap.png")
png(output_file, width = 10, height = 8, units = "in", res = 300, bg = "white")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-Sample Distance Heatmap")
dev.off()
cat(paste0("  Saved: ", output_file, "\n"))# --- 5. Plot 2: Dispersion Plot ---
cat("Generating dispersion plot...\n")

output_file <- file.path(opt$output_dir, "dispersion_plot.png")
png(output_file, width = 10, height = 8, units = "in", res = 300, bg = "white")
plotDispEsts(dds, 
             main = "Dispersion Estimates",
             ylab = "Dispersion",
             xlab = "Mean of normalized counts")
dev.off()
cat(paste0("  Saved: ", output_file, "\n"))

# --- 6. Plot 3: PCA Plot (enhanced version) ---
cat("Generating PCA plot...\n")

pca_ntop <- if(!is.null(adv_opts$pca_ntop)) adv_opts$pca_ntop else 500
pca_data <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE, ntop = pca_ntop)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create enhanced PCA plot
p <- ggplot(pca_data, aes_string(x = "PC1", y = "PC2", color = intgroup)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - Sample Clustering") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  coord_fixed()

output_file <- file.path(opt$output_dir, "pca_plot.png")
ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
cat(paste0("  Saved: ", output_file, "\n"))

# --- 7. Plot 4: Scree Plot (PCA variance explained) ---
cat("Generating PCA scree plot...\n")

# Calculate PCA for scree plot
pca_res <- prcomp(t(assay(vsd)))
variance_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

scree_data <- data.frame(
  PC = paste0("PC", 1:min(10, length(variance_explained))),
  Variance = variance_explained[1:min(10, length(variance_explained))]
)
scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)

p_scree <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Variance)), 
            vjust = -0.5, size = 3.5) +
  ylab("Variance Explained (%)") +
  xlab("Principal Component") +
  ggtitle("PCA Scree Plot - Variance Explained") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

output_file <- file.path(opt$output_dir, "pca_scree_plot.png")
ggsave(output_file, plot = p_scree, width = 10, height = 6, dpi = 300, bg = "white")
cat(paste0("  Saved: ", output_file, "\n"))

# --- 8. Plot 5: Count Distribution (boxplot) ---
cat("Generating count distribution boxplot...\n")

# Get normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Prepare data for plotting
count_data_long <- reshape2::melt(log2(norm_counts + 1))
colnames(count_data_long) <- c("Gene", "Sample", "Log2Count")
count_data_long$Group <- meta_data[count_data_long$Sample, intgroup]

p_boxplot <- ggplot(count_data_long, aes(x = Sample, y = Log2Count, fill = Group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  ylab("Log2(Normalized Count + 1)") +
  xlab("Sample") +
  ggtitle("Distribution of Normalized Counts Across Samples") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

output_file <- file.path(opt$output_dir, "count_distribution_boxplot.png")
ggsave(output_file, plot = p_boxplot, width = 12, height = 8, dpi = 300, bg = "white")
cat(paste0("  Saved: ", output_file, "\n"))

# --- 9. Summary ---
cat("\n=== QC Plots Generation Complete ===\n")
cat("Generated plots:\n")
cat("  1. sample_distance_heatmap.png - Sample similarity matrix\n")
cat("  2. dispersion_plot.png - Gene dispersion estimates\n")
cat("  3. pca_plot.png - Principal component analysis\n")
cat("  4. pca_scree_plot.png - PCA variance explained\n")
cat("  5. count_distribution_boxplot.png - Normalized count distribution\n")
cat(paste0("\nAll plots saved to: ", opt$output_dir, "\n\n"))
