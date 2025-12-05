# 파일 경로: src/analysis/02b_generate_pairwise_qc_plots.R
# 목적: DESeq2 Pairwise QC 플롯 생성 (각 비교군별)
# 사용법: Rscript 02b_generate_pairwise_qc_plots.R --config config.yml --comparison "H2O2_vs_Control" --de_results output/pairwise/H2O2_vs_Control/final_de_results.csv --output_dir output/pairwise/H2O2_vs_Control/qc_plots

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
  make_option(c("--comparison"), type = "character", 
              help = "Comparison name (e.g., 'H2O2_vs_Control')", metavar = "character"),
  make_option(c("-i", "--de_results"), type = "character", 
              help = "Path to DE results CSV file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", 
              help = "Directory to save QC plots", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$comparison) || is.null(opt$de_results) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("comparison, de_results, and output_dir arguments must be supplied.", call. = FALSE)
}

# Load config
config <- yaml.load_file(opt$config)

# Create output directory if it doesn't exist
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== DESeq2 Pairwise QC Plots Generation ===\n")
cat(paste0("Comparison: ", opt$comparison, "\n"))
cat(paste0("DE results: ", opt$de_results, "\n"))
cat(paste0("Output directory: ", opt$output_dir, "\n\n"))

# --- 2. Load DE results ---
cat("Loading DE results...\n")
de_results <- read.csv(opt$de_results, row.names = 1)

# Get QC plot configuration
qc_config <- config$qc_plots
if (is.null(qc_config)) {
  # Default values
  qc_config <- list(
    top_genes_n = 30,
    ma_plot_ylim = c(-5, 5),
    heatmap_cluster_rows = TRUE,
    heatmap_cluster_cols = TRUE
  )
}

top_genes_n <- if(!is.null(qc_config$top_genes_n)) qc_config$top_genes_n else 30
ma_ylim <- if(!is.null(qc_config$ma_plot_ylim)) qc_config$ma_plot_ylim else c(-5, 5)

# --- 3. Load DESeq2 object for additional plots ---
cat("Loading DESeq2 object...\n")
source(here("src", "utils", "load_data.R"))
dds <- create_de_object(config_path = opt$config)

# Run DESeq2 if not already done
if (!"dispersion" %in% names(mcols(dds))) {
  cat("Running DESeq2 analysis...\n")
  dds <- DESeq(dds)
}

# Get VST-transformed data
adv_opts <- config$de_analysis$advanced_options
blind_mode <- if(!is.null(adv_opts$vst_blind)) adv_opts$vst_blind else FALSE
vsd <- vst(dds, blind = blind_mode)

# --- 4. Plot 1: MA Plot ---
cat("Generating MA plot...\n")

# Prepare data for MA plot
ma_data <- data.frame(
  baseMean = de_results$baseMean,
  log2FoldChange = de_results$log2FoldChange,
  padj = de_results$padj
)
ma_data$significant <- ifelse(is.na(ma_data$padj), FALSE, ma_data$padj < 0.05)

# Create MA plot
p_ma <- ggplot(ma_data, aes(x = log10(baseMean), y = log2FoldChange, color = significant)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("gray50", "red"), 
                     name = "Significant",
                     labels = c("No", "Yes (padj < 0.05)")) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
  ylim(ma_ylim) +
  xlab("Log10(Mean Expression)") +
  ylab("Log2(Fold Change)") +
  ggtitle(paste0("MA Plot - ", opt$comparison)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

output_file <- file.path(opt$output_dir, "ma_plot.png")
ggsave(output_file, plot = p_ma, width = 10, height = 8, dpi = 300, bg = "white")
cat(paste0("  Saved: ", output_file, "\n"))

# --- 5. Plot 2: P-value Histogram ---
cat("Generating p-value histogram...\n")

pval_data <- data.frame(pvalue = de_results$pvalue)
pval_data <- pval_data[!is.na(pval_data$pvalue), , drop = FALSE]

p_hist <- ggplot(pval_data, aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  xlab("P-value") +
  ylab("Frequency") +
  ggtitle(paste0("P-value Distribution - ", opt$comparison)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

output_file <- file.path(opt$output_dir, "pvalue_histogram.png")
ggsave(output_file, plot = p_hist, width = 10, height = 6, dpi = 300, bg = "white")
cat(paste0("  Saved: ", output_file, "\n"))

# --- 6. Plot 3: Adjusted P-value Histogram ---
cat("Generating adjusted p-value histogram...\n")

padj_data <- data.frame(padj = de_results$padj)
padj_data <- padj_data[!is.na(padj_data$padj), , drop = FALSE]

p_padj_hist <- ggplot(padj_data, aes(x = padj)) +
  geom_histogram(bins = 50, fill = "coral", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = 0.05, y = Inf, label = "padj = 0.05", 
           vjust = 1.5, hjust = -0.1, color = "red", size = 4) +
  xlab("Adjusted P-value") +
  ylab("Frequency") +
  ggtitle(paste0("Adjusted P-value Distribution - ", opt$comparison)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

output_file <- file.path(opt$output_dir, "padj_histogram.png")
ggsave(output_file, plot = p_padj_hist, width = 10, height = 6, dpi = 300, bg = "white")
cat(paste0("  Saved: ", output_file, "\n"))

# --- 7. Plot 4: Top DE Genes Heatmap ---
cat(paste0("Generating heatmap of top ", top_genes_n, " DE genes...\n"))

# Sort by adjusted p-value and select top N
de_results_sorted <- de_results[order(de_results$padj), ]
top_genes <- rownames(de_results_sorted)[1:min(top_genes_n, nrow(de_results_sorted))]

# Get VST counts for top genes
top_genes_vst <- assay(vsd)[top_genes, ]

# Z-score normalization for better visualization
top_genes_scaled <- t(scale(t(top_genes_vst)))

# Save heatmap - simplified without annotations to avoid compatibility issues
output_file <- file.path(opt$output_dir, "top_genes_heatmap.png")
png(output_file, width = 12, height = 10, units = "in", res = 300, bg = "white")

pheatmap(top_genes_scaled,
         cluster_rows = if(!is.null(qc_config$heatmap_cluster_rows)) qc_config$heatmap_cluster_rows else TRUE,
         cluster_cols = if(!is.null(qc_config$heatmap_cluster_cols)) qc_config$heatmap_cluster_cols else TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = paste0("Top ", length(top_genes), " DE Genes - ", opt$comparison),
         fontsize_row = 8,
         fontsize_col = 10)

dev.off()
cat(paste0("  Saved: ", output_file, "\n"))

# --- 8. Plot 5: Log2FC Distribution ---
cat("Generating log2 fold change distribution...\n")

fc_data <- data.frame(
  log2FC = de_results$log2FoldChange,
  significant = ifelse(is.na(de_results$padj), FALSE, de_results$padj < 0.05)
)
fc_data <- fc_data[!is.na(fc_data$log2FC), ]

p_fc <- ggplot(fc_data, aes(x = log2FC, fill = significant)) +
  geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("gray50", "red"),
                    name = "Significant",
                    labels = c("No", "Yes (padj < 0.05)")) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  xlab("Log2(Fold Change)") +
  ylab("Frequency") +
  ggtitle(paste0("Log2 Fold Change Distribution - ", opt$comparison)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

output_file <- file.path(opt$output_dir, "log2fc_distribution.png")
ggsave(output_file, plot = p_fc, width = 10, height = 6, dpi = 300)
cat(paste0("  Saved: ", output_file, "\n"))

# --- 9. Plot 6: Volcano plot-style effect size vs significance ---
cat("Generating effect size vs significance plot...\n")

volcano_data <- data.frame(
  log2FC = abs(de_results$log2FoldChange),
  neglog10padj = -log10(de_results$padj),
  significant = ifelse(is.na(de_results$padj), FALSE, de_results$padj < 0.05)
)
volcano_data <- volcano_data[!is.na(volcano_data$log2FC) & !is.na(volcano_data$neglog10padj), ]

p_effect <- ggplot(volcano_data, aes(x = log2FC, y = neglog10padj, color = significant)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("gray50", "red"),
                     name = "Significant",
                     labels = c("No", "Yes (padj < 0.05)")) +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
  xlab("|Log2(Fold Change)|") +
  ylab("-Log10(Adjusted P-value)") +
  ggtitle(paste0("Effect Size vs Significance - ", opt$comparison)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

output_file <- file.path(opt$output_dir, "effect_size_vs_significance.png")
ggsave(output_file, plot = p_effect, width = 10, height = 8, dpi = 300)
cat(paste0("  Saved: ", output_file, "\n"))

# --- 10. Summary Statistics ---
cat("\n=== Summary Statistics ===\n")
total_genes <- nrow(de_results)
sig_genes <- sum(de_results$padj < 0.05, na.rm = TRUE)
up_genes <- sum(de_results$padj < 0.05 & de_results$log2FoldChange > 0, na.rm = TRUE)
down_genes <- sum(de_results$padj < 0.05 & de_results$log2FoldChange < 0, na.rm = TRUE)

cat(paste0("Total genes tested: ", total_genes, "\n"))
cat(paste0("Significant genes (padj < 0.05): ", sig_genes, " (", 
           round(sig_genes/total_genes*100, 2), "%)\n"))
cat(paste0("  Up-regulated: ", up_genes, "\n"))
cat(paste0("  Down-regulated: ", down_genes, "\n"))

# --- 11. Summary ---
cat("\n=== Pairwise QC Plots Generation Complete ===\n")
cat("Generated plots:\n")
cat("  1. ma_plot.png - Mean expression vs log2FC\n")
cat("  2. pvalue_histogram.png - Raw p-value distribution\n")
cat("  3. padj_histogram.png - Adjusted p-value distribution\n")
cat("  4. top_genes_heatmap.png - Heatmap of top DE genes\n")
cat("  5. log2fc_distribution.png - Distribution of fold changes\n")
cat("  6. effect_size_vs_significance.png - Effect size vs significance\n")
cat(paste0("\nAll plots saved to: ", opt$output_dir, "\n\n"))
