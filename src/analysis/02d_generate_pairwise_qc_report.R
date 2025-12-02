# 파일 경로: src/analysis/02d_generate_pairwise_qc_report.R
# 목적: Pairwise QC 플롯들을 모아 HTML 리포트 생성
# 사용법: Rscript 02d_generate_pairwise_qc_report.R --config config.yml --comparison "H2O2_vs_Control" --de_results output/pairwise/H2O2_vs_Control/final_de_results.csv --qc_plots_dir output/pairwise/H2O2_vs_Control/qc_plots --output_file output/pairwise/H2O2_vs_Control/qc_plots/pairwise_qc_report.html

# --- 1. Setup: Load libraries and parse arguments ---
suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(optparse)
})

# Argument parsing
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yml", 
              help = "Path to the config YAML file", metavar = "character"),
  make_option(c("--comparison"), type = "character", 
              help = "Comparison name (e.g., 'H2O2_vs_Control')", metavar = "character"),
  make_option(c("-i", "--de_results"), type = "character", 
              help = "Path to DE results CSV file", metavar = "character"),
  make_option(c("-q", "--qc_plots_dir"), type = "character", 
              help = "Directory containing QC plots", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", 
              help = "Output HTML file path", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$comparison) || is.null(opt$de_results) || is.null(opt$qc_plots_dir) || is.null(opt$output_file)) {
  print_help(opt_parser)
  stop("comparison, de_results, qc_plots_dir, and output_file arguments must be supplied.", call. = FALSE)
}

# Load config
config <- yaml.load_file(opt$config)

cat("\n=== Pairwise QC Report Generation ===\n")
cat(paste0("Comparison: ", opt$comparison, "\n"))
cat(paste0("QC plots directory: ", opt$qc_plots_dir, "\n"))
cat(paste0("Output HTML file: ", opt$output_file, "\n\n"))

# --- 2. Load DE results for summary ---
de_results <- read.csv(opt$de_results, row.names = 1)

# Calculate statistics
total_genes <- nrow(de_results)
sig_genes <- sum(de_results$padj < 0.05, na.rm = TRUE)
up_genes <- sum(de_results$padj < 0.05 & de_results$log2FoldChange > 0, na.rm = TRUE)
down_genes <- sum(de_results$padj < 0.05 & de_results$log2FoldChange < 0, na.rm = TRUE)

# Get top DE genes
de_sorted <- de_results[order(de_results$padj), ]
top_n <- min(20, nrow(de_sorted))
top_genes <- de_sorted[1:top_n, ]

# Parse comparison name
comp_parts <- strsplit(opt$comparison, "_vs_")[[1]]
compare_group <- comp_parts[1]
base_group <- comp_parts[2]

# --- 3. Define plot files ---
plot_files <- c(
  "ma_plot.png",
  "pvalue_histogram.png",
  "padj_histogram.png",
  "top_genes_heatmap.png",
  "log2fc_distribution.png",
  "effect_size_vs_significance.png"
)

plot_titles <- c(
  "MA Plot",
  "P-value Distribution",
  "Adjusted P-value Distribution",
  "Top DE Genes Heatmap",
  "Log2 Fold Change Distribution",
  "Effect Size vs Significance"
)

plot_descriptions <- c(
  "MA plot shows the relationship between mean expression (x-axis) and log2 fold change (y-axis). Red points indicate significant genes (padj < 0.05). The plot should be symmetric around y=0 with no expression-dependent bias.",
  "Distribution of raw p-values across all tested genes. A good analysis shows enrichment of small p-values (left side) and a relatively uniform distribution for non-significant genes. A flat distribution suggests low statistical power.",
  "Distribution of adjusted p-values after multiple testing correction. The vertical red line indicates the significance threshold (padj = 0.05). Most genes should have high adjusted p-values, with a clear separation of significant genes.",
  paste0("Heatmap showing normalized expression of the top ", top_n, " differentially expressed genes across all samples. Rows are genes, columns are samples. Colors represent z-scores of expression values. Samples should cluster by treatment group."),
  "Distribution of log2 fold changes for all genes. Red indicates significant genes (padj < 0.05). The distribution should be centered around 0, with significant genes showing larger fold changes in both directions.",
  "Scatter plot of absolute log2 fold change vs -log10(adjusted p-value). Points in the upper right have both large effect sizes and high significance. The horizontal line indicates padj = 0.05 threshold."
)

# --- 4. Generate HTML content ---
html_content <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pairwise QC Report - ', opt$comparison, '</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .header {
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .header h1 {
            margin: 0 0 10px 0;
            font-size: 2.5em;
        }
        .header p {
            margin: 5px 0;
            font-size: 1.1em;
            opacity: 0.9;
        }
        .comparison-badge {
            display: inline-block;
            background: rgba(255,255,255,0.2);
            padding: 8px 16px;
            border-radius: 20px;
            margin-top: 10px;
            font-weight: bold;
        }
        .summary {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .summary h2 {
            margin-top: 0;
            color: #f5576c;
            border-bottom: 2px solid #f5576c;
            padding-bottom: 10px;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }
        .summary-item {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #f5576c;
        }
        .summary-item.up {
            border-left-color: #e74c3c;
        }
        .summary-item.down {
            border-left-color: #3498db;
        }
        .summary-item .label {
            font-size: 0.9em;
            color: #6c757d;
            margin-bottom: 5px;
        }
        .summary-item .value {
            font-size: 1.5em;
            font-weight: bold;
            color: #212529;
        }
        .summary-item .percentage {
            font-size: 0.9em;
            color: #6c757d;
            margin-top: 5px;
        }
        .plot-section {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .plot-section h2 {
            margin-top: 0;
            color: #f5576c;
            border-bottom: 2px solid #f5576c;
            padding-bottom: 10px;
        }
        .plot-description {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin: 15px 0;
            font-size: 0.95em;
            color: #495057;
            border-left: 4px solid #f093fb;
        }
        .plot-container {
            text-align: center;
            margin: 20px 0;
        }
        .plot-container img {
            max-width: 100%;
            height: auto;
            border-radius: 5px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        .footer {
            text-align: center;
            margin-top: 40px;
            padding: 20px;
            color: #6c757d;
            font-size: 0.9em;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            font-size: 0.9em;
        }
        th, td {
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #dee2e6;
        }
        th {
            background-color: #f5576c;
            color: white;
            font-weight: 600;
            position: sticky;
            top: 0;
        }
        tr:hover {
            background-color: #f8f9fa;
        }
        .gene-table-container {
            max-height: 500px;
            overflow-y: auto;
            border: 1px solid #dee2e6;
            border-radius: 5px;
        }
        .up-regulated {
            color: #e74c3c;
            font-weight: bold;
        }
        .down-regulated {
            color: #3498db;
            font-weight: bold;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>📊 Pairwise QC Report</h1>
        <p>Differential Expression Analysis</p>
        <div class="comparison-badge">
            ', compare_group, ' vs ', base_group, '
        </div>
        <p style="margin-top: 15px;">Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
    </div>

    <div class="summary">
        <h2>📈 Analysis Summary</h2>
        <div class="summary-grid">
            <div class="summary-item">
                <div class="label">Total Genes Tested</div>
                <div class="value">', format(total_genes, big.mark=","), '</div>
            </div>
            <div class="summary-item">
                <div class="label">Significant Genes</div>
                <div class="value">', format(sig_genes, big.mark=","), '</div>
                <div class="percentage">', round(sig_genes/total_genes*100, 2), '% of total</div>
            </div>
            <div class="summary-item up">
                <div class="label">Up-regulated</div>
                <div class="value">', format(up_genes, big.mark=","), '</div>
                <div class="percentage">', round(up_genes/sig_genes*100, 1), '% of significant</div>
            </div>
            <div class="summary-item down">
                <div class="label">Down-regulated</div>
                <div class="value">', format(down_genes, big.mark=","), '</div>
                <div class="percentage">', round(down_genes/sig_genes*100, 1), '% of significant</div>
            </div>
        </div>

        <h3 style="margin-top: 25px; color: #495057;">Top ', top_n, ' Differentially Expressed Genes</h3>
        <div class="gene-table-container">
            <table>
                <thead>
                    <tr>
                        <th>Rank</th>
                        <th>Gene</th>
                        <th>Log2 Fold Change</th>
                        <th>Adjusted P-value</th>
                        <th>Regulation</th>
                    </tr>
                </thead>
                <tbody>
')

# Add top genes table
for (i in 1:top_n) {
  gene_name <- rownames(top_genes)[i]
  log2fc <- top_genes$log2FoldChange[i]
  padj <- top_genes$padj[i]
  regulation <- ifelse(log2fc > 0, "Up", "Down")
  regulation_class <- ifelse(log2fc > 0, "up-regulated", "down-regulated")
  
  html_content <- paste0(html_content, '
                    <tr>
                        <td>', i, '</td>
                        <td><strong>', gene_name, '</strong></td>
                        <td>', sprintf("%.2f", log2fc), '</td>
                        <td>', sprintf("%.2e", padj), '</td>
                        <td class="', regulation_class, '">', regulation, '</td>
                    </tr>
  ')
}

html_content <- paste0(html_content, '
                </tbody>
            </table>
        </div>
    </div>
')

# --- 5. Add plot sections ---
for (i in seq_along(plot_files)) {
  plot_path <- file.path(opt$qc_plots_dir, plot_files[i])
  
  # Convert to relative path for HTML
  rel_path <- plot_files[i]
  
  html_content <- paste0(html_content, '
    <div class="plot-section">
        <h2>', i, '. ', plot_titles[i], '</h2>
        <div class="plot-description">
            <strong>💡 Interpretation:</strong> ', plot_descriptions[i], '
        </div>
        <div class="plot-container">
            <img src="', rel_path, '" alt="', plot_titles[i], '">
        </div>
    </div>
  ')
}

# --- 6. Add footer ---
html_content <- paste0(html_content, '
    <div class="footer">
        <p>Generated by RNA-Seq DE Analysis Pipeline</p>
        <p>Configuration: ', opt$config, ' | Comparison: ', opt$comparison, '</p>
    </div>
</body>
</html>
')

# --- 7. Write HTML file ---
writeLines(html_content, opt$output_file)

cat("\n=== Pairwise QC Report Generation Complete ===\n")
cat(paste0("HTML report saved to: ", opt$output_file, "\n"))
cat(paste0("\nTo view the report, open it in a web browser:\n"))
cat(paste0("  file://", normalizePath(opt$output_file), "\n\n"))
