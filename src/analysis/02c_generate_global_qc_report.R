# 파일 경로: src/analysis/02c_generate_global_qc_report.R
# 목적: Global QC 플롯들을 모아 HTML 리포트 생성
# 사용법: Rscript 02c_generate_global_qc_report.R --config config.yml --qc_plots_dir output/qc_plots --output_file output/qc_plots/global_qc_report.html

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
  make_option(c("-q", "--qc_plots_dir"), type = "character", 
              help = "Directory containing QC plots", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", 
              help = "Output HTML file path", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$qc_plots_dir) || is.null(opt$output_file)) {
  print_help(opt_parser)
  stop("qc_plots_dir and output_file arguments must be supplied.", call. = FALSE)
}

# Load config
config <- yaml.load_file(opt$config)

cat("\n=== Global QC Report Generation ===\n")
cat(paste0("QC plots directory: ", opt$qc_plots_dir, "\n"))
cat(paste0("Output HTML file: ", opt$output_file, "\n\n"))

# --- 2. Load metadata for summary ---
meta_data <- read.csv(here(config$metadata_path), row.names = 1)
count_data <- read.csv(here(config$count_data_path), row.names = 1, check.names = FALSE)

# Get analysis info
group_var <- config$de_analysis$group_variable
groups <- unique(meta_data[[group_var]])
n_samples <- ncol(count_data)
n_genes <- nrow(count_data)

# --- 3. Define plot files ---
plot_files <- c(
  "sample_distance_heatmap.png",
  "dispersion_plot.png",
  "pca_plot.png",
  "pca_scree_plot.png",
  "count_distribution_boxplot.png"
)

plot_titles <- c(
  "Sample Distance Heatmap",
  "Dispersion Estimates",
  "PCA - Sample Clustering",
  "PCA Scree Plot",
  "Count Distribution Across Samples"
)

plot_descriptions <- c(
  "Hierarchical clustering of samples based on Euclidean distance. Samples from the same group should cluster together. Outliers may appear as samples that don't cluster with their group.",
  "Gene-wise dispersion estimates from DESeq2. Points should follow the fitted trend line. Genes with unusually high dispersion may indicate technical artifacts or biological heterogeneity.",
  "Principal Component Analysis showing the first two principal components. Samples should separate by treatment groups. PC1 and PC2 capture the largest sources of variation in the data.",
  "Variance explained by each principal component. Higher PCs explain progressively less variance. The first few PCs should capture most of the variation if the signal is strong.",
  "Distribution of normalized gene expression across all samples. Samples should have similar distributions after normalization. Large differences may indicate technical issues."
)

# --- 4. Generate HTML content ---
html_content <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Global QC Report - ', basename(config$output_dir), '</title>
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
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
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
        .summary {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .summary h2 {
            margin-top: 0;
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }
        .summary-item {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
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
        .plot-section {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .plot-section h2 {
            margin-top: 0;
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }
        .plot-description {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin: 15px 0;
            font-size: 0.95em;
            color: #495057;
            border-left: 4px solid #764ba2;
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
        }
        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #dee2e6;
        }
        th {
            background-color: #667eea;
            color: white;
            font-weight: 600;
        }
        tr:hover {
            background-color: #f8f9fa;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>🔬 Global QC Report</h1>
        <p>RNA-Seq Differential Expression Analysis</p>
        <p>Project: ', basename(config$output_dir), '</p>
        <p>Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
    </div>

    <div class="summary">
        <h2>📊 Dataset Summary</h2>
        <div class="summary-grid">
            <div class="summary-item">
                <div class="label">Total Samples</div>
                <div class="value">', n_samples, '</div>
            </div>
            <div class="summary-item">
                <div class="label">Total Genes</div>
                <div class="value">', format(n_genes, big.mark=","), '</div>
            </div>
            <div class="summary-item">
                <div class="label">Groups</div>
                <div class="value">', length(groups), '</div>
            </div>
            <div class="summary-item">
                <div class="label">Analysis Method</div>
                <div class="value">', config$de_analysis$method, '</div>
            </div>
        </div>

        <h3 style="margin-top: 25px; color: #495057;">Sample Groups</h3>
        <table>
            <thead>
                <tr>
                    <th>Group</th>
                    <th>Number of Samples</th>
                    <th>Sample IDs</th>
                </tr>
            </thead>
            <tbody>
')

# Add group information to table
for (grp in groups) {
  samples_in_group <- rownames(meta_data)[meta_data[[group_var]] == grp]
  html_content <- paste0(html_content, '
                <tr>
                    <td><strong>', grp, '</strong></td>
                    <td>', length(samples_in_group), '</td>
                    <td>', paste(samples_in_group, collapse=", "), '</td>
                </tr>
  ')
}

html_content <- paste0(html_content, '
            </tbody>
        </table>
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
        <p>Configuration: ', opt$config, '</p>
    </div>
</body>
</html>
')

# --- 7. Write HTML file ---
writeLines(html_content, opt$output_file)

cat("\n=== Global QC Report Generation Complete ===\n")
cat(paste0("HTML report saved to: ", opt$output_file, "\n"))
cat(paste0("\nTo view the report, open it in a web browser:\n"))
cat(paste0("  file://", normalizePath(opt$output_file), "\n\n"))
