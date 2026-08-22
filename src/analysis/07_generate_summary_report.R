# 파일 경로: src/analysis/07_generate_summary_report.R
# 목적: DE-GO 파이프라인 완료 후 모든 pairwise 비교군 결과를 통합 HTML 요약 리포트로 생성
# 사용법:
#   Rscript 07_generate_summary_report.R \
#     --config configs/config_PROJECT.yml \
#     --output-dir output/PROJECT \
#     --output output/PROJECT/summary_report.html

# --- 1. Setup ---
suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
})

option_list <- list(
  make_option(c("-c", "--config"), type = "character",
              help = "Path to the config YAML file", metavar = "character"),
  make_option(c("-d", "--output-dir"), type = "character",
              help = "Pipeline output directory (e.g. output/PROJECT)", metavar = "character"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output HTML file path", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config) || is.null(opt$`output-dir`) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("--config, --output-dir and --output are required.", call. = FALSE)
}

config      <- yaml.load_file(opt$config)
output_dir  <- normalizePath(opt$`output-dir`, mustWork = TRUE)
output_file <- opt$output

cat("\n=== Summary Report Generation ===\n")
cat("Config:     ", opt$config, "\n")
cat("Output dir: ", output_dir, "\n")
cat("HTML:       ", output_file, "\n\n")

# --- 2. Derive pairs from config ---
get_pairs <- function(cfg) {
  pairs <- character(0)
  for (p in cfg$de_analysis$pairwise_comparisons) {
    pairs <- c(pairs, paste0(p[[1]], "_vs_", p[[2]]))
  }
  pairs
}

# R doesn't have %||% built-in before 4.4 — define it
`%||%` <- function(a, b) if (!is.null(a)) a else b

PAIRS      <- get_pairs(config)
padj_cut   <- config$de_analysis$padj_cutoff  %||% 0.05
lfc_cut    <- config$de_analysis$log2fc_cutoff %||% 0.0
project_id <- basename(config$output_dir)
method     <- config$de_analysis$method %||% "DESeq2"
species    <- config$species %||% "unknown"
group_var  <- config$de_analysis$group_variable %||% "condition"

cat(sprintf("Project: %s | Pairs: %d | padj < %.2f | |log2FC| > %.1f\n\n",
            project_id, length(PAIRS), padj_cut, lfc_cut))

# --- Load metadata for sample info ---
condition_samples <- list()   # condition -> character vector of sample_ids
n_total_samples   <- NA_integer_

meta_path <- config$metadata_path %||% NULL
if (!is.null(meta_path) && file.exists(meta_path)) {
  tryCatch({
    meta <- read.csv(meta_path, stringsAsFactors = FALSE)
    id_col <- intersect(c("sample_id", "sample", "SampleID", "ID"), colnames(meta))[1]
    if (is.na(id_col)) id_col <- colnames(meta)[1]
    if (group_var %in% colnames(meta)) {
      condition_samples <- split(meta[[id_col]], meta[[group_var]])
      n_total_samples   <- nrow(meta)
      cat(sprintf("  Metadata loaded: %d samples, %d conditions\n\n",
                  n_total_samples, length(condition_samples)))
    }
  }, error = function(e) cat(sprintf("  [WARN] Could not read metadata: %s\n\n", e$message)))
}

# --- 3. Collect per-pair results ---
pair_data <- list()

for (pair in PAIRS) {
  pair_dir <- file.path(output_dir, "pairwise", pair)
  de_file  <- file.path(pair_dir, "final_de_results.csv")

  # DE results
  if (!file.exists(de_file)) {
    cat(sprintf("  [WARN] %s: final_de_results.csv not found, skipping.\n", pair))
    next
  }
  de <- read.csv(de_file, check.names = FALSE)

  # Identify padj column (DESeq2: padj, edgeR: FDR, limma: adj.P.Val)
  padj_col <- intersect(c("padj", "FDR", "adj.P.Val"), colnames(de))[1]
  lfc_col  <- intersect(c("log2FoldChange", "logFC"), colnames(de))[1]
  gene_col <- intersect(c("symbol", "gene_symbol", "SYMBOL",
                           "gene", "Gene", "gene_id", "GeneID"), colnames(de))[1]
  # If no named gene column found, use the first non-empty column name (skip "" row-name col)
  if (is.na(gene_col)) {
    non_empty <- colnames(de)[nchar(colnames(de)) > 0]
    gene_col  <- if (length(non_empty) > 0) non_empty[1] else colnames(de)[1]
  }

  sig <- de[!is.na(de[[padj_col]]) & de[[padj_col]] < padj_cut &
              abs(de[[lfc_col]]) > lfc_cut, ]
  up   <- sig[sig[[lfc_col]] > 0, ]
  down <- sig[sig[[lfc_col]] < 0, ]

  # Top genes (sorted by |log2FC|)
  up_sorted   <- up[order(-abs(up[[lfc_col]])), ]
  down_sorted <- down[order(-abs(down[[lfc_col]])), ]
  top_n <- 10
  top_up   <- head(up_sorted,   top_n)
  top_down <- head(down_sorted, top_n)

  # GO BP enrichment
  read_go <- function(geneset) {
    f <- file.path(pair_dir, paste0("go_enrichment_", geneset, "_BP.csv"))
    if (!file.exists(f)) return(data.frame())
    tryCatch(read.csv(f, check.names = FALSE), error = function(e) data.frame())
  }
  go_up   <- read_go("up")
  go_down <- read_go("down")

  top_go_up   <- if (nrow(go_up)   > 0) head(go_up,   5) else data.frame()
  top_go_down <- if (nrow(go_down) > 0) head(go_down, 5) else data.frame()

  # Volcano plot (relative path from HTML location = output_dir root)
  volcano_rel <- file.path("pairwise", pair, "volcano_plot.png")

  pair_data[[pair]] <- list(
    pair        = pair,
    n_up        = nrow(up),
    n_down      = nrow(down),
    n_total     = nrow(sig),
    n_genes     = nrow(de),
    top_up      = top_up,
    top_down    = top_down,
    go_up       = top_go_up,
    go_down     = top_go_down,
    gene_col    = gene_col,
    lfc_col     = lfc_col,
    padj_col    = padj_col,
    volcano_rel = volcano_rel
  )

  cat(sprintf("  %s: %d up / %d down / %d total\n",
              pair, nrow(up), nrow(down), nrow(sig)))
}

cat("\n")

# --- 4. HTML helpers ---

# Inline CSS (reuses 02c_generate_global_qc_report.R style)
css <- '
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
        .header h1 { margin: 0 0 10px 0; font-size: 2.2em; }
        .header p  { margin: 4px 0; font-size: 1.05em; opacity: 0.9; }
        .section {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 28px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .section h2 {
            margin-top: 0;
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 16px;
            margin-top: 20px;
        }
        .summary-item {
            background: #f8f9fa;
            padding: 14px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }
        .summary-item .label { font-size: 0.88em; color: #6c757d; margin-bottom: 4px; }
        .summary-item .value { font-size: 1.5em; font-weight: bold; color: #212529; }
        table { width: 100%; border-collapse: collapse; margin: 12px 0; font-size: 0.92em; }
        th, td { padding: 10px 12px; text-align: left; border-bottom: 1px solid #dee2e6; }
        th { background-color: #667eea; color: white; font-weight: 600; }
        tr:hover { background-color: #f8f9fa; }
        .up   { color: #c0392b; font-weight: 600; }
        .down { color: #2980b9; font-weight: 600; }
        .badge {
            display: inline-block;
            padding: 3px 10px;
            border-radius: 12px;
            font-size: 0.82em;
            font-weight: 600;
        }
        .badge-up   { background: #fdecea; color: #c0392b; }
        .badge-down { background: #eaf4fb; color: #2980b9; }
        details {
            background: white;
            border-radius: 10px;
            margin-bottom: 16px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.08);
            overflow: hidden;
        }
        details > summary {
            padding: 18px 25px;
            cursor: pointer;
            font-size: 1.1em;
            font-weight: 600;
            color: #495057;
            background: #f8f9fa;
            border-radius: 10px;
            list-style: none;
            display: flex;
            align-items: center;
            gap: 12px;
        }
        details > summary::-webkit-details-marker { display: none; }
        details > summary::before {
            content: "▶";
            font-size: 0.75em;
            color: #667eea;
            transition: transform 0.2s;
        }
        details[open] > summary::before { transform: rotate(90deg); }
        details[open] > summary { border-radius: 10px 10px 0 0; }
        .detail-body { padding: 20px 25px; }
        .two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
        .col-block h4 { margin: 0 0 8px 0; color: #495057; }
        .plot-container { text-align: center; margin: 16px 0; }
        .plot-container img {
            max-width: 80%;
            height: auto;
            border-radius: 6px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        .note {
            background: #f8f9fa;
            padding: 12px 16px;
            border-left: 4px solid #764ba2;
            border-radius: 4px;
            font-size: 0.9em;
            color: #495057;
            margin: 12px 0;
        }
        .footer {
            text-align: center;
            margin-top: 40px;
            padding: 20px;
            color: #6c757d;
            font-size: 0.88em;
        }
        @media (max-width: 768px) {
            .two-col { grid-template-columns: 1fr; }
            .summary-grid { grid-template-columns: 1fr 1fr; }
        }
'

html_table <- function(df, gene_col, lfc_col, padj_col, direction) {
  if (nrow(df) == 0) return('<p style="color:#888;font-size:0.9em;">No significant genes.</p>')
  colour_class <- if (direction == "up") "up" else "down"
  rows <- apply(df, 1, function(r) {
    lfc_val  <- as.numeric(r[lfc_col])
    padj_val <- as.numeric(r[padj_col])
    gene_val <- r[gene_col]
    sprintf('<tr><td>%s</td><td class="%s">%+.3f</td><td>%.2e</td></tr>',
            gene_val, colour_class, lfc_val, padj_val)
  })
  paste0('<table><thead><tr><th>Gene</th><th>log2FC</th><th>padj</th></tr></thead><tbody>',
         paste(rows, collapse = "\n"), '</tbody></table>')
}

go_table <- function(go_df) {
  if (nrow(go_df) == 0) return('<p style="color:#888;font-size:0.9em;">No enriched terms.</p>')
  desc_col  <- intersect(c("Description", "description", "Term"), colnames(go_df))[1]
  padj_col  <- intersect(c("p.adjust", "qvalue", "FDR"), colnames(go_df))[1]
  count_col <- intersect(c("Count", "count"), colnames(go_df))[1]
  if (is.na(desc_col) || is.na(padj_col)) return('<p style="color:#888;font-size:0.9em;">(columns not found)</p>')
  rows <- apply(head(go_df, 5), 1, function(r) {
    cnt <- if (!is.na(count_col)) r[count_col] else "—"
    sprintf('<tr><td>%s</td><td>%s</td><td>%.3g</td></tr>',
            r[desc_col], cnt, as.numeric(r[padj_col]))
  })
  paste0('<table><thead><tr><th>GO Term</th><th>Count</th><th>p.adjust</th></tr></thead><tbody>',
         paste(rows, collapse = "\n"), '</tbody></table>')
}

# --- 5. Build HTML ---

# 5a. Header
total_pairs  <- length(pair_data)
total_degs   <- sum(sapply(pair_data, function(x) x$n_total))

html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DE-GO Summary Report - ', project_id, '</title>
    <style>', css, '</style>
</head>
<body>

<div class="header">
    <h1>&#128202; DE-GO Analysis Summary</h1>
    <p>Project: <strong>', project_id, '</strong></p>
    <p>Species: ', species, ' &nbsp;|&nbsp; Method: ', method,
    ' &nbsp;|&nbsp; padj &lt; ', padj_cut, ' &nbsp;|&nbsp; |log2FC| &gt; ', lfc_cut, '</p>
    <p>Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
</div>

<div class="section">
    <h2>&#128203; Analysis Overview</h2>
    <div class="summary-grid">
        <div class="summary-item">
            <div class="label">Comparisons</div>
            <div class="value">', total_pairs, '</div>
        </div>
        <div class="summary-item">
            <div class="label">Total Samples</div>
            <div class="value">', ifelse(is.na(n_total_samples), "—", n_total_samples), '</div>
        </div>
        <div class="summary-item">
            <div class="label">Total DEGs (all pairs)</div>
            <div class="value">', format(total_degs, big.mark = ","), '</div>
        </div>
        <div class="summary-item">
            <div class="label">padj cutoff</div>
            <div class="value">', padj_cut, '</div>
        </div>
        <div class="summary-item">
            <div class="label">|log2FC| cutoff</div>
            <div class="value">', lfc_cut, '</div>
        </div>
    </div>',
  # Sample breakdown table (if metadata available)
  if (length(condition_samples) > 0) {
    sample_rows <- paste(sapply(names(condition_samples), function(cond) {
      sids <- condition_samples[[cond]]
      sprintf('<tr><td><strong>%s</strong></td><td>%d</td><td style="font-size:0.88em;color:#555;">%s</td></tr>',
              cond, length(sids), paste(sids, collapse = ", "))
    }), collapse = "\n")
    paste0('
    <h3 style="color:#495057;margin-top:24px;margin-bottom:8px;">Sample Groups</h3>
    <table>
        <thead><tr><th>Condition</th><th>N</th><th>Sample IDs</th></tr></thead>
        <tbody>', sample_rows, '</tbody>
    </table>')
  } else "",
'
</div>
')

# 5b. Overview table
overview_rows <- sapply(pair_data, function(d) {
  go_up_term   <- if (nrow(d$go_up)   > 0) {
    desc_col <- intersect(c("Description","description","Term"), colnames(d$go_up))[1]
    if (!is.na(desc_col)) d$go_up[[desc_col]][1] else "—"
  } else "—"
  go_down_term <- if (nrow(d$go_down) > 0) {
    desc_col <- intersect(c("Description","description","Term"), colnames(d$go_down))[1]
    if (!is.na(desc_col)) d$go_down[[desc_col]][1] else "—"
  } else "—"

  sprintf(
    '<tr>
      <td><strong>%s</strong></td>
      <td class="up">%d</td>
      <td class="down">%d</td>
      <td><strong>%d</strong></td>
      <td style="font-size:0.88em;color:#555;">%s</td>
      <td style="font-size:0.88em;color:#555;">%s</td>
    </tr>',
    d$pair, d$n_up, d$n_down, d$n_total, go_up_term, go_down_term
  )
})

html <- paste0(html, '
<div class="section">
    <h2>&#128200; Results Overview</h2>
    <table>
        <thead>
            <tr>
                <th>Comparison</th>
                <th>Up</th>
                <th>Down</th>
                <th>Total DEGs</th>
                <th>Top GO UP (BP)</th>
                <th>Top GO DOWN (BP)</th>
            </tr>
        </thead>
        <tbody>
            ', paste(overview_rows, collapse = "\n"), '
        </tbody>
    </table>
    <div class="note">
        Cutoffs: padj &lt; ', padj_cut, ' &amp; |log2FC| &gt; ', lfc_cut, '.
        GO terms: top-ranked Biological Process pathway (p.adjust &lt; 0.05).
    </div>
</div>
')

# 5c. Global PCA (if available)
pca_rel <- "global_pca_plot.png"
pca_abs <- file.path(output_dir, pca_rel)
if (file.exists(pca_abs)) {
  html <- paste0(html, '
<div class="section">
    <h2>&#127760; Global PCA</h2>
    <div class="plot-container">
        <img src="', pca_rel, '" alt="Global PCA">
    </div>
</div>
')
}

# 5d. Per-comparison <details> sections
html <- paste0(html, '
<div class="section">
    <h2>&#128300; Per-Comparison Details</h2>
    <p style="color:#666;font-size:0.93em;">Click a comparison to expand its results.</p>
')

for (d in pair_data) {
  parts <- strsplit(d$pair, "_vs_")[[1]]
  compare_label <- parts[1]
  base_label    <- parts[2]

  # Sample info for this comparison
  compare_samples <- condition_samples[[compare_label]]
  base_samples    <- condition_samples[[base_label]]
  sample_info_html <- if (!is.null(compare_samples) && !is.null(base_samples)) {
    sprintf('
            <div class="note" style="margin-bottom:14px;">
                <strong>%s</strong> (n=%d): %s<br>
                <strong>%s</strong> (n=%d): %s
            </div>',
      compare_label, length(compare_samples), paste(compare_samples, collapse = ", "),
      base_label,    length(base_samples),    paste(base_samples,    collapse = ", "))
  } else ""

  html <- paste0(html, sprintf('
    <details>
        <summary>
            %s
            &nbsp;<span class="badge badge-up">&#8593; %d up</span>
            <span class="badge badge-down">&#8595; %d down</span>
        </summary>
        <div class="detail-body">

            %s

            <div class="plot-container">
                <img src="%s" alt="Volcano plot %s">
            </div>

            <div class="two-col">
                <div class="col-block">
                    <h4>&#128308; Top Up-regulated Genes</h4>
                    %s
                </div>
                <div class="col-block">
                    <h4>&#128309; Top Down-regulated Genes</h4>
                    %s
                </div>
            </div>

            <div class="two-col" style="margin-top:16px;">
                <div class="col-block">
                    <h4>GO Biological Process &mdash; UP</h4>
                    %s
                </div>
                <div class="col-block">
                    <h4>GO Biological Process &mdash; DOWN</h4>
                    %s
                </div>
            </div>

        </div>
    </details>
',
    d$pair, d$n_up, d$n_down,
    sample_info_html,
    d$volcano_rel, d$pair,
    html_table(d$top_up,   d$gene_col, d$lfc_col, d$padj_col, "up"),
    html_table(d$top_down, d$gene_col, d$lfc_col, d$padj_col, "down"),
    go_table(d$go_up),
    go_table(d$go_down)
  ))
}

html <- paste0(html, '\n</div>\n')  # close .section

# 5e. Time-series analysis (maSigPro) — de_analysis.time_series.enabled 시에만
ts_csv <- file.path(output_dir, "time_series", "time_series_significant_genes.csv")
if (isTRUE(config$de_analysis$time_series$enabled) && file.exists(ts_csv)) {
  ts_df <- read.csv(ts_csv, check.names = FALSE)
  ts_plot_rel <- file.path("time_series", "time_series_pattern_plot.png")
  ts_plot_abs <- file.path(output_dir, ts_plot_rel)

  ts_gene_col <- intersect(c("gene_symbol", "gene_id"), colnames(ts_df))[1]
  ts_top <- head(ts_df[order(ts_df$p_value), ], 10)
  ts_rows <- if (nrow(ts_top) > 0) {
    paste(apply(ts_top, 1, function(r) {
      sprintf('<tr><td>%s</td><td>%.2e</td><td>%.3f</td><td>%s</td></tr>',
              r[ts_gene_col], as.numeric(r["p_value"]), as.numeric(r["r_squared"]), r["cluster_id"])
    }), collapse = "\n")
  } else ""
  ts_table <- if (nrow(ts_top) > 0) {
    paste0('<table><thead><tr><th>Gene</th><th>p-value</th><th>R-squared</th><th>Cluster</th></tr></thead><tbody>',
           ts_rows, '</tbody></table>')
  } else '<p style="color:#888;font-size:0.9em;">No significant genes.</p>'

  html <- paste0(html, sprintf('
<div class="section">
    <h2>&#8987; Time-Series Analysis (maSigPro)</h2>
    <div class="summary-grid">
        <div class="summary-item">
            <div class="label">Significant genes</div>
            <div class="value">%d</div>
        </div>
        <div class="summary-item">
            <div class="label">Pattern clusters</div>
            <div class="value">%d</div>
        </div>
    </div>
    %s
    <h3 style="color:#495057;margin-top:20px;">Top genes (by p-value)</h3>
    %s
</div>
', nrow(ts_df),
   if (nrow(ts_df) > 0) length(unique(ts_df$cluster_id)) else 0,
   if (file.exists(ts_plot_abs)) sprintf('<div class="plot-container"><img src="%s" alt="Time-series pattern plot"></div>', ts_plot_rel) else "",
   ts_table
  ))
}

# 5f. Coexpression module analysis — de_analysis.coexpression_modules.enabled 시에만
cm_csv <- file.path(output_dir, "coexpression_modules", "coexpression_module_assignments.csv")
if (isTRUE(config$de_analysis$coexpression_modules$enabled) && file.exists(cm_csv)) {
  cm_df <- read.csv(cm_csv, check.names = FALSE)
  cm_heatmap_rel <- file.path("coexpression_modules", "coexpression_module_heatmap.png")
  cm_heatmap_abs <- file.path(output_dir, cm_heatmap_rel)
  cm_pattern_rel <- file.path("coexpression_modules", "coexpression_pattern_plot.png")
  cm_pattern_abs <- file.path(output_dir, cm_pattern_rel)

  module_counts <- if (nrow(cm_df) > 0) table(cm_df$module_id) else integer(0)
  module_rows <- paste(sapply(names(module_counts), function(m) {
    sprintf('<tr><td>Module %s</td><td>%d</td></tr>', m, module_counts[[m]])
  }), collapse = "\n")
  module_table <- if (length(module_counts) > 0) {
    paste0('<table><thead><tr><th>Module</th><th>Gene count</th></tr></thead><tbody>', module_rows, '</tbody></table>')
  } else '<p style="color:#888;font-size:0.9em;">No modules found.</p>'

  html <- paste0(html, sprintf('
<div class="section">
    <h2>&#129504; Coexpression Modules</h2>
    <div class="summary-grid">
        <div class="summary-item">
            <div class="label">Genes clustered</div>
            <div class="value">%d</div>
        </div>
        <div class="summary-item">
            <div class="label">Modules</div>
            <div class="value">%d</div>
        </div>
    </div>
    %s
    <h3 style="color:#495057;margin-top:20px;">Module sizes</h3>
    %s
</div>
', nrow(cm_df), length(module_counts),
   if (file.exists(cm_heatmap_abs)) sprintf('<div class="plot-container"><img src="%s" alt="Coexpression module heatmap"></div>', cm_heatmap_rel) else "",
   if (file.exists(cm_pattern_abs)) paste0(module_table, sprintf('<div class="plot-container"><img src="%s" alt="Coexpression module pattern plot"></div>', cm_pattern_rel)) else module_table
  ))
}

# 5g. Footer
html <- paste0(html, '
<div class="footer">
    <p>Generated by RNA-Seq DE-GO Analysis Pipeline</p>
    <p>Configuration: ', opt$config, '</p>
</div>

</body>
</html>
')

# --- 6. Write output ---
writeLines(html, output_file)

cat("=== Summary Report Generation Complete ===\n")
cat(paste0("HTML report saved to: ", output_file, "\n"))
cat(paste0("View in browser: file://", normalizePath(output_file), "\n\n"))
