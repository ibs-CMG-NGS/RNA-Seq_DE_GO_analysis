#!/usr/bin/env Rscript
# 파일 경로: src/analysis/08_generate_methods_section.R
# 목적: DE-GO 분석에 사용된 모든 도구, 파라미터, 버전 정보를 담은 Methods 섹션 HTML 생성
# 사용법:
#   Rscript src/analysis/08_generate_methods_section.R \
#     --config configs/config_PROJECT.yml \
#     --output-dir output/PROJECT \
#     --output output/PROJECT/methods_section.html

suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !identical(a, "")) a else b

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

option_list <- list(
  make_option(c("-c", "--config"),     type = "character",
              help = "Path to config YAML file", metavar = "FILE"),
  make_option(c("-d", "--output-dir"), type = "character",
              help = "Pipeline output directory (for reading result files)", metavar = "DIR"),
  make_option(c("-o", "--output"),     type = "character",
              help = "Output HTML file path", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$config) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("--config and --output are required.", call. = FALSE)
}

cat("\n=== DE-GO Methods Section Generator ===\n")
cat(paste0("Config:     ", opt$config, "\n"))
cat(paste0("Output:     ", opt$output, "\n\n"))

# ---------------------------------------------------------------------------
# Load config
# ---------------------------------------------------------------------------

cfg <- yaml.load_file(opt$config)

project_id   <- basename(cfg$output_dir %||% opt$`output-dir` %||% "PROJECT")
species      <- cfg$species %||% "unknown"
gene_id_type <- cfg$gene_id_type %||% "ENSEMBL"

de_cfg       <- cfg$de_analysis  %||% list()
enr_cfg      <- cfg$enrichment   %||% list()
qc_cfg       <- cfg$qc_plots     %||% list()

de_method    <- de_cfg$method              %||% "DESeq2"
design_f     <- de_cfg$design_formula      %||% "~ condition"
group_var    <- de_cfg$group_variable      %||% "condition"
padj_cut     <- de_cfg$padj_cutoff         %||% 0.05
lfc_cut      <- de_cfg$log2fc_cutoff       %||% 1.0
norm_strat   <- (de_cfg$advanced_options   %||% list())$pairwise_normalization %||% "global"
prefilter_n  <- (de_cfg$advanced_options   %||% list())$prefilter_min_count    %||% 1
prefilter_s  <- (de_cfg$advanced_options   %||% list())$prefilter_min_samples  %||% "auto"
vst_blind    <- (de_cfg$advanced_options   %||% list())$vst_blind              %||% FALSE
pca_ntop     <- (de_cfg$advanced_options   %||% list())$pca_ntop               %||% 500

pairs_raw    <- de_cfg$pairwise_comparisons %||% list()
pairs        <- vapply(pairs_raw, function(p) paste0(p[[1]], "_vs_", p[[2]]), character(1))

go_ontologies  <- paste(enr_cfg$go_ontologies  %||% c("BP","CC","MF"), collapse = ", ")
gene_lists     <- paste(enr_cfg$gene_lists     %||% c("total","up","down"), collapse = ", ")
go_pval        <- enr_cfg$pvalue_cutoff %||% 0.05
go_qval        <- enr_cfg$qvalue_cutoff %||% 0.25
min_gs         <- enr_cfg$min_gs_size   %||% 5
max_gs         <- enr_cfg$max_gs_size   %||% 500
min_gene_count <- enr_cfg$min_gene_count %||% 2
plot_top_n     <- enr_cfg$plot_top_n    %||% 15

# organism DB
species_lower <- tolower(species)
db_cfg        <- (cfg$databases %||% list())[[species_lower]] %||% list()
organism_db   <- db_cfg$organism_db %||% if (species_lower == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
kegg_code     <- db_cfg$kegg_code   %||% if (species_lower == "mouse") "mmu" else "hsa"

# ---------------------------------------------------------------------------
# Dataset dimensions (from count matrix)
# ---------------------------------------------------------------------------

n_samples  <- NA_integer_
n_genes    <- NA_integer_
conditions <- list()

if (!is.null(cfg$count_data_path) && file.exists(cfg$count_data_path)) {
  tryCatch({
    count_data <- read.csv(cfg$count_data_path, row.names = 1,
                           check.names = FALSE, nrows = 5)
    n_samples  <- ncol(count_data)
    # get full row count quickly
    n_genes    <- length(readLines(cfg$count_data_path)) - 1L
    cat(paste0("[INFO] Count data: ", n_genes, " genes × ", n_samples, " samples\n"))
  }, error = function(e) {
    cat(paste0("[WARN] Could not read count data: ", e$message, "\n"))
  })
}

if (!is.null(cfg$metadata_path) && file.exists(cfg$metadata_path)) {
  tryCatch({
    meta <- read.csv(cfg$metadata_path, row.names = 1)
    if (group_var %in% colnames(meta)) {
      conditions <- split(rownames(meta), meta[[group_var]])
    }
    if (is.na(n_samples)) n_samples <- nrow(meta)
    cat(paste0("[INFO] Conditions: ", paste(names(conditions), collapse = ", "), "\n"))
  }, error = function(e) {
    cat(paste0("[WARN] Could not read metadata: ", e$message, "\n"))
  })
}

# ---------------------------------------------------------------------------
# Package versions (actual installed)
# ---------------------------------------------------------------------------

pkg_version <- function(pkg) {
  tryCatch(as.character(packageVersion(pkg)), error = function(e) "N/A")
}

r_ver         <- paste0(R.version$major, ".", R.version$minor)
deseq2_ver    <- pkg_version("DESeq2")
edger_ver     <- pkg_version("edgeR")
limma_ver     <- pkg_version("limma")
clusterp_ver  <- pkg_version("clusterProfiler")
annodbi_ver   <- pkg_version("AnnotationDbi")
orgdb_ver     <- pkg_version(organism_db)
enrichplot_ver<- pkg_version("enrichplot")
ggplot2_ver   <- pkg_version("ggplot2")
snakemake_ver <- tryCatch({
  out <- system2("snakemake", "--version", stdout = TRUE, stderr = FALSE)
  trimws(out[1])
}, error = function(e) "N/A")

# Primary DE package version
de_pkg_ver <- switch(de_method,
  "DESeq2" = deseq2_ver,
  "edgeR"  = edger_ver,
  "limma"  = limma_ver,
  pkg_version(de_method)
)

# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

html_table <- function(headers, rows, col_widths = NULL) {
  th <- paste0("<th>", headers, "</th>", collapse = "\n          ")
  tr <- vapply(rows, function(r) {
    tds <- paste0("<td>", r, "</td>", collapse = "\n          ")
    paste0("        <tr>\n          ", tds, "\n        </tr>")
  }, character(1))
  paste0(
    '        <table>\n',
    '          <thead><tr>\n          ', th, '\n          </tr></thead>\n',
    '          <tbody>\n',
    paste(tr, collapse = "\n"),
    '\n          </tbody>\n        </table>'
  )
}

summary_card <- function(label, value, color = "#667eea") {
  paste0(
    '<div class="summary-item" style="border-left-color:', color, ';">\n',
    '  <div class="label">', label, '</div>\n',
    '  <div class="value">', value, '</div>\n',
    '</div>'
  )
}

section <- function(title, icon = "", content) {
  paste0(
    '<div class="section">\n',
    '  <h2>', icon, ' ', title, '</h2>\n',
    content,
    '\n</div>'
  )
}

detail_block <- function(summary_text, content) {
  paste0(
    '<details>\n<summary>', summary_text, '</summary>\n',
    '<div class="details-content">\n', content, '\n</div>\n</details>'
  )
}

note_box <- function(text) {
  paste0('<div class="note">', text, '</div>')
}

# ---------------------------------------------------------------------------
# CSS (same palette as 02c)
# ---------------------------------------------------------------------------

css <- '
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      line-height: 1.6;
      max-width: 1100px;
      margin: 0 auto;
      padding: 20px;
      background-color: #f5f5f5;
      color: #212529;
    }
    .header {
      background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
      color: white;
      padding: 30px;
      border-radius: 10px;
      margin-bottom: 30px;
      box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    .header h1 { margin: 0 0 8px 0; font-size: 2em; }
    .header p  { margin: 4px 0; font-size: 1em; opacity: 0.9; }
    .summary-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
      gap: 16px;
      margin: 20px 0;
    }
    .summary-item {
      background: #f8f9fa;
      padding: 14px;
      border-radius: 8px;
      border-left: 4px solid #667eea;
    }
    .summary-item .label { font-size: 0.85em; color: #6c757d; margin-bottom: 4px; }
    .summary-item .value { font-size: 1.4em; font-weight: bold; }
    .section {
      background: white;
      padding: 25px;
      border-radius: 10px;
      margin-bottom: 24px;
      box-shadow: 0 2px 4px rgba(0,0,0,0.08);
    }
    .section h2 {
      margin-top: 0;
      color: #667eea;
      border-bottom: 2px solid #667eea;
      padding-bottom: 8px;
      font-size: 1.3em;
    }
    .section h3 { color: #495057; font-size: 1.05em; margin-top: 20px; }
    table {
      width: 100%;
      border-collapse: collapse;
      margin: 14px 0;
      font-size: 0.95em;
    }
    th, td { padding: 10px 12px; text-align: left; border-bottom: 1px solid #dee2e6; }
    th { background-color: #667eea; color: white; font-weight: 600; }
    tr:hover { background-color: #f8f9fa; }
    code {
      background: #f0f0f0;
      padding: 2px 6px;
      border-radius: 3px;
      font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
      font-size: 0.9em;
    }
    .note {
      background: #f8f9fa;
      border-left: 4px solid #764ba2;
      padding: 12px 16px;
      border-radius: 4px;
      margin: 14px 0;
      font-size: 0.93em;
      color: #495057;
    }
    details { margin: 10px 0; }
    summary {
      cursor: pointer;
      padding: 10px 14px;
      background: #f8f9fa;
      border-radius: 6px;
      font-weight: 600;
      color: #495057;
      border: 1px solid #dee2e6;
      list-style: none;
    }
    summary::-webkit-details-marker { display: none; }
    summary::before { content: "▶  "; font-size: 0.8em; }
    details[open] summary::before { content: "▼  "; }
    .details-content { padding: 14px 4px 4px 4px; }
    .footer {
      text-align: center;
      margin-top: 40px;
      padding: 20px;
      color: #6c757d;
      font-size: 0.88em;
    }
    .ref-list { padding-left: 20px; }
    .ref-list li { margin-bottom: 6px; font-size: 0.93em; }
    .badge {
      display: inline-block;
      padding: 2px 8px;
      border-radius: 4px;
      font-size: 0.85em;
      font-weight: 600;
      background: #667eea;
      color: white;
      margin-right: 4px;
    }
'

# ---------------------------------------------------------------------------
# Build HTML sections
# ---------------------------------------------------------------------------

# ── Header ──────────────────────────────────────────────────────────────────
header_html <- paste0(
  '<div class="header">\n',
  '  <h1>Methods</h1>\n',
  '  <p>Differential Expression &amp; Gene Ontology Analysis</p>\n',
  '  <p>Project: <strong>', project_id, '</strong></p>\n',
  '  <p>Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>\n',
  '</div>'
)

# ── Overview cards ───────────────────────────────────────────────────────────
card_species  <- summary_card("Species", paste0(toupper(substring(species,1,1)), substring(species,2)))
card_method   <- summary_card("DE Method", de_method, "#764ba2")
card_samples  <- if (!is.na(n_samples)) summary_card("Samples", n_samples) else NULL
card_genes    <- if (!is.na(n_genes))   summary_card("Genes (pre-filter)", format(n_genes, big.mark=","), "#28a745") else NULL
card_pairs    <- summary_card("Comparisons", length(pairs), "#fd7e14")
card_padj     <- summary_card("padj cutoff", padj_cut, "#dc3545")
card_lfc      <- summary_card("|log2FC| cutoff", lfc_cut, "#17a2b8")

cards <- paste(Filter(Negate(is.null), list(
  card_species, card_method, card_samples, card_genes, card_pairs, card_padj, card_lfc
)), collapse = "\n")

overview_section <- section("Dataset Overview", "&#128202;", paste0(
  '<div class="summary-grid">\n', cards, '\n</div>\n',
  if (length(conditions) > 0) {
    cond_rows <- lapply(names(conditions), function(g)
      c(g, length(conditions[[g]]), paste(conditions[[g]], collapse = ", ")))
    html_table(c("Condition", "N samples", "Sample IDs"), cond_rows)
  } else "",
  if (length(pairs) > 0) {
    pair_rows <- lapply(pairs, function(p) {
      parts <- strsplit(p, "_vs_")[[1]]
      c(parts[1], parts[2], p)
    })
    paste0("<h3>Pairwise Comparisons</h3>",
           html_table(c("Treatment", "Control", "Comparison ID"), pair_rows))
  } else ""
))

# ── DE Analysis ─────────────────────────────────────────────────────────────
de_note <- switch(de_method,
  "DESeq2" = paste0(
    "Raw counts were modelled using a negative binomial distribution. ",
    "Size factors for library normalization were estimated using the median-of-ratios method. ",
    "Gene-wise dispersion was estimated using maximum-likelihood and then empirically Bayes-shrunk toward a fitted trend. ",
    if (norm_strat == "global") "Normalization was performed globally across all samples (single DESeqDataSet). " else
      "Normalization was performed separately for each pairwise comparison. ",
    "Log2 fold changes were shrunk using the apeglim estimator."
  ),
  "edgeR" = paste0(
    "Raw counts were modelled using a negative binomial distribution (edgeR). ",
    "Normalization factors were computed using the TMM method. ",
    "Dispersion was estimated using the Cox-Reid adjusted profile likelihood."
  ),
  "limma" = paste0(
    "Counts were transformed using voom (limma-voom), which estimates the mean-variance ",
    "relationship and applies precision weights. Linear models were fit per gene."
  ),
  paste0("Differential expression was performed using ", de_method, ".")
)

prefilter_str <- if (identical(prefilter_s, "auto")) {
  paste0("Genes with fewer than ", prefilter_n, " count in at least [n/group] samples were removed prior to analysis.")
} else {
  paste0("Genes with fewer than ", prefilter_n, " count in at least ", prefilter_s, " samples were removed prior to analysis.")
}

de_params_table <- html_table(
  c("Parameter", "Value", "Description"),
  list(
    c("Method",                 de_method,      "Statistical framework for DE testing"),
    c("Design formula",         paste0("<code>", design_f, "</code>"), "Model formula"),
    c("Group variable",         group_var,       "Metadata column for grouping"),
    c("Normalization strategy", norm_strat,      "How samples are normalized across comparisons"),
    c("Pre-filter (min counts)",prefilter_n,     "Minimum count per gene to retain"),
    c("Pre-filter (min samples)",prefilter_s,    "Minimum number of samples passing count threshold"),
    c("VST blind",              as.character(vst_blind), "Blind to design when estimating VST dispersion"),
    c("PCA top genes",          pca_ntop,        "Number of most-variable genes used for PCA")
  )
)

de_cutoff_table <- html_table(
  c("Cutoff", "Value", "Applied to"),
  list(
    c("Adjusted p-value (padj)", padj_cut, "FDR-adjusted Wald test p-value (Benjamini-Hochberg)"),
    c("log2 fold change",        lfc_cut,  "Absolute log2(treatment / control)")
  )
)

de_section <- section("Differential Expression Analysis", "&#128202;", paste0(
  '<p>', de_note, '</p>',
  '<p>', prefilter_str, '</p>',
  '<h3>Analysis Parameters</h3>', de_params_table,
  '<h3>Significance Cutoffs</h3>', de_cutoff_table
))

# ── GO/KEGG Enrichment ───────────────────────────────────────────────────────
enr_note <- paste0(
  "Over-representation analysis (ORA) was performed on differentially expressed gene sets using ",
  "<strong>clusterProfiler</strong> (v", clusterp_ver, ") [Wu 2021]. ",
  "Gene symbols were mapped to Entrez IDs via <strong>", organism_db, "</strong> (v", orgdb_ver, "). ",
  "GO enrichment was tested across all three ontologies (Biological Process, Cellular Component, ",
  "Molecular Function). KEGG pathway enrichment was additionally performed for each gene set."
)

enr_params_table <- html_table(
  c("Parameter", "Value", "Description"),
  list(
    c("Gene ontologies",     go_ontologies,       "GO sub-ontologies tested"),
    c("Gene lists",          gene_lists,          "Gene set directions tested"),
    c("Organism DB",         organism_db,         "Annotation database for ID mapping"),
    c("KEGG organism code",  kegg_code,           "KEGG species code"),
    c("p-value cutoff",      go_pval,             "Maximum raw p-value for terms to report"),
    c("q-value cutoff",      go_qval,             "Maximum FDR-adjusted p-value for terms to report"),
    c("Min gene set size",   min_gs,              "Minimum number of annotated genes in a term"),
    c("Max gene set size",   max_gs,              "Maximum number of annotated genes in a term"),
    c("Min DE gene count",   min_gene_count,      "Minimum overlap between DE genes and term genes"),
    c("Plot top N terms",    plot_top_n,          "Number of top terms shown in dot/bar plots")
  )
)

enr_section <- section("Functional Enrichment Analysis", "&#127381;", paste0(
  '<p>', enr_note, '</p>',
  '<h3>Enrichment Parameters</h3>', enr_params_table
))

# ── Software Versions ────────────────────────────────────────────────────────
sw_rows <- list(
  c("R",               r_ver,         "Statistical computing environment"),
  c(de_method,         de_pkg_ver,    "Differential expression testing")
)
if (de_method != "DESeq2" && deseq2_ver != "N/A")
  sw_rows <- c(sw_rows, list(c("DESeq2",  deseq2_ver,     "Normalization / VST")))
if (de_method != "edgeR"  && edger_ver  != "N/A")
  sw_rows <- c(sw_rows, list(c("edgeR",   edger_ver,      "Normalization support")))
if (de_method != "limma"  && limma_ver  != "N/A")
  sw_rows <- c(sw_rows, list(c("limma",   limma_ver,      "Normalization support")))
sw_rows <- c(sw_rows, list(
  c("clusterProfiler", clusterp_ver,  "GO and KEGG enrichment analysis"),
  c("AnnotationDbi",   annodbi_ver,   "Bioconductor annotation interface"),
  c(organism_db,       orgdb_ver,     paste0("Gene annotation (", species, ")")),
  c("enrichplot",      enrichplot_ver,"Enrichment result visualization"),
  c("ggplot2",         ggplot2_ver,   "Data visualization"),
  c("Snakemake",       snakemake_ver, "Workflow management")
))

sw_section <- section("Software Versions", "&#128295;", paste0(
  html_table(c("Package / Tool", "Version", "Purpose"), sw_rows),
  note_box(paste0(
    "All R packages were managed via conda (environment: <code>rna-seq-de-go-analysis</code>). ",
    "The full environment specification is available at <code>environment.yml</code> ",
    "in the pipeline root directory."
  ))
))

# ── References ───────────────────────────────────────────────────────────────
refs <- c(
  paste0("Love MI, et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. <em>Genome Biology</em>, 15:550."),
  paste0("Robinson MD, et al. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. <em>Bioinformatics</em>, 26(1):139–140."),
  paste0("Ritchie ME, et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. <em>Nucleic Acids Research</em>, 43(7):e47."),
  paste0("Wu T, et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. <em>Innovation</em>, 2(3):100141."),
  paste0("Zhu A, et al. (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. <em>Bioinformatics</em>, 35(12):2084–2092. [apeglim shrinkage]"),
  paste0("Köster J &amp; Rahmann S. (2012). Snakemake — a scalable bioinformatics workflow engine. <em>Bioinformatics</em>, 28(19):2520–2522.")
)
ref_li <- paste0("<li>", refs, "</li>", collapse = "\n")
ref_section <- section("References", "&#128218;",
  paste0('<ol class="ref-list">\n', ref_li, '\n</ol>')
)

# ── Appendix: Full Config ────────────────────────────────────────────────────
flatten_config <- function(x, prefix = "") {
  rows <- list()
  for (k in names(x)) {
    full_key <- if (nchar(prefix) > 0) paste0(prefix, ".", k) else k
    v <- x[[k]]
    if (is.list(v)) {
      rows <- c(rows, flatten_config(v, full_key))
    } else {
      display <- paste(as.character(v), collapse = ", ")
      if (nchar(display) > 100) display <- paste0(substr(display, 1, 97), "...")
      rows <- c(rows, list(c(
        paste0("<code>", full_key, "</code>"),
        paste0("<code>", display, "</code>")
      )))
    }
  }
  rows
}

appendix_rows   <- flatten_config(cfg)
appendix_table  <- html_table(c("Parameter", "Value"), appendix_rows)
appendix_section <- section("Appendix: Full Configuration", "&#128196;",
  detail_block(
    paste0("Show all parameters (", length(appendix_rows), " entries)"),
    paste0(
      note_box(paste0(
        "Configuration file: <code>", opt$config, "</code>. ",
        "This appendix provides a complete record of all parameters for reproducibility."
      )),
      appendix_table
    )
  )
)

# ---------------------------------------------------------------------------
# Assemble full HTML
# ---------------------------------------------------------------------------

html <- paste0(
'<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Methods — ', project_id, '</title>
  <style>', css, '</style>
</head>
<body>
', header_html, '\n',
overview_section, '\n',
de_section, '\n',
enr_section, '\n',
sw_section, '\n',
ref_section, '\n',
appendix_section, '\n',
'<div class="footer">
  <p>Generated by <code>08_generate_methods_section.R</code> &mdash;
  <strong>Note to researcher:</strong> Review all parameters and retain only those
  relevant to your manuscript. Adjust descriptions as needed for your target journal.</p>
  <p>Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
</div>
</body>
</html>'
)

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------

output_path <- opt$output
output_dir  <- dirname(output_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

writeLines(html, output_path)
cat(paste0("[INFO] Methods section written to: ", output_path, "\n"))
