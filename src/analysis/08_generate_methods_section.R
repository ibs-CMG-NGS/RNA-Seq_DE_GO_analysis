#!/usr/bin/env Rscript
# 파일 경로: src/analysis/08_generate_methods_section.R
# 목적: DE-GO 분석에 사용된 모든 도구, 파라미터, 버전 정보를 담은 Methods 섹션 Markdown 생성
# 사용법:
#   Rscript src/analysis/08_generate_methods_section.R \
#     --config configs/config_PROJECT.yml \
#     --output-dir output/PROJECT \
#     --output output/PROJECT/methods_section.md

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
              help = "Pipeline output directory", metavar = "DIR"),
  make_option(c("-o", "--output"),     type = "character",
              help = "Output Markdown file path", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$config) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("--config and --output are required.", call. = FALSE)
}

cat("\n=== DE-GO Methods Section Generator ===\n")
cat(paste0("Config: ", opt$config, "\n"))
cat(paste0("Output: ", opt$output, "\n\n"))

# ---------------------------------------------------------------------------
# Load config
# ---------------------------------------------------------------------------

cfg <- yaml.load_file(opt$config)

project_id   <- basename(cfg$output_dir %||% opt$`output-dir` %||% "PROJECT")
project_name <- cfg$project_name %||% project_id
species      <- cfg$species      %||% "unknown"
gene_id_type <- cfg$gene_id_type %||% "ENSEMBL"

de_cfg  <- cfg$de_analysis %||% list()
enr_cfg <- cfg$enrichment  %||% list()
adv     <- de_cfg$advanced_options %||% list()

de_method    <- de_cfg$method           %||% "DESeq2"
design_f     <- de_cfg$design_formula   %||% "~ condition"
group_var    <- de_cfg$group_variable   %||% "condition"
padj_cut     <- de_cfg$padj_cutoff      %||% 0.05
lfc_cut      <- de_cfg$log2fc_cutoff    %||% 1.0
norm_strat   <- adv$pairwise_normalization %||% "global"
prefilter_n  <- adv$prefilter_min_count    %||% 1
prefilter_s  <- adv$prefilter_min_samples  %||% "auto"
vst_blind    <- adv$vst_blind              %||% FALSE
pca_ntop     <- adv$pca_ntop               %||% 500

pairs_raw <- de_cfg$pairwise_comparisons %||% list()
pairs     <- vapply(pairs_raw, function(p) paste0(p[[1]], "_vs_", p[[2]]), character(1))

go_ontologies  <- paste(enr_cfg$go_ontologies %||% c("BP","CC","MF"), collapse = ", ")
gene_lists     <- paste(enr_cfg$gene_lists    %||% c("total","up","down"), collapse = ", ")
go_pval        <- enr_cfg$pvalue_cutoff %||% 0.05
go_qval        <- enr_cfg$qvalue_cutoff %||% 0.25
min_gs         <- enr_cfg$min_gs_size   %||% 5
max_gs         <- enr_cfg$max_gs_size   %||% 500
min_gene_count <- enr_cfg$min_gene_count %||% 2
plot_top_n     <- enr_cfg$plot_top_n     %||% 15

species_lower <- tolower(species)
db_cfg        <- (cfg$databases %||% list())[[species_lower]] %||% list()
organism_db   <- db_cfg$organism_db %||% if (species_lower == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
kegg_code     <- db_cfg$kegg_code   %||% if (species_lower == "mouse") "mmu" else "hsa"

# ---------------------------------------------------------------------------
# Dataset dimensions
# ---------------------------------------------------------------------------

n_samples  <- NA_integer_
n_genes    <- NA_integer_
conditions <- list()

if (!is.null(cfg$count_data_path) && file.exists(cfg$count_data_path)) {
  tryCatch({
    n_genes   <- length(readLines(cfg$count_data_path)) - 1L
    hdr       <- read.csv(cfg$count_data_path, nrows = 1, check.names = FALSE)
    n_samples <- ncol(hdr)
    cat(paste0("[INFO] Count data: ", n_genes, " genes x ", n_samples, " samples\n"))
  }, error = function(e) cat(paste0("[WARN] count data: ", e$message, "\n")))
}

if (!is.null(cfg$metadata_path) && file.exists(cfg$metadata_path)) {
  tryCatch({
    meta <- read.csv(cfg$metadata_path, row.names = 1)
    if (group_var %in% colnames(meta))
      conditions <- split(rownames(meta), meta[[group_var]])
    if (is.na(n_samples)) n_samples <- nrow(meta)
    cat(paste0("[INFO] Conditions: ", paste(names(conditions), collapse = ", "), "\n"))
  }, error = function(e) cat(paste0("[WARN] metadata: ", e$message, "\n")))
}

# ---------------------------------------------------------------------------
# Package versions (actual installed)
# ---------------------------------------------------------------------------

pkg_ver <- function(pkg) tryCatch(as.character(packageVersion(pkg)), error = function(e) "N/A")

r_ver          <- paste0(R.version$major, ".", R.version$minor)
deseq2_ver     <- pkg_ver("DESeq2")
edger_ver      <- pkg_ver("edgeR")
limma_ver      <- pkg_ver("limma")
clusterp_ver   <- pkg_ver("clusterProfiler")
annodbi_ver    <- pkg_ver("AnnotationDbi")
orgdb_ver      <- pkg_ver(organism_db)
enrichplot_ver <- pkg_ver("enrichplot")
ggplot2_ver    <- pkg_ver("ggplot2")
snakemake_ver  <- tryCatch(
  trimws(system2("snakemake", "--version", stdout = TRUE, stderr = FALSE)[1]),
  error = function(e) "N/A"
)

de_pkg_ver <- switch(de_method,
  "DESeq2" = deseq2_ver, "edgeR" = edger_ver, "limma" = limma_ver,
  pkg_ver(de_method)
)

# ---------------------------------------------------------------------------
# Markdown builder helpers
# ---------------------------------------------------------------------------

lines <- character(0)

h <- function(level, text) {
  lines <<- c(lines, paste0(strrep("#", level), " ", text), "")
}
p <- function(text = "") lines <<- c(lines, text)
blank <- function() lines <<- c(lines, "")
li <- function(text) lines <<- c(lines, paste0("- ", text))

# Contiguous Markdown table (no blank lines between rows)
md_table <- function(headers, rows) {
  ncol <- length(headers)
  widths <- nchar(headers)
  for (row in rows)
    for (i in seq_along(row))
      if (i <= ncol) widths[i] <- max(widths[i], nchar(as.character(row[[i]])), 3)

  fmt_row <- function(cells) {
    padded <- vapply(seq_len(ncol), function(i) {
      formatC(as.character(cells[[i]]), width = -widths[i], flag = "-")
    }, character(1))
    paste0("| ", paste(padded, collapse = " | "), " |")
  }
  sep <- paste0("| ", paste(strrep("-", widths), collapse = " | "), " |")

  lines <<- c(lines, fmt_row(headers), sep)
  for (row in rows) lines <<- c(lines, fmt_row(row))
  lines <<- c(lines, "")
}

# ---------------------------------------------------------------------------
# Build Markdown document
# ---------------------------------------------------------------------------

today <- format(Sys.Date(), "%Y-%m-%d")

# Species capitalised
species_cap <- paste0(toupper(substring(species, 1, 1)), substring(species, 2))

# Sample/condition string
if (length(conditions) > 0) {
  cond_str <- paste(
    vapply(names(conditions), function(g) paste0(g, " (n=", length(conditions[[g]]), ")"), character(1)),
    collapse = ", "
  )
} else {
  cond_str <- if (!is.na(n_samples)) paste0("n=", n_samples) else "unknown"
}

# ── Title ──────────────────────────────────────────────────────────────────
h(1, "RNA-seq Differential Expression & GO Analysis Methods")
p(paste0("**Project:** ", project_name, " (`", project_id, "`)  "))
p(paste0("**Generated:** ", today, "  "))
p(paste0("**Species:** *", species_cap, "*  "))
if (!is.na(n_samples))
  p(paste0("**Samples:** ", n_samples, " (", cond_str, ")  "))
if (!is.na(n_genes))
  p(paste0("**Genes (pre-filter):** ", format(n_genes, big.mark = ","), "  "))
blank()
p("> **Note to researcher:** This document describes all tools, versions, and parameter")
p("> settings used in the analysis. Please review each section and retain only the")
p("> information relevant to your manuscript. Suggested citation keys are in the References section.")
blank()

# ── 1. Input Data ──────────────────────────────────────────────────────────
h(2, "1. Input Data")
p(paste0(
  "Gene-level read count matrices and sample metadata were used as input. ",
  "Gene identifiers are in **", gene_id_type, "** format."
))
blank()
if (length(conditions) > 0) {
  cond_rows <- lapply(names(conditions), function(g)
    list(g, length(conditions[[g]]), paste(conditions[[g]], collapse = ", ")))
  md_table(c("Condition", "N samples", "Sample IDs"), cond_rows)
}
if (length(pairs) > 0) {
  p("**Pairwise comparisons:**")
  blank()
  pair_rows <- lapply(pairs, function(p_str) {
    parts <- strsplit(p_str, "_vs_")[[1]]
    list(parts[1], parts[2], p_str)
  })
  md_table(c("Treatment", "Control", "Comparison ID"), pair_rows)
}

# ── 2. Differential Expression Analysis ────────────────────────────────────
h(2, "2. Differential Expression Analysis")

de_narrative <- switch(de_method,
  "DESeq2" = paste0(
    "Differential expression analysis was performed using **DESeq2** (v", deseq2_ver,
    ") [Love 2014], which models raw counts with a negative binomial distribution. ",
    "Library size normalization used the median-of-ratios method. ",
    "Gene-wise dispersion estimates were maximum-likelihood fitted and then empirically Bayes-shrunk ",
    "toward a fitted trend. ",
    if (norm_strat == "global")
      "Normalization was performed globally across all samples in a single DESeqDataSet. "
    else
      "Normalization was performed independently for each pairwise comparison. ",
    "Log2 fold changes were shrunk using the **apeglim** estimator [Zhu 2019] to reduce noise for low-count genes. ",
    "Statistical significance was assessed using the Wald test."
  ),
  "edgeR" = paste0(
    "Differential expression analysis was performed using **edgeR** (v", edger_ver,
    ") [Robinson 2010]. ",
    "Counts were normalised using the TMM method. ",
    "Dispersion was estimated using the Cox-Reid adjusted profile likelihood. ",
    "Statistical significance was assessed using the quasi-likelihood F-test."
  ),
  "limma" = paste0(
    "Differential expression analysis was performed using **limma-voom** (v", limma_ver,
    ") [Ritchie 2015]. ",
    "Counts were transformed using voom, which estimates the mean-variance relationship ",
    "and applies precision weights before fitting linear models per gene."
  ),
  paste0("Differential expression analysis was performed using **", de_method, "** (v", de_pkg_ver, ").")
)

p(de_narrative)
blank()

prefilter_str <- if (identical(as.character(prefilter_s), "auto")) {
  paste0(
    "Low-count genes were removed prior to fitting: genes with fewer than **",
    prefilter_n, "** count in at least one sample per group were excluded."
  )
} else {
  paste0(
    "Low-count genes were removed prior to fitting: genes with fewer than **",
    prefilter_n, "** count in at least **", prefilter_s, "** samples were excluded."
  )
}
p(prefilter_str)
blank()

p("**Analysis parameters:**")
blank()
md_table(
  c("Parameter", "Value", "Description"),
  list(
    list("DE method",                  de_method,              "Statistical framework"),
    list("Design formula",             paste0("`", design_f, "`"), "Model formula"),
    list("Group variable",             group_var,              "Metadata column for grouping"),
    list("Normalization strategy",     norm_strat,             "How samples are normalized across comparisons"),
    list("Pre-filter (min counts)",    prefilter_n,            "Minimum count per gene"),
    list("Pre-filter (min samples)",   prefilter_s,            "Minimum number of samples passing count threshold"),
    list("VST blind",                  as.character(vst_blind),"Blind to design when computing VST"),
    list("PCA top genes",              pca_ntop,               "Most-variable genes used for PCA")
  )
)

p("**Significance cutoffs:**")
blank()
md_table(
  c("Cutoff", "Value", "Description"),
  list(
    list("Adjusted p-value (padj)", padj_cut, "Benjamini-Hochberg FDR-corrected Wald test p-value"),
    list("log2 fold change",        lfc_cut,  "Minimum absolute log2(treatment / control)")
  )
)

# ── 3. Functional Enrichment Analysis ──────────────────────────────────────
h(2, "3. Functional Enrichment Analysis")
p(paste0(
  "Over-representation analysis (ORA) was performed on differentially expressed gene sets ",
  "using **clusterProfiler** (v", clusterp_ver, ") [Wu 2021]. ",
  "Gene symbols were mapped to Entrez IDs via **", organism_db, "** (v", orgdb_ver, ") ",
  "through **AnnotationDbi** (v", annodbi_ver, ") [Pagès 2023]. ",
  "GO enrichment was tested across ontologies: Biological Process (BP), Cellular Component (CC), ",
  "and Molecular Function (MF). ",
  "KEGG pathway enrichment was additionally performed for each gene set (KEGG organism code: `", kegg_code, "`)."
))
blank()

md_table(
  c("Parameter", "Value", "Description"),
  list(
    list("GO ontologies",      go_ontologies,       "GO sub-ontologies tested"),
    list("Gene set directions",gene_lists,          "Up-regulated, down-regulated, and total DEGs"),
    list("Organism DB",        organism_db,         "Bioconductor annotation database for ID mapping"),
    list("KEGG organism code", kegg_code,           "KEGG species identifier"),
    list("p-value cutoff",     go_pval,             "Maximum raw p-value for terms to report"),
    list("q-value cutoff",     go_qval,             "Maximum FDR-adjusted p-value (q-value)"),
    list("Min gene set size",  min_gs,              "Minimum annotated genes per GO/KEGG term"),
    list("Max gene set size",  max_gs,              "Maximum annotated genes per GO/KEGG term"),
    list("Min DE gene count",  min_gene_count,      "Minimum overlap between DEGs and term"),
    list("Plot top N",         plot_top_n,          "Top terms shown in dot/bar plots")
  )
)

# ── 4. Software Versions ───────────────────────────────────────────────────
h(2, "4. Software Versions")

sw_rows <- list(
  list("R",              r_ver,         "Statistical computing environment"),
  list(de_method,        de_pkg_ver,    "Differential expression testing")
)
# Include other DE packages if installed (may be used for normalization even if not primary)
if (de_method != "DESeq2" && deseq2_ver != "N/A")
  sw_rows <- c(sw_rows, list(list("DESeq2",  deseq2_ver,  "Normalization / VST transformation")))
if (de_method != "edgeR"  && edger_ver  != "N/A")
  sw_rows <- c(sw_rows, list(list("edgeR",   edger_ver,   "Normalization support")))
if (de_method != "limma"  && limma_ver  != "N/A")
  sw_rows <- c(sw_rows, list(list("limma",   limma_ver,   "Normalization support")))
sw_rows <- c(sw_rows, list(
  list("clusterProfiler", clusterp_ver,  "GO and KEGG enrichment analysis"),
  list("AnnotationDbi",   annodbi_ver,   "Bioconductor annotation interface"),
  list(organism_db,       orgdb_ver,     paste0("Gene annotation (", species_cap, ")")),
  list("enrichplot",      enrichplot_ver,"Enrichment result visualization"),
  list("ggplot2",         ggplot2_ver,   "Data visualization"),
  list("Snakemake",       snakemake_ver, "Workflow management")
))

md_table(c("Tool / Package", "Version", "Purpose"), sw_rows)
p(paste0(
  "All R packages were managed via conda (environment: `rna-seq-de-go-analysis`). ",
  "The full environment specification is available at `environment.yml` in the pipeline root directory."
))
blank()

# ── 5. References ──────────────────────────────────────────────────────────
h(2, "5. References")
li(paste0("Love MI, et al. (2014). Moderated estimation of fold change and dispersion for ",
          "RNA-seq data with DESeq2. *Genome Biology*, 15:550."))
li(paste0("Robinson MD, et al. (2010). edgeR: a Bioconductor package for differential expression analysis ",
          "of digital gene expression data. *Bioinformatics*, 26(1):139–140."))
li(paste0("Ritchie ME, et al. (2015). limma powers differential expression analyses for RNA-sequencing ",
          "and microarray studies. *Nucleic Acids Research*, 43(7):e47."))
li(paste0("Wu T, et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. ",
          "*Innovation*, 2(3):100141."))
li(paste0("Zhu A, et al. (2019). Heavy-tailed prior distributions for sequence count data: removing the noise ",
          "and preserving large differences. *Bioinformatics*, 35(12):2084–2092."))
li(paste0("Pagès H, et al. (2023). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. ",
          "R package version ", annodbi_ver, "."))
li(paste0("Köster J & Rahmann S. (2012). Snakemake — a scalable bioinformatics workflow engine. ",
          "*Bioinformatics*, 28(19):2520–2522."))
blank()

# ── Appendix: Full Configuration ───────────────────────────────────────────
h(2, "Appendix: Full Configuration Parameters")
p(paste0(
  "The following parameters were specified in the project configuration file ",
  "(`", basename(opt$config), "`). This appendix provides a complete record for reproducibility."
))
blank()

flatten_yaml <- function(x, prefix = "") {
  rows <- list()
  for (k in names(x)) {
    key <- if (nchar(prefix) > 0) paste0(prefix, ".", k) else k
    v   <- x[[k]]
    if (is.list(v) && !is.null(names(v))) {
      rows <- c(rows, flatten_yaml(v, key))
    } else {
      display <- paste(as.character(v), collapse = ", ")
      if (nchar(display) > 90) display <- paste0(substr(display, 1, 87), "...")
      rows <- c(rows, list(list(paste0("`", key, "`"), paste0("`", display, "`"))))
    }
  }
  rows
}

md_table(c("Parameter", "Value"), flatten_yaml(cfg))

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------

output_path <- opt$output
output_dir  <- dirname(output_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

writeLines(lines, output_path)
cat(paste0("[INFO] Methods section written to: ", output_path, "\n"))
