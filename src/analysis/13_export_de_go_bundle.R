# 파일 경로: src/analysis/13_export_de_go_bundle.R
# fig-atlas 그림 렌더링 파이프라인이 바로 소비할 수 있는 "번들" 폴더
# (inputs/data.csv + metadata/metadata.yaml)로 pairwise DE/GO/KEGG/rrvgo 결과를
# 재포장한다 (계약서: de_go_bundle_contract.md, 사용자 제공 — CMG-SeqViewer에서 매번
# 같은 파라미터로 반복 export하던 작업을 대체).
#
# 대부분의 번들은 이미 있는 CSV를 컬럼명만 바꾸거나(rrvgo는 변환조차 없이) 그대로
# 넣는 수준이라, 새로 계산하는 값은 fold_enrichment(KEGG bar chart)와
# _cluster_size/_fe_median(GO cluster dot) 뿐이다.
#
# 사용법: Rscript 13_export_de_go_bundle.R [config_path] [compare_group] [base_group] [pair_output_dir]

suppressPackageStartupMessages({
  library(yaml)
  library(openxlsx)
  library(AnnotationDbi)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 & config 로드 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 13_export_de_go_bundle.R [config_path] [compare_group] [base_group] [pair_output_dir]")
}
config_path     <- args[1]
compare_group   <- args[2]
base_group      <- args[3]
pair_output_dir <- args[4]

config <- yaml.load_file(config_path)
fb_cfg <- config$export$fig_atlas_bundle

if (is.null(fb_cfg) || !isTRUE(fb_cfg$enabled)) {
  cat("[13_export_de_go_bundle] export.fig_atlas_bundle.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

bundle_root <- file.path(config$output_dir, "fig_bundles", paste0(compare_group, "_vs_", base_group))
dir.create(bundle_root, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[13_export_de_go_bundle] %s vs %s -> %s\n", compare_group, base_group, bundle_root))

# --- 2. Organism DB (Entrez -> Symbol 변환용, 실패해도 export 자체는 계속) ---
organism_db <- NULL
species_info <- config$databases[[config$species]]
if (!is.null(species_info$organism_db) &&
    require(species_info$organism_db, character.only = TRUE, quietly = TRUE)) {
  organism_db <- get(species_info$organism_db)
}

convert_entrez_column_to_symbols <- function(entrez_ids_strings, organism_db) {
  if (is.null(organism_db)) return(entrez_ids_strings)
  non_empty <- entrez_ids_strings[!is.na(entrez_ids_strings) & entrez_ids_strings != ""]
  all_ids <- unique(unlist(strsplit(as.character(non_empty), "/")))
  lookup <- tryCatch({
    mapIds(organism_db, keys = all_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  }, error = function(e) setNames(all_ids, all_ids))
  vapply(entrez_ids_strings, function(x) {
    if (is.na(x) || x == "") return("")
    ids <- strsplit(as.character(x), "/")[[1]]
    symbols <- lookup[ids]
    symbols[is.na(symbols)] <- ids[is.na(symbols)]
    paste(symbols, collapse = "/")
  }, character(1), USE.NAMES = FALSE)
}

parse_ratio <- function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1] / parts[2]
}
compute_fold_enrichment <- function(df) {
  mapply(function(gr, br) parse_ratio(gr) / parse_ratio(br), df$GeneRatio, df$BgRatio)
}

write_bundle_csv <- function(df, bundle_name) {
  dir <- file.path(bundle_root, bundle_name, "inputs")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(dir, "data.csv"), row.names = FALSE)
  cat(sprintf("  [%s] %d rows\n", bundle_name, nrow(df)))
}

write_bundle_yaml <- function(params_list, bundle_name) {
  dir <- file.path(bundle_root, bundle_name, "metadata")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write_yaml(params_list, file.path(dir, "metadata.yaml"))
}

# ─────────────────────────────────────────────────────────────
# 1) volcano_plot_bundle/ — final_de_results.xlsx의 DE_Results 시트가 이미
#    gene_id/symbol/baseMean/log2FoldChange/lfcSE/stat/pvalue/padj + 샘플 컬럼까지
#    한 시트에 합쳐져 있어(01b_run_pairwise_de.R) 컬럼명만 바꾸면 된다.
# ─────────────────────────────────────────────────────────────
de_xlsx <- file.path(pair_output_dir, "final_de_results.xlsx")
if (file.exists(de_xlsx)) {
  de <- read.xlsx(de_xlsx, sheet = "DE_Results", check.names = FALSE, sep.names = " ")

  rename_map <- c(baseMean = "base_mean", log2FoldChange = "log2FC", lfcSE = "lfcse")
  for (old in names(rename_map)) {
    if (old %in% colnames(de)) colnames(de)[colnames(de) == old] <- rename_map[[old]]
  }
  write_bundle_csv(de, "volcano_plot_bundle")

  padj_cut    <- config$de_analysis$padj_cutoff %||% 0.05
  lfc_cut     <- config$de_analysis$log2fc_cutoff %||% 1.0
  volcano_cfg <- config$plot_aesthetics$volcano %||% list()

  x_range <- suppressWarnings(range(de$log2FC, na.rm = TRUE))
  y_max   <- suppressWarnings(max(-log10(de$padj), na.rm = TRUE))
  if (!all(is.finite(x_range))) x_range <- c(-1, 1)
  if (!is.finite(y_max)) y_max <- 1

  write_bundle_yaml(list(
    figure_id = "volcano_plot",
    plot_type = "volcano",
    plot_params = list(
      padj_threshold = padj_cut,
      log2fc_threshold = lfc_cut,
      down_color = volcano_cfg$down_color %||% "#0000ff",
      up_color = volcano_cfg$up_color %||% "#ff0000",
      ns_color = "#808080",
      dot_size = 20,
      x_min = floor(x_range[1] * 10) / 10,
      x_max = ceiling(x_range[2] * 10) / 10,
      y_min = 0.0,
      y_max = round(y_max + 0.5, 1),
      annotation_mode = "top_n",
      annotation_top_n = 30,
      annotation_label_size = 8,
      labels_title = "Volcano Plot",
      labels_xlabel = "Log2 Fold Change",
      labels_ylabel = "-Log10(Padj)"
    )
  ), "volcano_plot_bundle")
} else {
  cat("  [volcano_plot_bundle] final_de_results.xlsx not found — skipped.\n")
}

# ─────────────────────────────────────────────────────────────
# 2) go_bp_{up,down}_cluster_dot_bundle/ — go_termcluster_{gs}_BP.csv
#    (03_enrichment_analysis.R이 유의 term>=3일 때만 생성 -> 파일 존재 여부가 곧
#    "해당 방향에 유의 BP term이 있을 때만"이라는 계약서 규칙과 일치)
# ─────────────────────────────────────────────────────────────
for (gs in c("up", "down")) {
  tc_path <- file.path(pair_output_dir, "enrichment", sprintf("go_termcluster_%s_BP.csv", gs))
  if (!file.exists(tc_path)) next
  d <- read.csv(tc_path, stringsAsFactors = FALSE)
  if (nrow(d) == 0) next

  cluster_size <- table(d$cluster)
  fe_median    <- tapply(d$FoldEnrichment, d$cluster, median, na.rm = TRUE)

  out <- data.frame(
    cluster_id      = sprintf("%03d", d$cluster),
    gene_set        = toupper(gs),
    ontology        = "BP",
    term_id         = d$ID,
    description     = d$Description,
    gene_ratio      = d$GeneRatio,
    bg_ratio        = d$BgRatio,
    pvalue          = d$pvalue,
    fdr             = d$p.adjust,
    qvalue          = d$qvalue,
    gene_count      = d$Count,
    gene_symbols    = convert_entrez_column_to_symbols(d$geneID, organism_db),
    fold_enrichment = d$FoldEnrichment,
    direction       = toupper(gs),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  out[["_cluster_size"]] <- as.integer(cluster_size[as.character(d$cluster)])
  out[["_fe_median"]]    <- as.numeric(fe_median[as.character(d$cluster)])

  bundle_name <- sprintf("go_bp_%s_cluster_dot_bundle", gs)
  write_bundle_csv(out, bundle_name)
  write_bundle_yaml(list(
    plot_type = "go_cluster_dot",
    plot_params = list(
      x_axis = "-log10(FDR)", color_by = "Fold Enrichment", size_by = "Cluster size",
      sort_by = "FDR", top_n = 30, top_n_by = "FDR",
      dot_size_min = 40.0, dot_size_scale = 20.0, cmap = "YlOrRd_r",
      colorbar_loc = "right", colorbar_shrink = 0.6, show_size_legend = TRUE
    )
  ), bundle_name)
}

# ─────────────────────────────────────────────────────────────
# 3) go_kegg_{up,down}_bar_chart_bundle/ — kegg_enrichment_{gs}.csv (raw enrichResult,
#    FoldEnrichment 없음 -> compute_fold_enrichment()로 새로 계산)
# ─────────────────────────────────────────────────────────────
for (gs in c("up", "down")) {
  kegg_path <- file.path(pair_output_dir, "enrichment", sprintf("kegg_enrichment_%s.csv", gs))
  if (!file.exists(kegg_path)) next
  d <- read.csv(kegg_path, stringsAsFactors = FALSE)
  if (nrow(d) == 0) next

  out <- data.frame(
    gene_set        = toupper(gs),
    ontology        = "KEGG",
    term_id         = d$ID,
    description     = d$Description,
    gene_ratio      = d$GeneRatio,
    bg_ratio        = d$BgRatio,
    pvalue          = d$pvalue,
    fdr             = d$p.adjust,
    qvalue          = d$qvalue,
    gene_count      = d$Count,
    gene_symbols    = convert_entrez_column_to_symbols(d$geneID, organism_db),
    fold_enrichment = compute_fold_enrichment(d),
    direction       = toupper(gs),
    stringsAsFactors = FALSE
  )

  bundle_name <- sprintf("go_kegg_%s_bar_chart_bundle", gs)
  write_bundle_csv(out, bundle_name)
  write_bundle_yaml(list(
    plot_type = "go_bar",
    plot_params = list(
      top_n = 30, x_axis = "-log10(FDR)", sort_by = "FDR (ascending)",
      bar_color = "#4682b4", horizontal = TRUE, xlabel_text = "-log10(FDR)"
    )
  ), bundle_name)
}

# ─────────────────────────────────────────────────────────────
# 4) go_{bp,cc,mf}_{up,down}_rrvgo_treemap/ — go_rrvgo_{gs}_{ont}.csv
#    (rrvgo::reduceSimMatrix() 원본 컬럼이 계약서 스펙과 이미 정확히 일치 -> 변환 없음)
# ─────────────────────────────────────────────────────────────
for (ont in c("BP", "CC", "MF")) {
  for (gs in c("up", "down")) {
    rr_path <- file.path(pair_output_dir, "enrichment", sprintf("go_rrvgo_%s_%s.csv", gs, ont))
    if (!file.exists(rr_path)) next
    d <- read.csv(rr_path, stringsAsFactors = FALSE)
    if (nrow(d) == 0) next
    write_bundle_csv(d, sprintf("go_%s_%s_rrvgo_treemap", tolower(ont), gs))
  }
}

cat("[13_export_de_go_bundle] Done.\n")
