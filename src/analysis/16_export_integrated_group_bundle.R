# 파일 경로: src/analysis/16_export_integrated_group_bundle.R
# fig-atlas 그림 렌더링 파이프라인이 바로 소비할 수 있는 "번들" 폴더로, time-series
# 클러스터 / co-expression 모듈 결과를 재포장한다
# (계약서: integrated_bundle_contract.md §1/§4 — 그룹-주석 히트맵 + GO BP dot/KEGG bar).
#
# 11_run_group_enrichment.R과 동일한 인자 스타일에 bundle_prefix 하나를 추가했다.
#
# 사용법: Rscript 16_export_integrated_group_bundle.R [config_path] [input_csv] [group_col] [group_prefix] [bundle_prefix] [output_dir]
#   예) time-series : ... time_series/time_series_significant_genes.csv cluster_id cluster ts_cluster   .../time_series
#       coexpression: ... coexpression_modules/coexpression_module_assignments.csv module_id module coexp_module .../coexpression_modules

suppressPackageStartupMessages({
  library(yaml)
  library(AnnotationDbi)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 & config 로드 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript 16_export_integrated_group_bundle.R [config_path] [input_csv] [group_col] [group_prefix] [bundle_prefix] [output_dir]")
}
config_path   <- args[1]
input_csv     <- args[2]
group_col     <- args[3]   # cluster_id | module_id
group_prefix  <- args[4]   # cluster | module (11_run_group_enrichment.R의 원본 라벨 접두어)
bundle_prefix <- args[5]   # ts_cluster | coexp_module (fig-atlas 번들 폴더 접두어)
output_dir    <- args[6]   # time_series/ 또는 coexpression_modules/ (11이 GO CSV를 쓴 곳)

config <- yaml.load_file(config_path)
fb_cfg <- config$export$fig_atlas_bundle

if (is.null(fb_cfg) || !isTRUE(fb_cfg$enabled)) {
  cat("[16_export_integrated_group_bundle] export.fig_atlas_bundle.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}
if (!file.exists(input_csv)) {
  cat(sprintf("[16_export_integrated_group_bundle] %s not found — skipping.\n", input_csv))
  quit(save = "no", status = 0)
}

# 계약서 폴더 계층: fig_bundles/{time_series,coexpression_modules}/...
bundle_root <- file.path(config$output_dir, "fig_bundles", basename(output_dir))
dir.create(bundle_root, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("[16_export_integrated_group_bundle] %s (group_col=%s) -> %s\n", input_csv, group_col, bundle_root))

# coexp_module만 2자리 zero-pad, ts_cluster는 무패딩 — 계약서 예시(ts_cluster1 vs
# coexp_module01)와 일치하는 폴더명 규칙.
zero_pad_folder <- bundle_prefix == "coexp_module"
folder_suffix <- function(raw_value) if (zero_pad_folder) sprintf("%02d", as.integer(raw_value)) else as.character(raw_value)

# 11_run_group_enrichment.R의 Gene Set 라벨 변환과 동일 — GO 번들의 gene_set 컬럼 값은
# 폴더명 패딩 규칙과 무관하게 항상 이 형식("Cluster01"/"Module02")을 쓴다(계약서 §4).
format_gene_set_label <- function(raw_value, group_prefix) {
  prefix_title <- paste0(toupper(substr(group_prefix, 1, 1)), substr(group_prefix, 2, nchar(group_prefix)))
  sprintf("%s%02d", prefix_title, as.integer(raw_value))
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
py_list_str <- function(x) paste0("[", paste(sprintf("'%s'", x), collapse = ", "), "]")
py_dict_of_list_str <- function(named_list) {
  entries <- sprintf("'%s': %s", names(named_list), vapply(named_list, py_list_str, character(1)))
  paste0("{", paste(entries, collapse = ", "), "}")
}
py_dict_of_str_str <- function(named_vec) {
  entries <- sprintf("'%s': '%s'", names(named_vec), named_vec)
  paste0("{", paste(entries, collapse = ", "), "}")
}

# --- Organism DB (Entrez -> Symbol 변환용) ---
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

# --- 샘플 -> 그룹 매핑 (15번 스크립트와 동일 패턴) ---
meta <- read.csv(config$metadata_path, row.names = 1)
group_var <- config$de_analysis$group_variable
sample_groups_all <- split(rownames(meta), meta[[group_var]])
group_of_sample <- function(s) {
  hit <- names(sample_groups_all)[vapply(sample_groups_all, function(x) s %in% x, logical(1))]
  if (length(hit) == 0) NA_character_ else hit[1]
}
palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# ─────────────────────────────────────────────────────────────
# 그룹(cluster_id/module_id)마다 반복
# ─────────────────────────────────────────────────────────────
d <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
group_values <- sort(unique(d[[group_col]]))

# input_csv에 이미 내장된 gene x sample VST 행렬(01c/10 확인 완료) — 통계/식별자 컬럼만
# 제외하면 나머지가 전부 샘플 컬럼.
known_stat_cols <- c("gene_id", "gene_symbol", "p_value", "r_squared", "padj", group_col)

for (gv in group_values) {
  raw_label   <- paste0(group_prefix, gv)          # 11_run_group_enrichment.R의 원본 파일명 라벨(무패딩)
  suffix      <- folder_suffix(gv)                  # 번들 폴더명용 N/NN
  gene_set_lb <- format_gene_set_label(gv, group_prefix)  # "Cluster01"/"Module02"

  # --- 1) Filtered_{bundle_prefix}{N|NN}_genes_heatmap_bundle/ ---
  sub_d <- d[d[[group_col]] == gv, , drop = FALSE]
  if (nrow(sub_d) > 0) {
    gene_label  <- if ("gene_symbol" %in% colnames(sub_d)) sub_d$gene_symbol else sub_d$gene_id
    sample_cols <- setdiff(colnames(sub_d), known_stat_cols)

    out <- data.frame(gene_label = gene_label, sub_d[, sample_cols, drop = FALSE],
                       check.names = FALSE, stringsAsFactors = FALSE)
    heatmap_bundle <- sprintf("Filtered_%s%s_genes_heatmap_bundle", bundle_prefix, suffix)
    write_bundle_csv(out, heatmap_bundle)

    group_order <- unique(vapply(sample_cols, group_of_sample, character(1)))
    group_order <- group_order[!is.na(group_order)]
    sample_groups <- sample_groups_all[group_order]
    sample_groups <- lapply(sample_groups, function(x) intersect(x, sample_cols))
    group_colors <- setNames(palette[((seq_along(sample_groups) - 1) %% length(palette)) + 1], names(sample_groups))

    write_bundle_yaml(list(
      plot_type = "multi_group_heatmap",
      plot_params = list(
        gene_label_col = "gene_label",
        sample_columns = py_list_str(sample_cols),
        sample_groups  = py_dict_of_list_str(sample_groups),
        group_colors   = py_dict_of_str_str(group_colors),
        cmap = "bwr", linkage = "ward", metric = "euclidean",
        cluster_rows = FALSE, cluster_cols = FALSE,
        z_auto = TRUE, z_min = -2.0, z_max = 2.0,
        show_gene_labels = FALSE, gene_fontsize = 7, show_col_labels = TRUE,
        title = sprintf("Filtered: %s (%d genes) | Z-score", gene_set_lb, nrow(out))
      )
    ), heatmap_bundle)
  }

  # --- 2) go_bp_{group_prefix}{N|NN}_dot_bundle/ — 계약서 1 §2와 동일 스키마 ---
  tc_path <- file.path(output_dir, sprintf("go_termcluster_%s_BP.csv", raw_label))
  if (file.exists(tc_path)) {
    td <- read.csv(tc_path, stringsAsFactors = FALSE)
    if (nrow(td) > 0) {
      cluster_size <- table(td$cluster)
      fe_median    <- tapply(td$FoldEnrichment, td$cluster, median, na.rm = TRUE)

      go_out <- data.frame(
        cluster_id      = sprintf("%03d", td$cluster),
        gene_set        = gene_set_lb,
        ontology        = "BP",
        term_id         = td$ID,
        description     = td$Description,
        gene_ratio      = td$GeneRatio,
        bg_ratio        = td$BgRatio,
        pvalue          = td$pvalue,
        fdr             = td$p.adjust,
        qvalue          = td$qvalue,
        gene_count      = td$Count,
        gene_symbols    = convert_entrez_column_to_symbols(td$geneID, organism_db),
        fold_enrichment = td$FoldEnrichment,
        direction       = gene_set_lb,
        check.names = FALSE, stringsAsFactors = FALSE
      )
      go_out[["_cluster_size"]] <- as.integer(cluster_size[as.character(td$cluster)])
      go_out[["_fe_median"]]    <- as.numeric(fe_median[as.character(td$cluster)])

      dot_bundle <- sprintf("go_bp_%s%s_dot_bundle", group_prefix, suffix)
      write_bundle_csv(go_out, dot_bundle)
      write_bundle_yaml(list(
        plot_type = "go_cluster_dot",
        plot_params = list(
          x_axis = "-log10(FDR)", color_by = "Fold Enrichment", size_by = "Cluster size",
          sort_by = "FDR", top_n = 30, top_n_by = "FDR",
          dot_size_min = 40.0, dot_size_scale = 20.0, cmap = "YlOrRd_r",
          colorbar_loc = "right", colorbar_shrink = 0.6, show_size_legend = TRUE
        )
      ), dot_bundle)
    }
  }

  # --- 3) go_kegg_{group_prefix}{N|NN}_bar_chart_bundle/ — 계약서 1 §3와 동일 스키마 ---
  kegg_path <- file.path(output_dir, sprintf("kegg_enrichment_%s.csv", raw_label))
  if (file.exists(kegg_path)) {
    kd <- read.csv(kegg_path, stringsAsFactors = FALSE)
    if (nrow(kd) > 0) {
      kegg_out <- data.frame(
        gene_set        = gene_set_lb,
        ontology        = "KEGG",
        term_id         = kd$ID,
        description     = kd$Description,
        gene_ratio      = kd$GeneRatio,
        bg_ratio        = kd$BgRatio,
        pvalue          = kd$pvalue,
        fdr             = kd$p.adjust,
        qvalue          = kd$qvalue,
        gene_count      = kd$Count,
        gene_symbols    = convert_entrez_column_to_symbols(kd$geneID, organism_db),
        fold_enrichment = compute_fold_enrichment(kd),
        direction       = gene_set_lb,
        stringsAsFactors = FALSE
      )
      bar_bundle <- sprintf("go_kegg_%s%s_bar_chart_bundle", group_prefix, suffix)
      write_bundle_csv(kegg_out, bar_bundle)
      write_bundle_yaml(list(
        plot_type = "go_bar",
        plot_params = list(
          top_n = 30, x_axis = "-log10(FDR)", sort_by = "FDR (ascending)",
          bar_color = "#4682b4", horizontal = TRUE, xlabel_text = "-log10(FDR)"
        )
      ), bar_bundle)
    }
  }
}

cat("[16_export_integrated_group_bundle] Done.\n")
