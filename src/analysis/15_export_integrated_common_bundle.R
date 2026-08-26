# 파일 경로: src/analysis/15_export_integrated_common_bundle.R
# fig-atlas 그림 렌더링 파이프라인이 바로 소비할 수 있는 "번들" 폴더로, 프로젝트 전체
# (여러 pairwise 비교조건을 가로지르는) 결과를 재포장한다
# (계약서: integrated_bundle_contract.md §2/§3 — count_summary/venn/common DEG heatmap).
#
# 사용법: Rscript 15_export_integrated_common_bundle.R [config_path]

suppressPackageStartupMessages({
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 & config 로드 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 15_export_integrated_common_bundle.R [config_path]")
}
config_path <- args[1]
config <- yaml.load_file(config_path)

fb_cfg <- config$export$fig_atlas_bundle
if (is.null(fb_cfg) || !isTRUE(fb_cfg$enabled)) {
  cat("[15_export_integrated_common_bundle] export.fig_atlas_bundle.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

pairs_cfg <- config$de_analysis$pairwise_comparisons %||% list()
if (length(pairs_cfg) < 2) {
  cat("[15_export_integrated_common_bundle] pairwise_comparisons < 2 — skipping (count_summary/venn/common DEG need >=2 comparisons).\n")
  quit(save = "no", status = 0)
}

bundle_root <- file.path(config$output_dir, "fig_bundles")
dir.create(bundle_root, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("[15_export_integrated_common_bundle] -> %s\n", bundle_root))

padj_cut <- config$de_analysis$padj_cutoff %||% 0.05
lfc_cut  <- config$de_analysis$log2fc_cutoff %||% 1.0

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
# 계약서: sample_columns/sample_groups/group_colors는 파이썬 literal을 문자열로 담은
# 형태(ast.literal_eval로 파싱) — 네이티브 YAML 시퀀스/맵으로 쓰면 안 됨.
py_list_str <- function(x) paste0("[", paste(sprintf("'%s'", x), collapse = ", "), "]")
py_dict_of_list_str <- function(named_list) {
  entries <- sprintf("'%s': %s", names(named_list), vapply(named_list, py_list_str, character(1)))
  paste0("{", paste(entries, collapse = ", "), "}")
}
py_dict_of_str_str <- function(named_vec) {
  entries <- sprintf("'%s': '%s'", names(named_vec), named_vec)
  paste0("{", paste(entries, collapse = ", "), "}")
}

# ─────────────────────────────────────────────────────────────
# 1) count_summary_bundle/ + venn_diagram_bundle/ (+ common DEG 교집합용 집합 수집)
# ─────────────────────────────────────────────────────────────
count_rows    <- list()
venn_rows     <- list()
sig_gene_sets <- list()

for (p in pairs_cfg) {
  compare <- p[[1]]; base <- p[[2]]
  dataset_label <- paste(compare, "vs", base)
  de_path <- file.path(config$output_dir, "pairwise", paste0(compare, "_vs_", base), "final_de_results.csv")
  if (!file.exists(de_path)) {
    cat(sprintf("  [common] %s not found — skipped.\n", de_path))
    next
  }
  de <- read.csv(de_path, row.names = 1, check.names = FALSE)
  sig <- de[!is.na(de$padj) & de$padj < padj_cut &
            !is.na(de$log2FoldChange) & abs(de$log2FoldChange) >= lfc_cut, ]
  if (nrow(sig) == 0) next

  count_rows[[dataset_label]] <- data.frame(
    dataset = dataset_label, log2fc = sig$log2FoldChange, adj_pvalue = sig$padj,
    stringsAsFactors = FALSE
  )
  venn_rows[[dataset_label]] <- data.frame(
    dataset = dataset_label, item = rownames(sig), stringsAsFactors = FALSE
  )
  sig_gene_sets[[dataset_label]] <- rownames(sig)
}

if (length(count_rows) > 0) {
  write_bundle_csv(do.call(rbind, count_rows), "count_summary_bundle")
}

if (length(venn_rows) >= 2 && length(venn_rows) <= 3) {
  write_bundle_csv(do.call(rbind, venn_rows), "venn_diagram_bundle")
} else if (length(venn_rows) > 3) {
  cat(sprintf("  [venn_diagram_bundle] %d comparisons (>3) — skipped (venn renderer targets 2-3 sets).\n", length(venn_rows)))
}

# ─────────────────────────────────────────────────────────────
# 2) Filtered_common_deg_genes_heatmap_bundle/ — 모든 pairwise 비교의 교집합
# ─────────────────────────────────────────────────────────────
if (length(sig_gene_sets) == length(pairs_cfg) && length(sig_gene_sets) >= 2) {
  common_genes <- Reduce(intersect, sig_gene_sets)
} else {
  common_genes <- character(0)
  cat("  [Filtered_common_deg_genes_heatmap_bundle] one or more comparisons had 0 significant genes or failed to load — no common DEG.\n")
}

mg_path <- file.path(config$output_dir, "multi_group_result.csv")
if (length(common_genes) == 0) {
  cat("  [Filtered_common_deg_genes_heatmap_bundle] no common DEG across all comparisons — skipped.\n")
} else if (!file.exists(mg_path)) {
  cat("  [Filtered_common_deg_genes_heatmap_bundle] multi_group_result.csv not found (de_analysis.multi_group_export.enabled false?) — skipped.\n")
} else {
  mg <- read.csv(mg_path, row.names = 1, check.names = FALSE)
  mg_common <- mg[rownames(mg) %in% common_genes, , drop = FALSE]

  if (nrow(mg_common) == 0) {
    cat("  [Filtered_common_deg_genes_heatmap_bundle] common DEG not found in multi_group_result.csv — skipped.\n")
  } else {
    stat_cols <- c("baseMean", "stat", "pvalue", "padj")
    if ("gene_symbol" %in% colnames(mg_common)) {
      gene_label  <- mg_common$gene_symbol
      sample_cols <- setdiff(colnames(mg_common), c("gene_symbol", stat_cols))
    } else {
      gene_label  <- rownames(mg_common)
      sample_cols <- setdiff(colnames(mg_common), stat_cols)
    }

    out <- data.frame(gene_label = gene_label, mg_common[, sample_cols, drop = FALSE],
                       check.names = FALSE, stringsAsFactors = FALSE)
    write_bundle_csv(out, "Filtered_common_deg_genes_heatmap_bundle")

    # 샘플 -> 그룹 매핑 (08_generate_methods_section.R의 group-map 패턴 재사용)
    meta <- read.csv(config$metadata_path, row.names = 1)
    group_var <- config$de_analysis$group_variable
    sample_groups_all <- split(rownames(meta), meta[[group_var]])

    # multi_group_result.csv의 샘플 컬럼 등장 순서를 그대로 따라 그룹 순서를 정한다
    # (09_export_multi_group.R이 이미 reference_group 우선 + 알파벳순으로 정렬해둠).
    group_of_sample <- function(s) {
      hit <- names(sample_groups_all)[vapply(sample_groups_all, function(x) s %in% x, logical(1))]
      if (length(hit) == 0) NA_character_ else hit[1]
    }
    group_order <- unique(vapply(sample_cols, group_of_sample, character(1)))
    group_order <- group_order[!is.na(group_order)]
    sample_groups <- sample_groups_all[group_order]
    # 이 번들에 실제로 없는 샘플(다른 비교조건 전용 등)은 그룹 목록에서 제외
    sample_groups <- lapply(sample_groups, function(x) intersect(x, sample_cols))

    palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
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
        title = sprintf("Filtered: Gene List (%d genes) | Z-score | padj<=%.3g, baseMean>=0, n=%d",
                         nrow(out), padj_cut, length(sample_cols))
      )
    ), "Filtered_common_deg_genes_heatmap_bundle")
  }
}

cat("[15_export_integrated_common_bundle] Done.\n")
