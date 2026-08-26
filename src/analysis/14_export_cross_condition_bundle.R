# 파일 경로: src/analysis/14_export_cross_condition_bundle.R
# fig-atlas 그림 렌더링 파이프라인이 바로 소비할 수 있는 "번들" 폴더로 cross-condition
# GO/KEGG 비교 결과를 재포장한다 (계약서: cross_condition_bundle_contract.md, 사용자 제공).
#
# heatmap/dotplot semantic 번들의 소스 표(cross_condition_plot_data_semantic_{ont}.csv)와
# upset 번들의 소스 표(upset_membership_{DIR}_{ont}.csv)는 12_run_cross_condition_comparison.R이
# 이미 계약서 컬럼과 정확히 같은 스키마로 만들어두므로 변환 없이 복사만 한다.
#
# 사용법: Rscript 14_export_cross_condition_bundle.R [config_path] [output_dir]
#   output_dir는 12_run_cross_condition_comparison.R과 동일하게 cross_condition/ 폴더

suppressPackageStartupMessages({
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 & config 로드 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 14_export_cross_condition_bundle.R [config_path] [output_dir]")
}
config_path <- args[1]
output_dir  <- args[2]

config <- yaml.load_file(config_path)
fb_cfg <- config$export$fig_atlas_bundle

if (is.null(fb_cfg) || !isTRUE(fb_cfg$enabled)) {
  cat("[14_export_cross_condition_bundle] export.fig_atlas_bundle.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

cc_cfg <- config$enrichment$cross_condition
if (is.null(cc_cfg) || !isTRUE(cc_cfg$enabled)) {
  cat("[14_export_cross_condition_bundle] enrichment.cross_condition.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

go_ontologies <- cc_cfg$go_ontologies %||% c("BP")
bundle_root <- file.path(config$output_dir, "fig_bundles", "cross_condition")
dir.create(bundle_root, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[14_export_cross_condition_bundle] -> %s\n", bundle_root))

write_bundle_csv <- function(df, bundle_name) {
  dir <- file.path(bundle_root, bundle_name, "inputs")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(dir, "data.csv"), row.names = FALSE)
  cat(sprintf("  [%s] %d rows\n", bundle_name, nrow(df)))
}

# ─────────────────────────────────────────────────────────────
# 1) cross_condition_heatmap_semantic-{ont}/ + cross_condition_dotplot_semantic-{ont}/
#    소스 표는 heatmap/dotplot 두 뷰가 완전히 동일(같은 원본을 다르게 시각화) — 두
#    번들 폴더에 동일 파일을 각각 복사.
#    KEGG는 12_run_cross_condition_comparison.R이 rrvgo semantic reduction을 원래
#    건너뛰므로(term 수가 적어 축약 실익이 낮다는 기존 설계) semantic 버전이 없다 —
#    이 경우 raw candidate 표(cross_condition_plot_data_{ont}.csv)를 대신 쓴다(KEGG는
#    애초에 추가 축약이 없어 raw candidate 표가 사실상 최종 표이기 때문).
# ─────────────────────────────────────────────────────────────
for (ont in go_ontologies) {
  is_kegg <- toupper(ont) == "KEGG"
  semantic_path <- file.path(output_dir, sprintf("cross_condition_plot_data_semantic_%s.csv", ont))
  raw_path      <- file.path(output_dir, sprintf("cross_condition_plot_data_%s.csv", ont))
  src_path <- if (is_kegg) raw_path else semantic_path

  if (!file.exists(src_path)) {
    cat(sprintf("  [cross_condition_*_semantic-%s] source not found (%s) — skipped.\n",
                tolower(ont), basename(src_path)))
    next
  }
  d <- read.csv(src_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(d) == 0) next

  ont_lower <- tolower(ont)
  write_bundle_csv(d, sprintf("cross_condition_heatmap_semantic-%s", ont_lower))
  write_bundle_csv(d, sprintf("cross_condition_dotplot_semantic-%s", ont_lower))
}

# ─────────────────────────────────────────────────────────────
# 2) upset-{ont}-{up,down}/ — upset_membership_{DIR}_{ont}.csv, 컬럼 이미 정확히 일치.
#    N Conditions>=2인 행이 하나도 없으면(구조적으로 겹치는 term이 없으면) 폴더 자체를
#    생성하지 않는다(계약서 규칙) — 행 자체는 필터링하지 않고 폴더 생성 여부만 이 기준.
# ─────────────────────────────────────────────────────────────
for (ont in go_ontologies) {
  for (dir_lab in c("UP", "DOWN")) {
    membership_path <- file.path(output_dir, sprintf("upset_membership_%s_%s.csv", dir_lab, ont))
    if (!file.exists(membership_path)) next
    d <- read.csv(membership_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(d) == 0) next
    if (!any(d[["N Conditions"]] >= 2, na.rm = TRUE)) {
      cat(sprintf("  [upset-%s-%s] no term with N Conditions>=2 — skipped.\n", tolower(ont), tolower(dir_lab)))
      next
    }
    write_bundle_csv(d, sprintf("upset-%s-%s", tolower(ont), tolower(dir_lab)))
  }
}

cat("[14_export_cross_condition_bundle] Done.\n")
