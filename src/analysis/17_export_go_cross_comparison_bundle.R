# 파일 경로: src/analysis/17_export_go_cross_comparison_bundle.R
#
# 계약서: go_cross_comparison_bundle_contract.md (사용자 제공, 외부 fig-atlas 레포
# raw_data/H2O2-human-mouse-mrna-atac-seq/... 로 최종 전달될 예정).
# 이 레포 자체는 "output/ 안에만 생성" 원칙(이전 계약서들과 동일)을 따르므로, 외부
# raw_data 경로에 직접 쓰지 않고 이 프로젝트의 output_dir 아래 계약서와 동일한 폴더
# 구조(go-cross-comparison/)로 생성한다. 최종 외부 레포 복사는 사용자가 직접 처리.
#
# 12_run_cross_condition_comparison.R이 이미 계약서와 동일한 컬럼 스키마로
# cross_condition_plot_data_semantic_{ont}.csv / cross_condition_plot_data_{ont}.csv를
# 생성해두므로 변환은 최소한(Condition 값의 base-group 라벨을 계약서가 요구하는 값으로
# 치환)만 필요 — 값 재계산이나 재필터링은 하지 않는다.
#
# 사용법: Rscript 17_export_go_cross_comparison_bundle.R <config_path> <cross_condition_dir> <base_group_label> <target_base_label>
#   예) Rscript 17_export_go_cross_comparison_bundle.R \
#         configs/config_2026-06-hiy-human-rna.yml \
#         output/2026-06-hiy-human-rna/cross_condition \
#         Cont Control

suppressPackageStartupMessages({
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 17_export_go_cross_comparison_bundle.R <config_path> <cross_condition_dir> <base_group_label> <target_base_label>")
}
config_path       <- args[1]
cc_dir            <- args[2]
base_group_label  <- args[3]
target_base_label <- args[4]

config <- yaml.load_file(config_path)
bundle_root <- file.path(config$output_dir, "go-cross-comparison")
dir.create(bundle_root, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[17_export_go_cross_comparison_bundle] source: %s\n", cc_dir))
cat(sprintf("[17_export_go_cross_comparison_bundle] target: %s\n", bundle_root))
cat(sprintf("[17_export_go_cross_comparison_bundle] Condition label: '_vs_%s' -> '_vs_%s'\n\n",
            base_group_label, target_base_label))

relabel_condition <- function(x) {
  suffix_old <- paste0("_vs_", base_group_label)
  suffix_new <- paste0("_vs_", target_base_label)
  if (!all(endsWith(x, suffix_old))) {
    bad <- unique(x[!endsWith(x, suffix_old)])
    stop("Condition 값 중 '_vs_", base_group_label, "'로 끝나지 않는 값이 있습니다: ",
         paste(bad, collapse = ", "))
  }
  sub(paste0(suffix_old, "$"), suffix_new, x)
}

# 계약서 §3a row-level invariant 검증
validate_bundle <- function(d, is_kegg) {
  problems <- character(0)

  dup_key <- paste(d[["GO ID"]], d[["Condition"]], sep = "|||")
  if (any(duplicated(dup_key))) {
    problems <- c(problems, sprintf("duplicate (GO ID, Condition) pairs: %d", sum(duplicated(dup_key))))
  }

  expected_sign <- ifelse(d[["Direction"]] == "UP", 1, -1)
  actual_sign <- sign(d[["Signed -log10(FDR)"]])
  if (!all(actual_sign == expected_sign)) {
    problems <- c(problems, sprintf("sign(Signed -log10(FDR)) mismatch in %d row(s)", sum(actual_sign != expected_sign)))
  }

  expected_abs <- -log10(d[["Adjusted P-value"]])
  if (!all(abs(abs(d[["Signed -log10(FDR)"]]) - expected_abs) < 1e-6)) {
    problems <- c(problems, "abs(Signed -log10(FDR)) does not match -log10(Adjusted P-value) within tolerance")
  }

  if (!all(d[["Adjusted P-value"]] <= 0.05)) {
    problems <- c(problems, sprintf("%d row(s) with Adjusted P-value > 0.05", sum(d[["Adjusted P-value"]] > 0.05)))
  }

  if (is_kegg) {
    prefix <- substr(tolower(config$databases[[config$species]]$kegg_code %||% ""), 1, 3)
    if (nzchar(prefix) && !all(grepl(paste0("^", prefix), tolower(d[["GO ID"]])))) {
      problems <- c(problems, sprintf("GO ID(=KEGG pathway ID) values not all prefixed with '%s'", prefix))
    }
  }

  # Direction Flip 재계산 후 기존 값과 일치하는지 확인 (재라벨링은 GO ID별 flip 여부에
  # 영향 없어야 함)
  flip_recalc <- ave(d[["Direction"]], d[["GO ID"]], FUN = function(x) length(unique(x)) > 1)
  flip_recalc <- as.logical(flip_recalc)
  if (!all(flip_recalc == d[["Direction Flip"]])) {
    problems <- c(problems, "Direction Flip recomputation mismatch")
  }

  problems
}

`%||%` <- function(x, y) if (is.null(x)) y else x

write_bundle <- function(src_path, bundle_name, is_kegg) {
  if (!file.exists(src_path)) {
    cat(sprintf("  [%s] source not found (%s) — skipped.\n", bundle_name, src_path))
    return(invisible(NULL))
  }
  d <- read.csv(src_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(d) == 0) {
    cat(sprintf("  [%s] source has 0 rows — skipped.\n", bundle_name))
    return(invisible(NULL))
  }
  d[["Condition"]] <- relabel_condition(d[["Condition"]])

  problems <- validate_bundle(d, is_kegg)
  if (length(problems) > 0) {
    stop(sprintf("[%s] contract invariant violation(s):\n  - %s", bundle_name, paste(problems, collapse = "\n  - ")))
  }

  out_dir <- file.path(bundle_root, bundle_name, "inputs")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(d, file.path(out_dir, "data.csv"), row.names = FALSE)

  n_terms <- length(unique(d[["GO ID"]]))
  n_conditions <- length(unique(d[["Condition"]]))
  cat(sprintf("  [%s] OK — %d rows, %d unique terms, %d conditions -> %s\n",
              bundle_name, nrow(d), n_terms, n_conditions, file.path(out_dir, "data.csv")))
}

write_bundle(file.path(cc_dir, "cross_condition_plot_data_semantic_BP.csv"),
             "cross_condition_heatmap_semantic-bp", is_kegg = FALSE)
write_bundle(file.path(cc_dir, "cross_condition_plot_data_KEGG.csv"),
             "cross_condition_heatmap_semantic-kegg", is_kegg = TRUE)

cat("\n[17_export_go_cross_comparison_bundle] Done.\n")
