# 파일 경로: src/analysis/01c_run_masigpro_timeseries.R
# Time-series DE 분석 (maSigPro, count 기반 negative binomial 회귀)
# + 유의 유전자 발현 패턴 클러스터링 (see.genes)
# + (선택) CMG-SeqViewer용 parquet + staging JSON export
#
# 사용법: Rscript 01c_run_masigpro_timeseries.R [config_path] [output_dir]

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(DESeq2)
  library(edgeR)
  library(MASS)
  library(maSigPro)
  library(openxlsx)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 01c_run_masigpro_timeseries.R [config_path] [output_dir]")
}
config_path <- args[1]
output_dir  <- args[2]

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. Config 로드 & time_series 설정 확인 ---
config <- yaml.load_file(config_path)
ts_cfg <- config$de_analysis$time_series

if (is.null(ts_cfg) || !isTRUE(ts_cfg$enabled)) {
  cat("[01c_run_masigpro_timeseries] time_series.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

time_var     <- ts_cfg$time_variable %||% "time"
series_var   <- ts_cfg$series_variable                 # NULL이면 단일 시리즈
degree_cfg   <- ts_cfg$degree %||% "auto"
q_value      <- ts_cfg$q_value %||% 0.05
rsq_cutoff   <- ts_cfg$rsq_cutoff %||% 0.6
pattern_k    <- ts_cfg$pattern_k %||% 6
export_sv    <- isTRUE(ts_cfg$export_seqviewer)

cat("[01c_run_masigpro_timeseries] Starting maSigPro time-series analysis\n")
cat(paste("  time_variable  :", time_var, "\n"))
cat(paste("  series_variable:", series_var %||% "(single series)", "\n"))
cat(paste("  degree         :", degree_cfg, "\n"))
cat(paste("  q_value        :", q_value, "\n"))
cat(paste("  rsq_cutoff     :", rsq_cutoff, "\n"))

# --- 3. 데이터 로드 & 샘플 정렬 (load_data.R create_de_object()와 동일한 정합 로직) ---
counts <- read.csv(here(config$count_data_path), row.names = 1, check.names = FALSE)
meta   <- read.csv(here(config$metadata_path),   row.names = 1)

meta_samples      <- rownames(meta)
missing_in_counts <- meta_samples[!meta_samples %in% colnames(counts)]
if (length(missing_in_counts) > 0) {
  stop(paste("Metadata samples not found in count data:", paste(missing_in_counts, collapse = ", ")))
}
counts <- counts[, meta_samples, drop = FALSE]

# --- 3b. (선택) include_groups로 시계열에 쓸 샘플만 제한 ---
# 한 프로젝트 안에서 시간축이 다른 여러 시계열(예: acute 0/1/3일 vs chronic 0/24/48시간)을
# 공존시킬 때, 같은 metadata의 group_variable 값으로 서브셋을 나눠 각각 별도 실행하기 위함
# (10_run_coexpression_modules.R의 동일 기능과 같은 패턴). null/미설정이면 전체 샘플 사용.
include_groups <- ts_cfg$include_groups
if (!is.null(include_groups)) {
  group_var <- config$de_analysis$group_variable
  if (is.null(group_var) || !group_var %in% colnames(meta)) {
    stop(paste0("[01c_run_masigpro_timeseries] include_groups를 쓰려면 de_analysis.group_variable('",
                group_var, "')이 metadata에 있어야 합니다."))
  }
  include_groups <- as.character(include_groups)
  all_groups <- unique(as.character(meta[[group_var]]))
  unknown_groups <- setdiff(include_groups, all_groups)
  if (length(unknown_groups) > 0) {
    stop(paste0("[01c_run_masigpro_timeseries] include_groups에 '", group_var, "' 컬럼에 없는 값이 있습니다: ",
                paste(unknown_groups, collapse = ", ")))
  }
  keep_samples <- rownames(meta)[as.character(meta[[group_var]]) %in% include_groups]
  if (length(keep_samples) < 3) {
    stop("[01c_run_masigpro_timeseries] include_groups 필터 적용 후 샘플이 3개 미만입니다.")
  }
  meta   <- meta[keep_samples, , drop = FALSE]
  counts <- counts[, keep_samples, drop = FALSE]
  cat(paste("[01c_run_masigpro_timeseries] include_groups filter:", paste(include_groups, collapse = ", "),
            "->", length(keep_samples), "samples\n"))
}

# --- 4. time / series 컬럼 검증 ---
if (!time_var %in% colnames(meta)) {
  stop(paste0("[01c_run_masigpro_timeseries] time_variable '", time_var,
              "' not found in metadata. Add a numeric time column to enable time_series analysis."))
}
meta[[time_var]] <- suppressWarnings(as.numeric(meta[[time_var]]))
if (any(is.na(meta[[time_var]]))) {
  stop(paste0("[01c_run_masigpro_timeseries] time_variable '", time_var,
              "' contains non-numeric values that could not be coerced."))
}

if (!is.null(series_var) && !series_var %in% colnames(meta)) {
  stop(paste0("[01c_run_masigpro_timeseries] series_variable '", series_var,
              "' not found in metadata."))
}
series_groups <- if (!is.null(series_var)) as.character(meta[[series_var]]) else rep("all", nrow(meta))

# 시리즈별 최소 3개의 서로 다른 time point 요구 (maSigPro 하드 제약)
n_time_per_series <- sapply(split(meta[[time_var]], series_groups), function(x) length(unique(x)))
bad_series <- names(n_time_per_series)[n_time_per_series < 3]
if (length(bad_series) > 0) {
  stop(paste0("[01c_run_masigpro_timeseries] The following series have fewer than 3 distinct time points ",
              "(maSigPro requires >=3): ", paste(bad_series, collapse = ", "),
              ". n_timepoints = ", paste(n_time_per_series[bad_series], collapse = ", ")))
}

min_n_time <- min(n_time_per_series)
if (identical(degree_cfg, "auto")) {
  degree <- min(min_n_time - 1, 2)
} else {
  degree <- as.integer(degree_cfg)
  if (degree > min_n_time - 1) {
    stop(paste0("[01c_run_masigpro_timeseries] degree=", degree, " exceeds max allowed (n_timepoints-1=",
                min_n_time - 1, ") for the series with the fewest time points."))
  }
}
cat(paste("[01c_run_masigpro_timeseries] Using degree =", degree,
          "(min distinct time points across series =", min_n_time, ")\n"))

# --- 5. edesign 행렬 구성 ---
# Replicates: 동일 Time x Series 조합에 속한 샘플은 동일한 값을 가짐 (개별 replicate ID 아님)
interaction_key <- paste(meta[[time_var]], series_groups, sep = "___")
replicates_id   <- as.integer(factor(interaction_key, levels = unique(interaction_key)))

if (is.null(series_var)) {
  edesign <- data.frame(Time = meta[[time_var]], Replicates = replicates_id, Group = 1L)
} else {
  series_levels <- sort(unique(series_groups))
  dummy_cols <- as.data.frame(sapply(series_levels, function(lv) as.integer(series_groups == lv)))
  colnames(dummy_cols) <- make.names(series_levels)
  edesign <- cbind(data.frame(Time = meta[[time_var]], Replicates = replicates_id), dummy_cols)
}
rownames(edesign) <- rownames(meta)
edesign <- edesign[colnames(counts), , drop = FALSE]  # counts 컬럼 순서와 정렬

design_mat_masigpro <- make.design.matrix(edesign, degree = degree)

# --- 6. Pre-filtering (edgeR::filterByExpr) ---
# design_mat_masigpro$dis는 시간 다항식(연속형) 컬럼을 포함하므로 filterByExpr에 그대로 넘기면
# 그룹 구조를 제대로 인식하지 못해 필터가 지나치게 느슨해진다. 대신 실제 Time x Series 조합
# (edesign$Replicates)을 이산 그룹으로 넘겨 각 조합 셀 기준으로 엄격하게 필터링한다.
# 발현이 극히 낮은/희소한 유전자를 남겨두면 count 모드 NB GLM이 수렴하지 않거나
# glm.fit에서 NA/Inf가 발생할 수 있음 (실제 확인된 실패 사례).
dge_all <- DGEList(counts = as.matrix(counts))
keep_genes <- filterByExpr(dge_all, group = factor(edesign$Replicates))
counts_filtered <- as.matrix(counts[keep_genes, , drop = FALSE])
cat(paste("[01c_run_masigpro_timeseries] Pre-filtering (filterByExpr, group=Time x Series): kept",
          sum(keep_genes), "of", length(keep_genes), "genes\n"))

# --- 7. Dispersion(theta) 외부 추정 (edgeR) ---
# maSigPro의 count 모드는 유전자별 dispersion을 자동 추정하지 않으므로
# edgeR의 common dispersion을 theta로 변환해 전달한다.
# 주의: 분산 추정에는 시간 다항식(연속형) design이 아니라 실제 실험 셀(Time x Series)을
# 나타내는 factor 모델을 사용해야 함. 다항식 design으로 추정하면 잔차 자유도가 적어
# dispersion이 비정상적으로 크게(theta가 매우 작게) 추정되어 이후 p.vector의 NB GLM이
# 수렴하지 않거나 NA/Inf가 발생할 수 있음 (실제 확인된 실패 사례).
dge <- DGEList(counts = counts_filtered)
dge <- calcNormFactors(dge)
cell_design_mat <- model.matrix(~ factor(edesign$Replicates))
dge <- estimateGLMCommonDisp(dge, cell_design_mat)
theta <- 1 / dge$common.dispersion
cat(paste("[01c_run_masigpro_timeseries] Estimated theta (1/common.dispersion) =", round(theta, 3), "\n"))

# --- 8. maSigPro 회귀 (p.vector -> T.fit -> get.siggenes) ---
min_obs <- min(6, nrow(edesign))
fit <- p.vector(counts_filtered, design_mat_masigpro, counts = TRUE,
                 family = negative.binomial(theta), theta = theta,
                 Q = q_value, MT.adjust = "BH", min.obs = min_obs)
cat(paste("[01c_run_masigpro_timeseries] p.vector:", fit$g, "genes tested,",
          nrow(fit$SELEC), "with significant model fit (Q =", q_value, ")\n"))

tstep <- T.fit(fit, step.method = "backward", alfa = q_value)
sigs  <- get.siggenes(tstep, rsq = rsq_cutoff, vars = "all")

sig_gene_ids <- unique(as.character(rownames(sigs$sig.genes$sig.profiles)))
cat(paste("[01c_run_masigpro_timeseries]", length(sig_gene_ids), "significant genes (R-squared >=", rsq_cutoff, ")\n"))

# --- 9. VST 값 (클러스터 패턴 시각화 + seqviewer export용) ---
dds_vst <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = meta, design = ~1)
dds_vst <- estimateSizeFactors(dds_vst)
vsd     <- vst(dds_vst, blind = TRUE)
vst_mat <- assay(vsd)

sample_order <- rownames(edesign)[order(series_groups[match(rownames(edesign), rownames(meta))],
                                         edesign$Time)]

if (length(sig_gene_ids) == 0) {
  cat("[01c_run_masigpro_timeseries] No significant genes found — writing empty outputs.\n")
  empty_df <- data.frame(gene_id = character(0), p_value = numeric(0),
                          r_squared = numeric(0), cluster_id = integer(0))
  write.csv(empty_df, file.path(output_dir, "time_series_significant_genes.csv"), row.names = FALSE)
  file.copy(config_path, file.path(output_dir, "config_used.yml"), overwrite = TRUE)
  cat("[01c_run_masigpro_timeseries] Done (no significant genes).\n")
  quit(save = "no", status = 0)
}

# --- 10. 패턴 클러스터링 (see.genes) ---
# maSigPro의 see.genes/PlotGroups는 유의 유전자가 1개일 때 행렬이 벡터로 축소되며
# "incorrect number of dimensions" 에러로 죽음 (실제 확인된 실패 사례) — 이 경우
# 클러스터링을 건너뛰고 단일 유전자를 cluster_id=1로 취급한다.
if (length(sig_gene_ids) < 2) {
  cat("[01c_run_masigpro_timeseries] Only 1 significant gene — skipping see.genes clustering.\n")
  cluster_assignment <- setNames(rep(1L, length(sig_gene_ids)), sig_gene_ids)
} else {
  pattern_plot_path <- file.path(output_dir, "time_series_pattern_plot.png")
  png(pattern_plot_path, width = 12, height = 10, units = "in", res = 300, bg = "white")
  sg_result <- see.genes(vst_mat[sig_gene_ids, sample_order, drop = FALSE],
                          edesign = edesign[sample_order, , drop = FALSE],
                          k = min(pattern_k, length(sig_gene_ids)),
                          cluster.method = "hclust", newX11 = FALSE)
  dev.off()
  cat(paste("[01c_run_masigpro_timeseries] Pattern plot saved:", pattern_plot_path, "\n"))
  cluster_assignment <- sg_result$cut[sig_gene_ids]
}

# --- 11. 통계값(p-value, R-squared) 추출 ---
pval_df <- sigs$sig.genes$sig.pvalues
pval_col <- grep("^p.value$|p-value", colnames(pval_df), value = TRUE, ignore.case = TRUE)[1]
rsq_col  <- grep("R-squared|r.squared", colnames(pval_df), value = TRUE, ignore.case = TRUE)[1]

final_df <- data.frame(
  gene_id     = sig_gene_ids,
  p_value     = pval_df[sig_gene_ids, pval_col],
  r_squared   = pval_df[sig_gene_ids, rsq_col],
  cluster_id  = as.integer(cluster_assignment[sig_gene_ids]),
  stringsAsFactors = FALSE
)

# gene_symbol 주석 (선택, 실패해도 파이프라인은 계속 진행)
species     <- config$species %||% "human"
org_db_name <- config$databases[[species]]$organism_db
if (!is.null(org_db_name)) {
  tryCatch({
    suppressPackageStartupMessages(requireNamespace(org_db_name, quietly = TRUE))
    org_db <- get(org_db_name, envir = loadNamespace(org_db_name))
    gene_id_type <- config$gene_id_type %||% "ENSEMBL"
    symbols <- AnnotationDbi::mapIds(org_db, keys = sig_gene_ids, column = "SYMBOL",
                                      keytype = gene_id_type, multiVals = "first")
    symbols[is.na(symbols)] <- names(symbols)[is.na(symbols)]
    final_df <- cbind(gene_symbol = symbols[final_df$gene_id], final_df)
  }, error = function(e) {
    cat(paste("[01c_run_masigpro_timeseries] gene_symbol annotation skipped:", conditionMessage(e), "\n"))
  })
}

vst_ordered <- as.data.frame(vst_mat[sig_gene_ids, sample_order, drop = FALSE])
final_df <- cbind(final_df, vst_ordered[final_df$gene_id, , drop = FALSE])

# --- 12. 결과 저장 (CSV + xlsx + config_used.yml) ---
output_csv_path <- file.path(output_dir, "time_series_significant_genes.csv")
write.csv(final_df, output_csv_path, row.names = FALSE)
cat(paste("[01c_run_masigpro_timeseries] CSV saved:", output_csv_path, "\n"))

if (isTRUE(config$export$export_to_excel)) {
  wb <- createWorkbook()
  addWorksheet(wb, "Significant_Genes")
  writeData(wb, "Significant_Genes", final_df, rowNames = FALSE)
  saveWorkbook(wb, file.path(output_dir, "time_series_significant_genes.xlsx"), overwrite = TRUE)
  cat("[01c_run_masigpro_timeseries] Excel file saved.\n")
}

file.copy(config_path, file.path(output_dir, "config_used.yml"), overwrite = TRUE)

# ─────────────────────────────────────────────────────────────
# 13. CMG-SeqViewer export (parquet + staging JSON) — export_seqviewer: true 일 때만
#     09_export_multi_group.R과 동일한 헬퍼/디렉토리 규격 사용
# ─────────────────────────────────────────────────────────────
if (export_sv) {
  suppressPackageStartupMessages({
    library(arrow)
    library(jsonlite)
    library(tibble)
  })

  make_alias_slug <- function(alias, max_len = 80) {
    slug <- gsub("[^\\w가-힣]+", "_", alias, perl = TRUE)
    slug <- gsub("^_|_$", "", slug)
    substr(slug, 1, max_len)
  }

  new_uuid <- function() {
    hex <- paste0(sample(c(0:9, letters[1:6]), 32, replace = TRUE), collapse = "")
    paste(substr(hex,1,8), substr(hex,9,12),
          paste0("4", substr(hex,14,16)),
          paste0(sample(c("8","9","a","b"),1), substr(hex,18,20)),
          substr(hex,21,32), sep = "-")
  }

  write_parquet_dataset <- function(df, alias, datasets_dir) {
    uid      <- new_uuid()
    slug     <- make_alias_slug(alias)
    filename <- paste0(slug, ".parquet")
    old_files <- list.files(datasets_dir, pattern = paste0("^", slug, "\\.parquet$"), full.names = TRUE)
    if (length(old_files) > 0) file.remove(old_files)
    write_parquet(df, file.path(datasets_dir, filename))
    list(uid = uid, filename = filename)
  }

  seqviewer_dir <- file.path(config$output_dir, "seqviewer")
  datasets_dir  <- file.path(seqviewer_dir, "datasets")
  staging_dir   <- file.path(seqviewer_dir, "staging")
  dir.create(datasets_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(staging_dir,  recursive = TRUE, showWarnings = FALSE)

  parquet_df <- tibble::rownames_to_column(final_df, var = "row_id")
  parquet_df$row_id <- NULL  # gene_id는 이미 컬럼으로 존재

  series_label <- if (!is.null(series_var)) paste(sort(unique(series_groups)), collapse = " / ") else "single series"

  # 같은 프로젝트 안에서 include_groups로 시간축이 다른 여러 시계열(acute/chronic 등)을
  # 따로 export할 수 있도록 variant_label이 있으면 alias/parquet slug/staging 파일명에
  # 반영해 서로 겹치지 않게 한다(10_run_coexpression_modules.R과 동일 패턴).
  variant_label <- ts_cfg$variant_label
  ts_alias <- if (!is.null(variant_label)) {
    paste(basename(config$output_dir), "Time-Series", variant_label)
  } else {
    paste(basename(config$output_dir), "Time-Series")
  }
  ts_info  <- write_parquet_dataset(parquet_df, ts_alias, datasets_dir)

  ts_entry <- list(
    dataset_id           = ts_info$uid,
    alias                = ts_alias,
    original_filename     = ts_info$filename,
    dataset_type         = "time_series",
    experiment_condition = series_label,
    organism             = config$species %||% "",
    cell_type            = "",
    tissue               = "",
    timepoint             = paste(sort(unique(meta[[time_var]])), collapse = ", "),
    row_count            = nrow(parquet_df),
    gene_count           = nrow(parquet_df),
    significant_genes    = nrow(parquet_df),
    import_date          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    file_path            = ts_info$filename,
    notes                = paste0("maSigPro count-based regression (degree=", degree,
                                    ", Q=", q_value, ", R2>=", rsq_cutoff, ") + see.genes clustering (k=", pattern_k, ")",
                                    if (!is.null(variant_label)) paste0(" [subset: ", variant_label, "]") else ""),
    tags                 = as.list(c("time_series", "maSigPro",
                                      if (!is.null(variant_label)) variant_label else NULL,
                                      sort(unique(series_groups))))
  )

  staging_filename <- if (!is.null(variant_label)) {
    paste0("time_series_", make_alias_slug(variant_label), "_entries.json")
  } else {
    "time_series_entries.json"
  }
  staging_path <- file.path(staging_dir, staging_filename)
  write_json(list(ts_entry), staging_path, pretty = TRUE, auto_unbox = TRUE)
  cat(paste("[01c_run_masigpro_timeseries] Seqviewer parquet saved:", ts_info$filename, "\n"))
  cat(paste("[01c_run_masigpro_timeseries] Seqviewer staging JSON saved:", staging_path, "\n"))
}

cat("[01c_run_masigpro_timeseries] Done.\n")
