# 파일 경로: src/analysis/09_export_multi_group.R
# Multi-group result export: omnibus LRT 통계 + size-factor normalized counts 통합
# CMG-SeqViewer MULTI_GROUP 타입용 CSV 생성
#
# 사용법:
#   Rscript 09_export_multi_group.R <config_path> <omnibus_csv_path> <output_csv_path>

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# NULL-coalescing 헬퍼
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- 1. Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript 09_export_multi_group.R <config_path> <omnibus_csv_path> <output_csv_path>")
}
config_path      <- args[1]
omnibus_csv_path <- args[2]
output_csv_path  <- args[3]

# --- 2. Load config ---
config <- yaml.load_file(config_path)
de_cfg <- config$de_analysis
mg_cfg <- de_cfg$multi_group_export

if (is.null(mg_cfg) || !isTRUE(mg_cfg$enabled)) {
  cat("[09_export_multi_group] multi_group_export.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

# abundance_type 분기 — VST는 미구현
abundance_type <- mg_cfg$abundance_type %||% "normalized"
if (abundance_type == "vst") {
  stop(
    "[09_export_multi_group] abundance_type='vst' is not yet implemented.\n",
    "Set abundance_type: 'normalized' in config.\n",
    "# TODO: vst 구현 시 아래 추가:\n",
    "#   vsd <- vst(dds, blind = config$de_analysis$advanced_options$vst_blind %||% FALSE)\n",
    "#   norm_mat <- assay(vsd)\n"
  )
}

include_gene_symbol <- isTRUE(mg_cfg$include_gene_symbol)
filter_padj         <- mg_cfg$filter_padj     # NULL이면 필터 없음
filter_basemean     <- mg_cfg$filter_basemean # NULL이면 필터 없음

cat("[09_export_multi_group] Starting multi-group result export\n")
cat(paste("  abundance_type     :", abundance_type, "\n"))
cat(paste("  include_gene_symbol:", include_gene_symbol, "\n"))
cat(paste("  filter_padj        :", filter_padj %||% "none", "\n"))
cat(paste("  filter_basemean    :", filter_basemean %||% "none", "\n"))

# --- 3. Load omnibus stats (LRT 재실행 없음) ---
omnibus_df <- read.csv(omnibus_csv_path, row.names = 1, check.names = FALSE)
stat_cols  <- intersect(c("baseMean", "stat", "pvalue", "padj"), colnames(omnibus_df))
stats_df   <- omnibus_df[, stat_cols, drop = FALSE]
cat(paste("[09_export_multi_group] Omnibus stats loaded:", nrow(stats_df), "genes,",
          "columns:", paste(stat_cols, collapse = ", "), "\n"))

# --- 4. Load count data & metadata ---
counts <- read.csv(here(config$count_data_path), row.names = 1, check.names = FALSE)
meta   <- read.csv(here(config$metadata_path),   row.names = 1)

# metadata에 있는 샘플만 사용 (공유 counts 파일 대응)
meta_samples      <- rownames(meta)
missing_in_counts <- meta_samples[!meta_samples %in% colnames(counts)]
if (length(missing_in_counts) > 0) {
  stop(paste("Metadata samples not found in count data:",
             paste(missing_in_counts, collapse = ", ")))
}
counts <- counts[, meta_samples, drop = FALSE]

group_var      <- de_cfg$group_variable
meta[[group_var]] <- as.factor(meta[[group_var]])
design_formula    <- as.formula(de_cfg$design_formula)

cat(paste("[09_export_multi_group] Loaded", nrow(counts), "genes,", ncol(counts), "samples\n"))

# --- 5. Normalized counts (estimateSizeFactors만 실행 — full DESeq() 불필요) ---
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
dds <- estimateSizeFactors(dds)
norm_mat <- counts(dds, normalized = TRUE)

# --- 6. 샘플 컬럼 그룹별 정렬 ---
group_levels  <- levels(meta[[group_var]])
sample_order  <- unlist(lapply(group_levels, function(grp)
  sort(rownames(meta)[meta[[group_var]] == grp])))
norm_ordered  <- norm_mat[, sample_order, drop = FALSE]

cat(paste("[09_export_multi_group] Sample order:",
          paste(sample_order, collapse = ", "), "\n"))

# --- 7. 통합 (omnibus stats + normalized counts) ---
common_genes <- intersect(rownames(stats_df), rownames(norm_ordered))
if (length(common_genes) == 0) {
  stop("[09_export_multi_group] No common genes between omnibus CSV and count matrix.")
}
final_df <- cbind(
  stats_df[common_genes, , drop = FALSE],
  as.data.frame(norm_ordered[common_genes, , drop = FALSE])
)
cat(paste("[09_export_multi_group] Merged:", nrow(final_df), "genes\n"))

# --- 8. 선택적 필터링 ---
if (!is.null(filter_padj) && !is.na(filter_padj) && "padj" %in% colnames(final_df)) {
  before <- nrow(final_df)
  final_df <- final_df[!is.na(final_df$padj) & final_df$padj <= filter_padj, ]
  cat(paste("[09_export_multi_group] padj filter (<=", filter_padj, "):",
            before, "->", nrow(final_df), "genes\n"))
}
if (!is.null(filter_basemean) && !is.na(filter_basemean) && "baseMean" %in% colnames(final_df)) {
  before <- nrow(final_df)
  final_df <- final_df[!is.na(final_df$baseMean) & final_df$baseMean >= filter_basemean, ]
  cat(paste("[09_export_multi_group] baseMean filter (>=", filter_basemean, "):",
            before, "->", nrow(final_df), "genes\n"))
}

# --- 9. gene_symbol annotation (선택) ---
if (include_gene_symbol) {
  species    <- config$species %||% "human"
  org_db_name <- config$databases[[species]]$organism_db
  if (is.null(org_db_name)) {
    warning("[09_export_multi_group] databases.", species,
            ".organism_db not found in config — skipping gene_symbol annotation.")
  } else {
    suppressPackageStartupMessages(
      requireNamespace(org_db_name, quietly = TRUE) ||
        stop("[09_export_multi_group] Package '", org_db_name, "' is not installed.")
    )
    org_db <- get(org_db_name, envir = loadNamespace(org_db_name))
    gene_id_type <- config$gene_id_type %||% "ENSEMBL"

    tryCatch({
      symbols <- AnnotationDbi::mapIds(
        org_db,
        keys      = rownames(final_df),
        column    = "SYMBOL",
        keytype   = gene_id_type,
        multiVals = "first"
      )
      # NA → gene_id fallback
      symbols[is.na(symbols)] <- names(symbols)[is.na(symbols)]
      final_df <- cbind(
        gene_symbol = symbols[rownames(final_df)],
        final_df
      )
      cat(paste("[09_export_multi_group] gene_symbol annotation added\n"))
    }, error = function(e) {
      warning("[09_export_multi_group] gene_symbol annotation failed: ", conditionMessage(e))
    })
  }
}

# --- 10. Save ---
write.csv(final_df, output_csv_path, row.names = TRUE)
cat(paste("[09_export_multi_group] Saved", nrow(final_df), "genes to:", output_csv_path, "\n"))
