# 파일 경로: src/analysis/09_export_multi_group.R
# Multi-group result export: omnibus LRT 통계 + size-factor normalized counts 통합
# CMG-SeqViewer MULTI_GROUP 타입용 CSV + Parquet 생성
#
# 사용법:
#   Rscript 09_export_multi_group.R <config_path> <omnibus_csv_path> <output_csv_path>
#   (parquet는 output_csv_path의 .csv → .parquet로 자동 저장)

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(arrow)
  library(jsonlite)
})

# NULL-coalescing 헬퍼
`%||%` <- function(a, b) if (!is.null(a)) a else b

# seqviewer 헬퍼 (06_export_seqviewer.R 동일 규격)
make_alias_slug <- function(alias, max_len = 40) {
  slug <- gsub("[^\\w\uac00-\ud7a3]+", "_", alias, perl = TRUE)
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
reference_group     <- mg_cfg$reference_group # NULL이면 자동 감지

cat("[09_export_multi_group] Starting multi-group result export\n")
cat(paste("  abundance_type     :", abundance_type, "\n"))
cat(paste("  include_gene_symbol:", include_gene_symbol, "\n"))
cat(paste("  filter_padj        :", filter_padj %||% "none", "\n"))
cat(paste("  filter_basemean    :", filter_basemean %||% "none", "\n"))
cat(paste("  reference_group    :", reference_group %||% "auto-detect", "\n"))

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
# reference_group을 맨 앞에, 나머지는 오름차순
all_groups <- sort(unique(as.character(meta[[group_var]])))

# reference_group 미지정 시 pairwise_comparisons base 그룹 최빈값으로 자동 감지
if (is.null(reference_group)) {
  pairs <- de_cfg$pairwise_comparisons
  if (!is.null(pairs) && length(pairs) > 0) {
    bases        <- sapply(pairs, function(p) p[[2]])
    reference_group <- names(sort(table(bases), decreasing = TRUE))[1]
    cat(paste("[09_export_multi_group] reference_group auto-detected:", reference_group, "\n"))
  }
}

if (!is.null(reference_group) && reference_group %in% all_groups) {
  ordered_groups <- c(reference_group, sort(setdiff(all_groups, reference_group)))
} else {
  if (!is.null(reference_group))
    warning("[09_export_multi_group] reference_group '", reference_group,
            "' not found in metadata — falling back to alphabetical order.")
  ordered_groups <- all_groups
}

sample_order <- unlist(lapply(ordered_groups, function(grp)
  sort(rownames(meta)[meta[[group_var]] == grp])))
norm_ordered <- norm_mat[, sample_order, drop = FALSE]

cat(paste("[09_export_multi_group] Group order:", paste(ordered_groups, collapse = " → "), "\n"))
cat(paste("[09_export_multi_group] Sample order:", paste(sample_order, collapse = ", "), "\n"))

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

# --- 10. Save CSV (root) + Parquet + Staging JSON (seqviewer) ---
write.csv(final_df, output_csv_path, row.names = TRUE)
cat(paste("[09_export_multi_group] CSV saved:", output_csv_path, "\n"))

# seqviewer 디렉토리 구조 생성
seqviewer_dir <- file.path(config$output_dir, "seqviewer")
datasets_dir  <- file.path(seqviewer_dir, "datasets")
staging_dir   <- file.path(seqviewer_dir, "staging")
dir.create(datasets_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir,  recursive = TRUE, showWarnings = FALSE)

# parquet → seqviewer/datasets/
parquet_df <- tibble::rownames_to_column(final_df, var = "gene_id")
mg_alias   <- paste(config$output_dir, "multi_group") # 프로젝트별 고유 alias
mg_alias   <- basename(config$output_dir)             # e.g. "kkj-rna-seq-ctx"
mg_alias   <- paste(mg_alias, "Multi-Group")
mg_info    <- write_parquet_dataset(parquet_df, mg_alias, datasets_dir)
cat(paste("[09_export_multi_group] Parquet saved:", file.path(datasets_dir, mg_info$filename), "\n"))

# staging JSON entry (06b_aggregate_seqviewer.R가 수집)
mg_entry <- list(
  dataset_id           = mg_info$uid,
  alias                = mg_alias,
  original_filename    = mg_info$filename,
  dataset_type         = "multi_group",
  experiment_condition = paste(ordered_groups, collapse = " / "),
  organism             = config$species %||% "",
  cell_type            = "",
  tissue               = "",
  timepoint            = "",
  row_count            = nrow(parquet_df),
  gene_count           = nrow(parquet_df),
  significant_genes    = if ("padj" %in% colnames(final_df))
                           sum(!is.na(final_df$padj) & final_df$padj <= 0.05, na.rm = TRUE)
                         else 0L,
  import_date          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  file_path            = mg_info$filename,
  notes                = paste0(de_cfg$method %||% "DESeq2", " LRT omnibus + normalized counts"),
  tags                 = as.list(c("multi_group", ordered_groups))
)

staging_path <- file.path(staging_dir, "multi_group_entries.json")
write_json(list(mg_entry), staging_path, pretty = TRUE, auto_unbox = TRUE)
cat(paste("[09_export_multi_group] Staging JSON saved:", staging_path, "\n"))
cat(paste("[09_export_multi_group] Done:", nrow(final_df), "genes exported\n"))
