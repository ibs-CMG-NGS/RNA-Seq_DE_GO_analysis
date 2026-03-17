# 파일 경로: src/analysis/06_export_seqviewer.R
# cmg-seqviewer용 parquet + staging metadata JSON 생성 (per pair)
# Usage: Rscript 06_export_seqviewer.R [config_path] [compare_group] [base_group] [pair_output_dir]

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(jsonlite)
  library(dplyr)
  library(openxlsx)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 06_export_seqviewer.R [config_path] [compare_group] [base_group] [pair_output_dir]")
}
config_path    <- args[1]
compare_group  <- args[2]
base_group     <- args[3]
pair_output_dir <- args[4]   # e.g. output/mouse-monSTIM-2026/pairwise/transgenic_vs_wildtype

config <- yaml.load_file(config_path)

pair_label   <- paste0(compare_group, "_vs_", base_group)
output_dir   <- config$output_dir                                # e.g. output/mouse-monSTIM-2026
seqviewer_dir <- file.path(output_dir, "seqviewer")
datasets_dir  <- file.path(seqviewer_dir, "datasets")
staging_dir   <- file.path(seqviewer_dir, "staging")

dir.create(datasets_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir,  recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────
# 헬퍼 함수 (pipeline-integration.md 규격)
# ─────────────────────────────────────────────────────────────

make_alias_slug <- function(alias, max_len = 40) {
  slug <- gsub("[^\\w\uac00-\ud7a3]+", "_", alias, perl = TRUE)
  slug <- gsub("^_|_$", "", slug)
  substr(slug, 1, max_len)
}

new_uuid <- function() {
  hex <- paste0(sample(c(0:9, letters[1:6]), 32, replace = TRUE), collapse = "")
  paste(
    substr(hex,  1,  8),
    substr(hex,  9, 12),
    paste0("4", substr(hex, 14, 16)),
    paste0(sample(c("8", "9", "a", "b"), 1), substr(hex, 18, 20)),
    substr(hex, 21, 32),
    sep = "-"
  )
}

write_parquet_dataset <- function(df, alias, datasets_dir) {
  uid      <- new_uuid()
  slug     <- make_alias_slug(alias)
  filename <- paste0(slug, ".parquet")
  # 같은 slug의 이전 파일 제거 (재실행 시 중복 방지)
  old_files <- list.files(datasets_dir, pattern = paste0("^", slug, "\\.parquet$"), full.names = TRUE)
  if (length(old_files) > 0) file.remove(old_files)
  write_parquet(df, file.path(datasets_dir, filename))
  list(uid = uid, filename = filename)
}

# ─────────────────────────────────────────────────────────────
# 2. DE 결과 → parquet  (final_de_results.xlsx DE_Results 시트)
# ─────────────────────────────────────────────────────────────
de_xlsx <- file.path(pair_output_dir, "final_de_results.xlsx")
if (!file.exists(de_xlsx)) {
  stop(paste("DE results not found:", de_xlsx))
}

# DE_Results 시트: gene_id, symbol, DESeq2 통계량, 샘플별 normalized counts 전체 포함
de_raw <- read.xlsx(de_xlsx, sheet = "DE_Results", check.names = FALSE)

de_std <- de_raw %>%
  rename(
    any_of(c(
      base_mean  = "baseMean",
      log2fc     = "log2FoldChange",
      lfcse      = "lfcSE",
      adj_pvalue = "padj"
    ))
  )

if (!"symbol" %in% colnames(de_std)) {
  de_std$symbol <- de_std$gene_id
}

# baseMean = 0인 유전자 제거 (전체 샘플 미발현 → 통계값 전부 NA)
if ("base_mean" %in% colnames(de_std)) {
  n_before <- nrow(de_std)
  de_std   <- de_std %>% filter(!is.na(base_mean) & base_mean > 0)
  cat(paste("Filtered", n_before - nrow(de_std), "zero-expression genes;",
            nrow(de_std), "genes retained\n"))
}

organism  <- config$species %||% ""
condition <- paste(compare_group, "vs", base_group)
de_alias  <- paste(condition, "DE")

de_info   <- write_parquet_dataset(de_std, de_alias, datasets_dir)

padj_cut <- config$de_analysis$padj_cutoff  %||% 0.05
lfc_cut  <- config$de_analysis$log2fc_cutoff %||% 1.0

sig_count <- sum(!is.na(de_std$adj_pvalue) & !is.na(de_std$log2fc) &
                   de_std$adj_pvalue < padj_cut &
                   abs(de_std$log2fc) >= lfc_cut,
                 na.rm = TRUE)

de_entry <- list(
  dataset_id           = de_info$uid,
  alias                = de_alias,
  original_filename    = de_info$filename,
  dataset_type         = "differential_expression",
  experiment_condition = condition,
  organism             = organism,
  cell_type            = "",
  tissue               = "",
  timepoint            = "",
  row_count            = nrow(de_std),
  gene_count           = nrow(de_std),
  significant_genes    = sig_count,
  import_date          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  file_path            = de_info$filename,
  notes                = paste0(config$de_analysis$method %||% "DESeq2"),
  tags                 = as.list(c("DE", compare_group, base_group))
)
cat(paste("DE parquet saved:", de_info$filename, "\n"))

# ─────────────────────────────────────────────────────────────
# 3. GO/KEGG 결과 → parquet (final_go_results.xlsx 전체를 단일 parquet)
# ─────────────────────────────────────────────────────────────
go_entries <- list()

go_xlsx <- file.path(pair_output_dir, "final_go_results.xlsx")
if (file.exists(go_xlsx)) {
  sheets <- getSheetNames(go_xlsx)
  data_sheets <- sheets[sheets != "Analysis_Info"]

  all_go <- lapply(data_sheets, function(sh) {
    df <- read.xlsx(go_xlsx, sheet = sh, check.names = FALSE)
    if (nrow(df) == 0) return(NULL)

    # GO 시트: GO.ID / GO.Term  |  KEGG 시트: KEGG.ID / KEGG.Pathway
    df <- df %>%
      rename(any_of(c(
        term_id      = "GO.ID",
        term_id      = "KEGG.ID",
        description  = "GO.Term",
        description  = "KEGG.Pathway",
        direction    = "Gene.Set",
        gene_ratio   = "Gene.Ratio",
        bg_ratio     = "Background.Ratio",
        pvalue       = "P-value",
        fdr          = "Adjusted.P-value",
        qvalue       = "Q-value",
        gene_count   = "Gene.Count",
        gene_symbols = "Gene.Symbols"
      )))

    # KEGG 시트는 Ontology 컬럼이 없으므로 추가
    if (!"Ontology" %in% colnames(df)) df$Ontology <- "KEGG"
    df <- rename(df, any_of(c(ontology = "Ontology")))
    df
  })

  all_go <- Filter(Negate(is.null), all_go)

  if (length(all_go) > 0) {
    go_combined <- bind_rows(all_go)
    go_alias <- paste(condition, "GO+KEGG")
    go_info  <- write_parquet_dataset(go_combined, go_alias, datasets_dir)

    go_entries[[1]] <- list(
      dataset_id           = go_info$uid,
      alias                = go_alias,
      original_filename    = go_info$filename,
      dataset_type         = "go_analysis",
      experiment_condition = condition,
      organism             = organism,
      cell_type            = "",
      tissue               = "",
      timepoint            = "",
      row_count            = nrow(go_combined),
      gene_count           = 0L,
      significant_genes    = 0L,
      import_date          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      file_path            = go_info$filename,
      notes                = "clusterProfiler enrichGO + enrichKEGG (all gene sets)",
      tags                 = as.list(c("GO", "KEGG", compare_group, base_group))
    )
    cat(paste("GO+KEGG parquet saved:", go_info$filename,
              "(", nrow(go_combined), "terms )\n"))
  }
}

# ─────────────────────────────────────────────────────────────
# 4. staging JSON 저장 (병렬 실행 안전)
# ─────────────────────────────────────────────────────────────
all_entries <- c(list(de_entry), go_entries)
staging_path <- file.path(staging_dir, paste0(pair_label, "_entries.json"))
write_json(all_entries, staging_path, pretty = TRUE, auto_unbox = TRUE)
cat(paste("Staging JSON saved:", staging_path,
          "(", length(all_entries), "entries )\n"))

# 완료 플래그
flag_path <- file.path(pair_output_dir, ".seqviewer_export_done.flag")
file.create(flag_path)
cat("seqviewer export done.\n")
