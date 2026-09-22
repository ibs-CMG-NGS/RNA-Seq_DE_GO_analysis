# 파일 경로: src/analysis/02e_check_celltype_purity.R
# 목적: cell-type marker gene set 기반 순도(purity) 체크 (모든 샘플 대상)
#   샘플별 VST 발현량 순위에 preranked GSEA(fgsea)를 적용해, 각 cell-type
#   marker set이 발현 상위에 유의하게 몰려 있는지(NES + BH-adjusted p-value)를
#   검정한다. 참조 매트릭스 없이 샘플 1개만으로도 통계적 유의성 판단이 가능.
# 사용법: Rscript 02e_check_celltype_purity.R --config config.yml --output_dir output/qc_plots

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(optparse)
  library(DESeq2)
  library(openxlsx)
  library(fgsea)
  library(pheatmap)
  library(AnnotationDbi)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- 1. Setup: Load libraries and parse arguments ---
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yml",
              help = "Path to the config YAML file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character",
              help = "Directory to save cell-type purity results", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Output directory must be supplied.", call. = FALSE)
}

config <- yaml.load_file(opt$config)
ct_cfg <- config$qc_plots$celltype_check

if (is.null(ct_cfg) || !isTRUE(ct_cfg$enabled)) {
  cat("[02e_check_celltype_purity] qc_plots.celltype_check.enabled이 false입니다. 종료.\n")
  quit(status = 0)
}
if (is.null(ct_cfg$marker_list_path)) {
  stop("[02e_check_celltype_purity] qc_plots.celltype_check.marker_list_path가 설정되지 않았습니다.")
}

out_dir <- file.path(opt$output_dir, "celltype_purity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

marker_id_type   <- toupper(ct_cfg$marker_id_type %||% config$gene_id_type %||% "SYMBOL")
marker_sheets    <- ct_cfg$marker_sheets
expected_type    <- ct_cfg$expected_cell_type
fdr_cutoff       <- ct_cfg$fdr_cutoff %||% 0.05
min_marker_genes <- ct_cfg$min_marker_genes %||% 5

cat("\n=== Cell-type Purity Check (preranked fgsea) ===\n")
cat(paste0("Marker list   : ", ct_cfg$marker_list_path, "\n"))
cat(paste0("Marker ID type: ", marker_id_type, "\n"))
cat(paste0("Expected type : ", expected_type %||% "(none specified)", "\n\n"))

# --- 2. Marker gene list 읽기 (xlsx, wide format: 1행 = cell type, 컬럼마다 marker gene 세로 나열) ---
# NOTE: openxlsx::read.xlsx()의 sep.names 인자는 check.names와 무관하게 공백을
# 마침표로 바꿔버리는 알려진 함정이 있음 — sep.names=" "로 명시해서 방지한다.
#
# 워크북에 cell-type 마커가 아닌 시트(예: 질병 위험유전자 목록)가 섞여 있을 수 있으므로
# marker_sheets로 명시된 시트만 읽는다. 시트 이름에 괄호가 있으면(예: "cell-type-specific
# (Kang)") 괄호 안 텍스트를 출처 라벨로 써서 "{컬럼명} [{출처}]" 형식으로 cell_type을
# 구분한다 — 여러 문헌 출처를 동시에 쓸 때 결과에서 서로 다른 마커셋임을 알 수 있도록.
read_celltype_markers <- function(xlsx_path, marker_sheets) {
  all_sheets <- getSheetNames(xlsx_path)

  if (is.null(marker_sheets)) {
    if (length(all_sheets) != 1) {
      stop("[02e_check_celltype_purity] xlsx에 시트가 여러 개입니다 (",
           paste(all_sheets, collapse = ", "),
           "). qc_plots.celltype_check.marker_sheets에 사용할 시트 이름을 명시해주세요.")
    }
    marker_sheets <- all_sheets
  }
  unknown <- setdiff(marker_sheets, all_sheets)
  if (length(unknown) > 0) {
    stop("[02e_check_celltype_purity] marker_sheets에 워크북에 없는 시트가 있습니다: ",
         paste(unknown, collapse = ", "), " (워크북 시트: ", paste(all_sheets, collapse = ", "), ")")
  }

  marker_list <- list()
  for (sh in marker_sheets) {
    df <- read.xlsx(xlsx_path, sheet = sh, colNames = TRUE, sep.names = " ")
    paren_match <- regmatches(sh, regexpr("(?<=\\()[^)]+(?=\\))", sh, perl = TRUE))
    source_label <- if (length(paren_match) == 1) paren_match else sh

    for (cn in colnames(df)) {
      genes <- trimws(as.character(df[[cn]]))
      genes <- unique(toupper(genes[!is.na(genes) & genes != ""]))
      if (length(genes) == 0) next
      label <- paste0(cn, " [", source_label, "]")
      marker_list[[label]] <- genes
    }
  }
  if (length(marker_list) == 0) {
    stop("[02e_check_celltype_purity] marker gene을 하나도 읽지 못했습니다: ", xlsx_path)
  }
  marker_list
}

marker_list <- read_celltype_markers(here(ct_cfg$marker_list_path), marker_sheets)
cat(paste0("Cell types in marker list (", length(marker_list), " total):\n"))
for (nm in names(marker_list)) {
  cat(paste0("  - ", nm, ": ", length(marker_list[[nm]]), " genes\n"))
}
cat("\n")

# --- 3. DESeq2 dds 생성 + VST ---
source(here("src", "utils", "load_data.R"))
dds <- create_de_object(config_path = opt$config)

if (!"dispersion" %in% names(mcols(dds))) {
  cat("Running DESeq2 analysis...\n")
  dds <- DESeq(dds)
}

adv_opts <- config$de_analysis$advanced_options
blind_mode <- if (!is.null(adv_opts$vst_blind)) adv_opts$vst_blind else FALSE
cat("\nApplying variance stabilizing transformation (VST)...\n")
vsd <- vst(dds, blind = blind_mode)
vst_mat <- assay(vsd)   # genes(rownames = gene_id_type) x samples

# 사실상 발현이 없는(거의 모든 샘플에서 raw count 0) 유전자를 VST 그대로 두면 이
# 유전자들이 랭킹 최하단에 큰 덩어리로 몰려서, 실제로 발현되는 유전자로만 구성된
# 큐레이션 marker gene set은 cell type과 무관하게 항상 유의하게 "상위"로 보이는
# 인위적 신호가 생긴다(실측 확인 — 필터링 전에는 37개 marker set 전체가 모든 샘플에서
# 극단적으로 유의했음). 절반 이상의 샘플에서 count>=1인 유전자만 남겨서 이 artifact를 줄인다.
detectable <- rowSums(counts(dds) >= 1) >= ceiling(ncol(dds) / 2)
cat(paste0("Filtering to detectably-expressed genes: ", sum(detectable), " / ", length(detectable), " genes kept\n"))
vst_mat <- vst_mat[detectable, , drop = FALSE]

# --- 4. gene_id_type -> SYMBOL 매핑 (marker_id_type이 SYMBOL이고 count matrix가 아닌 경우) ---
top_gene_id_type <- toupper(config$gene_id_type %||% "ENSEMBL")

# 같은 rowname(gene)이 여러 개면 평균으로 집계 (fgsea는 유일한 이름 필요)
aggregate_duplicate_rows <- function(mat) {
  if (!any(duplicated(rownames(mat)))) return(mat)
  sums <- rowsum(mat, group = rownames(mat))
  counts <- as.vector(table(rownames(mat))[rownames(sums)])
  sums / counts
}

if (marker_id_type == top_gene_id_type) {
  expr_mat <- vst_mat
  rownames(expr_mat) <- rownames(vst_mat)
} else if (marker_id_type == "SYMBOL") {
  species_info <- config$databases[[config$species]]
  if (is.null(species_info) || is.null(species_info$organism_db)) {
    stop("[02e_check_celltype_purity] databases[", config$species, "].organism_db 설정이 없습니다.")
  }
  organism_db_name <- species_info$organism_db
  if (!require(organism_db_name, character.only = TRUE)) {
    stop("[02e_check_celltype_purity] organism DB 패키지가 설치되어 있지 않습니다: ", organism_db_name)
  }
  organism_db <- get(organism_db_name)

  keytype <- if (top_gene_id_type == "ENSEMBL") "ENSEMBL" else top_gene_id_type
  cat(paste0("Mapping ", keytype, " -> SYMBOL for marker matching...\n"))
  symbols <- suppressMessages(
    mapIds(organism_db, keys = rownames(vst_mat), column = "SYMBOL",
           keytype = keytype, multiVals = "first")
  )
  keep <- !is.na(symbols)
  expr_mat <- vst_mat[keep, , drop = FALSE]
  rownames(expr_mat) <- unname(symbols[keep])
} else {
  stop("[02e_check_celltype_purity] marker_id_type '", marker_id_type,
       "'은 top-level gene_id_type('", top_gene_id_type,
       "')과 다르지만 SYMBOL이 아니라서 자동 매핑을 지원하지 않습니다.")
}

if (marker_id_type == "SYMBOL") {
  # marker 목록은 항상 대문자(uppercase)로 읽어들이므로(사람 marker set을 마우스/랫
  # 데이터에 적용하는 경우처럼 종간 대소문자 표기 차이 — 예: "Aqp4" vs "AQP4" — 를
  # 흡수하기 위해) 발현 매트릭스 쪽도 동일하게 대문자로 맞춰서 매칭한다.
  rownames(expr_mat) <- toupper(rownames(expr_mat))
}
expr_mat <- aggregate_duplicate_rows(expr_mat)

cat(paste0("Expression matrix ready: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples\n\n"))

# --- 5. 샘플별 preranked fgsea (절대 VST 값 기준 — "이 샘플 안에서" 이 marker set이
# 발현 상위에 몰려 있는가를 검정. 위에서 미발현 유전자를 걸러냈으므로 배경이
# "발현되는 유전자들"로 한정되어 있음) ---
run_fgsea_for_sample <- function(sample_id) {
  ranks <- expr_mat[, sample_id]
  ranks <- sort(ranks, decreasing = TRUE)
  res <- suppressWarnings(
    fgsea(pathways = marker_list, stats = ranks,
          minSize = min_marker_genes, eps = 0)
  )
  res$sample_id <- sample_id
  res
}

cat("Running preranked fgsea per sample...\n")
all_results <- lapply(colnames(expr_mat), run_fgsea_for_sample)
all_results <- do.call(rbind, all_results)

skipped <- setdiff(names(marker_list), unique(all_results$pathway))
if (length(skipped) > 0) {
  cat(paste0("[WARN] 다음 cell type은 matched marker gene이 min_marker_genes(",
             min_marker_genes, ") 미만이라 건너뜀: ", paste(skipped, collapse = ", "), "\n"))
}

all_results$leadingEdge <- vapply(all_results$leadingEdge, function(x) paste(x, collapse = ","), character(1))
result_df <- as.data.frame(all_results)[, c("sample_id", "pathway", "NES", "pval", "padj", "size", "leadingEdge")]
colnames(result_df)[colnames(result_df) == "pathway"] <- "cell_type"
result_df <- result_df[order(result_df$sample_id, -result_df$NES), ]

csv_path <- file.path(out_dir, "celltype_purity_results.csv")
write.csv(result_df, csv_path, row.names = FALSE)
cat(paste0("Saved: ", csv_path, "\n"))

# --- 6. Heatmap (cell type x sample, NES, 유의 셀은 *) ---
# min_marker_genes 미만이라 fgsea에서 건너뛴(=skipped) cell type은 전 샘플에서 NA라
# pheatmap의 거리 계산(hclust)이 깨지므로 실제 결과가 있는 cell type만 사용한다.
cell_types <- intersect(names(marker_list), unique(result_df$cell_type))
sample_ids <- colnames(expr_mat)

nes_mat <- matrix(NA_real_, nrow = length(cell_types), ncol = length(sample_ids),
                   dimnames = list(cell_types, sample_ids))
sig_mat <- matrix("", nrow = length(cell_types), ncol = length(sample_ids),
                   dimnames = list(cell_types, sample_ids))

for (i in seq_len(nrow(result_df))) {
  r <- result_df[i, ]
  if (r$cell_type %in% cell_types && r$sample_id %in% sample_ids) {
    nes_mat[r$cell_type, r$sample_id] <- r$NES
    if (!is.na(r$padj) && r$padj < fdr_cutoff) sig_mat[r$cell_type, r$sample_id] <- "*"
  }
}

png_path <- file.path(out_dir, "celltype_purity_heatmap.png")
max_abs <- max(abs(nes_mat), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(c("#3375FF", "white", "#FF5733"))(100)

png(png_path, width = max(8, 1 + 0.5 * length(sample_ids)), height = max(4, 1 + 0.4 * length(cell_types)),
    units = "in", res = 150)
pheatmap(nes_mat, color = colors, breaks = breaks,
         display_numbers = sig_mat, number_color = "black",
         cluster_rows = TRUE, cluster_cols = TRUE,
         na_col = "grey90",
         main = "Cell-type marker enrichment (NES, * = padj < fdr_cutoff)")
dev.off()
cat(paste0("Saved: ", png_path, "\n"))

# --- 7. 요약 텍스트 ---
summary_path <- file.path(out_dir, "celltype_purity_summary.txt")
lines <- c(
  "=== Cell-type Purity Check Summary ===",
  "",
  "NOTE: marker set이 크고(수십~수백 유전자) 이 검정은 '발현되는 유전자들' 배경 대비",
  "preranked GSEA이므로, 실제 조직/배양이라면 여러 cell type이 동시에 padj<fdr_cutoff를",
  "통과할 수 있습니다(marker 자체가 어느 정도 baseline 발현을 갖는 경우가 흔함). 절대적인",
  "유의/비유의보다 아래 'Top by NES' 순위(상대적 강도)를 1차 판단 근거로 보는 것을 권장합니다.",
  ""
)

for (sid in sample_ids) {
  sample_results <- result_df[result_df$sample_id == sid, ]
  sample_results <- sample_results[order(-sample_results$NES), ]
  top_n <- min(5, nrow(sample_results))
  sig <- sample_results[!is.na(sample_results$padj) & sample_results$padj < fdr_cutoff, ]

  lines <- c(lines, paste0("[", sid, "]"))
  lines <- c(lines, "  Top by NES (regardless of significance):")
  for (i in seq_len(top_n)) {
    lines <- c(lines, sprintf("    %-20s NES=%.2f  padj=%.3g  (n_marker_genes=%d)",
                               sample_results$cell_type[i], sample_results$NES[i],
                               sample_results$padj[i], sample_results$size[i]))
  }

  lines <- c(lines, sprintf("  Significant (padj < %.3g):", fdr_cutoff))
  if (nrow(sig) == 0) {
    lines <- c(lines, "    None.")
  } else {
    for (i in seq_len(nrow(sig))) {
      lines <- c(lines, sprintf("    %-20s NES=%.2f  padj=%.3g  (n_marker_genes=%d)",
                                 sig$cell_type[i], sig$NES[i], sig$padj[i], sig$size[i]))
    }
  }

  if (!is.null(expected_type)) {
    # cell_type은 "{base} [{source}]" 형식이므로 출처 태그를 떼고 base만 비교
    base_label <- sub(" \\[[^]]*\\]$", "", sig$cell_type)
    is_expected <- tolower(base_label) == tolower(expected_type)
    expected_hit <- sig[is_expected, ]
    other_hits   <- sig[!is_expected, ]
    if (nrow(expected_hit) == 0) {
      lines <- c(lines, paste0("  [FLAG] Expected type '", expected_type,
                                "' is NOT significantly enriched in this sample."))
    }
    if (nrow(other_hits) > 0) {
      lines <- c(lines, paste0("  [FLAG] Non-expected cell type(s) significantly enriched: ",
                                paste(other_hits$cell_type, collapse = ", "),
                                " -- possible contamination."))
    }
  }
  lines <- c(lines, "")
}

writeLines(lines, summary_path)
cat(paste0("Saved: ", summary_path, "\n"))

cat("\n=== Cell-type Purity Check Complete ===\n")
