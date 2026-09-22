# 파일 경로: src/analysis/02f_run_bretigea_deconvolution.R
# 목적: BRETIGEA(BRain cEll Type specific gene Expression Analysis) 패키지를 이용한
#   뇌세포타입 deconvolution 기반 순도(purity) 체크 (모든 샘플 대상)
#   02e_check_celltype_purity.R(preranked fgsea, 절대 발현량 기준 단일 샘플 검정)와
#   상호보완적인 방법 — BRETIGEA는 각 cell type의 top marker gene을 행 중심화(centering)한
#   뒤 SVD 1st component로 "surrogate proportion variable(SPV)"을 뽑는다. 이 값은 이
#   프로젝트 샘플 코호트 내에서의 상대적 크기(임의 단위, 코호트 평균이 0)이며, 절대적인
#   세포 비율(%)도 아니고 유의성 검정(p-value)도 제공하지 않는다 — 샘플 간 상대 비교/
#   이상치 탐지에 적합하고, 02e의 "이 샘플이 절대적으로 이 타입처럼 보이는가"와는
#   다른 질문에 답한다.
# 사용법: Rscript 02f_run_bretigea_deconvolution.R --config config.yml --output_dir output/qc_plots

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(optparse)
  library(DESeq2)
  library(BRETIGEA)
  library(pheatmap)
  library(AnnotationDbi)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- 1. Setup: Load libraries and parse arguments ---
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yml",
              help = "Path to the config YAML file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character",
              help = "Directory to save deconvolution results", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Output directory must be supplied.", call. = FALSE)
}

config <- yaml.load_file(opt$config)
bt_cfg <- config$qc_plots$celltype_deconvolution

if (is.null(bt_cfg) || !isTRUE(bt_cfg$enabled)) {
  cat("[02f_run_bretigea_deconvolution] qc_plots.celltype_deconvolution.enabled이 false입니다. 종료.\n")
  quit(status = 0)
}

out_dir <- file.path(opt$output_dir, "celltype_deconvolution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_marker  <- bt_cfg$nMarker %||% 50
celltypes <- bt_cfg$celltypes %||% c("ast", "end", "mic", "neu", "oli", "opc")
method    <- bt_cfg$method %||% "SVD"

# BRETIGEA는 human/mouse 전용 marker set만 있음 — 그 외 종은 "combined"(human-style
# uppercase 심볼, McKenzie et al. 통합 세트)로 대체하고 대소문자를 강제 정규화한다.
top_species <- tolower(config$species %||% "")
if (!is.null(bt_cfg$species)) {
  bretigea_species <- bt_cfg$species
} else if (top_species == "human") {
  bretigea_species <- "human"
} else if (top_species == "mouse") {
  bretigea_species <- "mouse"
} else {
  bretigea_species <- "combined"
  cat(paste0("[WARN] species '", top_species, "'에 대한 BRETIGEA 전용 marker set이 없어 ",
             "'combined'(human-style) marker set을 대소문자 정규화 후 사용합니다.\n"))
}

cat("\n=== BRETIGEA Deconvolution (SPV via ", method, ") ===\n", sep = "")
cat(paste0("BRETIGEA species: ", bretigea_species, "\n"))
cat(paste0("Cell types       : ", paste(celltypes, collapse = ", "), "\n"))
cat(paste0("nMarker per type : ", n_marker, "\n\n"))

# --- 2. DESeq2 dds 생성 + VST ---
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
vst_mat <- assay(vsd)

# 02e와 동일한 이유로 사실상 미발현 유전자를 제거 — SVD는 특히 분산이 큰(=미발현 유전자
# floor 근처의 인위적 분산 포함) 유전자에 민감하므로 이 필터가 더 중요하다.
detectable <- rowSums(counts(dds) >= 1) >= ceiling(ncol(dds) / 2)
cat(paste0("Filtering to detectably-expressed genes: ", sum(detectable), " / ", length(detectable), " genes kept\n"))
vst_mat <- vst_mat[detectable, , drop = FALSE]

# --- 3. gene_id_type -> SYMBOL 매핑 ---
top_gene_id_type <- toupper(config$gene_id_type %||% "ENSEMBL")

aggregate_duplicate_rows <- function(mat) {
  if (!any(duplicated(rownames(mat)))) return(mat)
  sums <- rowsum(mat, group = rownames(mat))
  counts <- as.vector(table(rownames(mat))[rownames(sums)])
  sums / counts
}

if (top_gene_id_type == "SYMBOL") {
  expr_mat <- vst_mat
} else {
  species_info <- config$databases[[config$species]]
  if (is.null(species_info) || is.null(species_info$organism_db)) {
    stop("[02f_run_bretigea_deconvolution] databases[", config$species, "].organism_db 설정이 없습니다.")
  }
  organism_db_name <- species_info$organism_db
  if (!require(organism_db_name, character.only = TRUE)) {
    stop("[02f_run_bretigea_deconvolution] organism DB 패키지가 설치되어 있지 않습니다: ", organism_db_name)
  }
  organism_db <- get(organism_db_name)

  keytype <- if (top_gene_id_type == "ENSEMBL") "ENSEMBL" else top_gene_id_type
  cat(paste0("Mapping ", keytype, " -> SYMBOL...\n"))
  symbols <- suppressMessages(
    mapIds(organism_db, keys = rownames(vst_mat), column = "SYMBOL",
           keytype = keytype, multiVals = "first")
  )
  keep <- !is.na(symbols)
  expr_mat <- vst_mat[keep, , drop = FALSE]
  rownames(expr_mat) <- unname(symbols[keep])
}

if (bretigea_species == "combined") {
  # combined marker set은 human-style uppercase 심볼 — mouse/rat 등 종간 대소문자
  # 표기 차이(예: "Aqp4" vs "AQP4")를 흡수하기 위해 발현 매트릭스도 대문자로 정규화.
  rownames(expr_mat) <- toupper(rownames(expr_mat))
}
expr_mat <- aggregate_duplicate_rows(expr_mat)

cat(paste0("Expression matrix ready: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples\n\n"))

# --- 4. BRETIGEA::brainCells() 실행 ---
# brainCells()/findCells() 내부에 print() 호출이 다수 있어(matched marker 목록,
# correlation 등) capture.output으로 표준출력을 삼키고 반환값만 사용한다.
cat("Running BRETIGEA::brainCells()...\n")
invisible(capture.output(
  spv_mat <- brainCells(inputMat = as.matrix(expr_mat), nMarker = n_marker,
                         species = bretigea_species, celltypes = celltypes,
                         method = method, scale = TRUE)
))
# spv_mat: sample(row) x celltype(col)

result_df <- as.data.frame(spv_mat)
result_df$sample_id <- rownames(result_df)
result_df <- result_df[, c("sample_id", celltypes)]

csv_path <- file.path(out_dir, "celltype_deconvolution_results.csv")
write.csv(result_df, csv_path, row.names = FALSE)
cat(paste0("Saved: ", csv_path, "\n"))

# --- 5. Heatmap (cell type x sample, SPV) ---
heat_mat <- t(spv_mat)   # celltype(row) x sample(col)

png_path <- file.path(out_dir, "celltype_deconvolution_heatmap.png")
max_abs <- max(abs(heat_mat), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(c("#3375FF", "white", "#FF5733"))(100)

png(png_path, width = max(8, 1 + 0.5 * ncol(heat_mat)), height = max(4, 1 + 0.4 * nrow(heat_mat)),
    units = "in", res = 150)
pheatmap(heat_mat, color = colors, breaks = breaks,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "BRETIGEA surrogate proportion variable (SPV, relative to this sample cohort)")
dev.off()
cat(paste0("Saved: ", png_path, "\n"))

# --- 6. 요약 텍스트 ---
summary_path <- file.path(out_dir, "celltype_deconvolution_summary.txt")
lines <- c(
  "=== BRETIGEA Cell-type Deconvolution Summary ===",
  "",
  "NOTE: 아래 SPV(surrogate proportion variable)는 이 프로젝트 샘플 코호트 안에서만",
  "의미 있는 상대값입니다(코호트 평균이 0이 되도록 행중심화한 marker 유전자의 SVD",
  "1st component). 절대적인 세포 비율(%)이 아니고 BRETIGEA 자체는 유의성 검정을",
  "제공하지 않습니다 — '이 샘플이 코호트 내 다른 샘플에 비해 상대적으로 이 타입",
  "signature가 강한가/약한가'로 해석하세요. 02e_check_celltype_purity.R(fgsea,",
  "절대 발현량 기준)과 함께 보는 것을 권장합니다.",
  ""
)

for (sid in result_df$sample_id) {
  row <- result_df[result_df$sample_id == sid, celltypes, drop = TRUE]
  ord <- order(-unlist(row))
  lines <- c(lines, paste0("[", sid, "]"))
  for (ct in celltypes[ord]) {
    lines <- c(lines, sprintf("    %-6s SPV=%.3f", ct, row[[ct]]))
  }
  lines <- c(lines, "")
}

writeLines(lines, summary_path)
cat(paste0("Saved: ", summary_path, "\n"))

cat("\n=== BRETIGEA Deconvolution Complete ===\n")
