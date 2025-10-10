# 파일 경로: src/analysis/01_run_de_analysis.R

# 필요한 라이브러리 로드
library(here)
library(AnnotationDbi)

# --- [수정] 1. 공통 데이터 로드 및 객체 생성 ---
# DGE 분석에 공통적으로 필요한 데이터와 객체를 먼저 준비합니다.
count_data <- read.csv(here(config$count_data_path), row.names = 1)
meta_data <- read.csv(here(config$metadata_path), row.names = 1)
design_formula <- as.formula(config$de_analysis$design_formula)
# edgeR, limma-voom을 위한 design matrix를 미리 생성합니다.
design_matrix <- model.matrix(design_formula, data = meta_data)


# --- [교체] 2. config 설정에 따른 DGE 분석 방법 분기 ---
# 기존의 DESeq2 분석 실행 부분을 이 if/else if/else 블록으로 완전히 교체합니다.
dge_method <- config$de_analysis$method
cat(paste("Running DGE analysis using method:", dge_method, "\n"))

if (dge_method == "DESeq2") {
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = meta_data, design = design_formula)
  dds <- DESeq(dds)
  res <- results(dds, alpha = config$de_analysis$padj_cutoff)
  res_df <- as.data.frame(res)
  normalized_counts <- counts(dds, normalized = TRUE) # DESeq2용 정규화 카운트

} else if (dge_method == "edgeR") {
  library(edgeR)
  dge <- DGEList(counts = count_data)
  keep <- filterByExpr(dge, design_matrix)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design_matrix)
  fit <- glmQLFit(dge, design_matrix)
  qlf <- glmQLFTest(fit)
  res <- topTags(qlf, n = Inf)$table
  # [중요] 컬럼 이름 표준화
  res$log2FoldChange <- res$logFC
  res$padj <- res$FDR
  res$pvalue <- res$PValue
  res$baseMean <- 2^res$logCPM # baseMean과 유사한 값으로 logCPM 사용
  res_df <- as.data.frame(res)
  normalized_counts <- cpm(dge) # edgeR/limma용 정규화 카운트 (CPM)

} else if (dge_method == "limma-voom") {
  library(limma)
  library(edgeR) # DGEList 생성을 위해 필요
  dge <- DGEList(counts = count_data)
  keep <- filterByExpr(dge, design_matrix)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design_matrix, plot = FALSE)
  fit <- lmFit(v, design_matrix)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = ncol(design_matrix), number = Inf, sort.by = "P")
  # [중요] 컬럼 이름 표준화
  res$log2FoldChange <- res$logFC
  res$padj <- res$adj.P.Val
  res$pvalue <- res$P.Value
  res$baseMean <- res$AveExpr # baseMean과 유사한 값으로 AveExpr 사용
  res_df <- as.data.frame(res)
  normalized_counts <- cpm(dge) # edgeR/limma용 정규화 카운트 (CPM)

} else {
  stop("Invalid DGE method in config.yml. Choose from 'DESeq2', 'edgeR', 'limma-voom'.")
}


# --- [유지] 3. 공통 후처리 과정 ---
# 이 부분은 어떤 분석법을 사용했든 표준화된 'res_df'와 'normalized_counts'를 사용하므로
# 기존 코드를 거의 그대로 사용할 수 있습니다.

# DGE 분석 결과 정렬
res_df <- res_df[order(res_df$padj), ]

# 데이터-설정 일치 여부 검사
first_gene_id <- rownames(res_df)[1]
if (config$species == "human" && !startsWith(first_gene_id, "ENSG")) {
  stop("ERROR: Species is set to 'human' but gene IDs do not start with 'ENSG'.")
}
if (config$species == "mouse" && !startsWith(first_gene_id, "ENSMUSG")) {
  stop("ERROR: Species is set to 'mouse' but gene IDs do not start with 'ENSMUSG'.")
}

# Gene Symbol Annotation
species_info <- config$databases[[config$species]]
organism_db_name <- species_info$organism_db
if (!require(organism_db_name, character.only = TRUE)) {
  stop(paste("유전체 패키지", organism_db_name, "가 설치되지 않았습니다."))
}
res_df$symbol <- mapIds(get(organism_db_name),
                        keys = rownames(res_df),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# [수정] 표준화된 컬럼 이름만 선택하도록 수정
# lfcSE, stat 등은 DESeq2에만 있으므로 공통 컬럼만 선택합니다.
standard_cols <- c("symbol", "baseMean", "log2FoldChange", "pvalue", "padj")
res_df <- res_df[, standard_cols]


# 최종 결과 테이블 생성 및 저장
final_results_df <- merge(res_df, as.data.frame(normalized_counts), by = "row.names", sort = FALSE)
rownames(final_results_df) <- final_results_df$Row.names
final_results_df$Row.names <- NULL

stat_cols <- colnames(res_df) # res_df의 모든 컬럼을 사용하도록 수정
sample_cols <- colnames(normalized_counts)
final_results_df <- final_results_df[, c(stat_cols, sample_cols)]

# CSV 파일 저장
output_csv_path <- file.path(output_path, paste0("final_de_results.csv"))
write.csv(final_results_df, file = output_csv_path, row.names = TRUE)
print(paste("Combined results saved to:", output_csv_path))

# Excel 파일 저장 (선택 사항)
if (isTRUE(config$export$export_to_excel)) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Please install 'openxlsx' package to export to Excel.")
  }
  library(openxlsx)
  output_xlsx_path <- file.path(output_path, paste0("final_de_results.xlsx"))
  write.xlsx(final_results_df, file = output_xlsx_path, rowNames = TRUE, overwrite = TRUE)
  print(paste("Excel file also saved to:", output_xlsx_path))
}

# [추가] 6. 재현성을 위해 사용된 config 파일을 결과 폴더에 복사합니다.
file.copy(from = here("config.yml"), 
          to = file.path(output_path, "config_used.yml"),
          overwrite = TRUE)
print("Copied config file to output directory for reproducibility.")

print(paste(dge_method, "analysis and result saving process finished."))