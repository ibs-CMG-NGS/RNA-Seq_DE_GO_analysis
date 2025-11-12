# 파일 경로: src/analysis/01a_run_omnibus_test.R
# Omnibus (ANOVA-like) test
# 사용법: Rscript 01a_run_omnibus_test.R [config_path] [output_csv_path]

# --- 1. Setup: Load config and libraries ---
suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(DESeq2)
  library(edgeR)
  library(limma)
})

# --- 2. Get Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 01a_run_omnibus_test.R [config_path] [output_csv_path]")
}
config_path <- args[1]
output_csv_path <- args[2]

# --- 3. Load Config and Data ---
config <- yaml.load_file(config_path)
counts <- read.csv(here(config$count_data_path), row.names = 1)
meta <- read.csv(here(config$metadata_path), row.names = 1)

# 그룹 변수 factor로 변환
group_var <- config$de_analysis$group_variable
meta[[group_var]] <- as.factor(meta[[group_var]])

# 디자인 공식
design_formula <- as.formula(config$de_analysis$design_formula)

# --- 4. Run Omnibus Test based on method ---
dge_method <- config$de_analysis$method
cat(paste("Running Omnibus test using method:", dge_method, "\n"))

if (dge_method == "DESeq2") {
  # DESeq2: Likelihood Ratio Test (LRT)
  # 'reduced' 모델은 그룹 변수만 제거한 모델입니다.
  reduced_formula <- as.formula(gsub(paste0("\\+ *", group_var), "", as.character(design_formula)[2]))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
  dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
  res <- results(dds)
  res_df <- as.data.frame(res[order(res$padj), ])

} else if (dge_method == "edgeR") {
  # edgeR: Quasi-Likelihood F-test
  design_matrix <- model.matrix(design_formula, data = meta)
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, design_matrix)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design_matrix)
  fit <- glmQLFit(dge, design_matrix)
  
  # 모든 그룹 계수(coefficients)에 대해 F-test 수행
  coef_indices <- (ncol(fit$design) - nlevels(meta[[group_var]]) + 2):ncol(fit$design)
  qlf <- glmQLFTest(fit, coef = coef_indices)
  res <- topTags(qlf, n = Inf)$table
  res_df <- as.data.frame(res)
  
} else if (dge_method == "limma-voom") {
  # limma: F-test
  design_matrix <- model.matrix(design_formula, data = meta)
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, design_matrix)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design_matrix, plot = FALSE)
  fit <- lmFit(v, design_matrix)
  fit <- eBayes(fit)

  # 모든 그룹 계수에 대해 F-test 수행
  coef_indices <- (ncol(fit$design) - nlevels(meta[[group_var]]) + 2):ncol(fit$design)
  res <- topTable(fit, coef = coef_indices, number = Inf, sort.by = "F")
  res_df <- as.data.frame(res)
  
} else {
  stop("Invalid DGE method.")
}

# --- 5. Save Results ---
write.csv(res_df, output_csv_path, row.names = TRUE)
cat(paste("Omnibus test results saved to:", output_csv_path, "\n"))