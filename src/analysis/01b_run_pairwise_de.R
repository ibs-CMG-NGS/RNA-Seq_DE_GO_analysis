# 파일 경로: src/analysis/01b_run_pairwise_de.R
# Pairwise DE analysis
# 사용법: Rscript 01b_run_pairwise_de.R [config_path] [compare_group] [base_group] [output_dir]

# --- 1. Setup: Load config and libraries ---
suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(openxlsx)
})

# --- 2. Get Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 01b_run_pairwise_de.R [config_path] [compare_group] [base_group] [output_dir]")
}
config_path <- args[1]
compare_group <- args[2] # 예: "treated"
base_group <- args[3]    # 예: "control"
output_dir <- args[4]    # 예: "output/pairwise/treated_vs_control"

# --- 3. Load Config and Data ---
config <- yaml.load_file(config_path)
counts <- read.csv(here(config$count_data_path), row.names = 1)
meta <- read.csv(here(config$metadata_path), row.names = 1)
group_var <- config$de_analysis$group_variable
meta[[group_var]] <- as.factor(meta[[group_var]])
design_formula <- as.formula(config$de_analysis$design_formula)

# --- 4. Run Pairwise Test based on method ---
dge_method <- config$de_analysis$method
cat(paste("Running Pairwise DE:", compare_group, "vs", base_group, "using", dge_method, "\n"))

if (dge_method == "DESeq2") {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
  dds <- DESeq(dds) # 표준 Wald test 실행
  res <- results(dds, contrast = c(group_var, compare_group, base_group), 
                 alpha = config$de_analysis$padj_cutoff)
  res_df <- as.data.frame(res[order(res$padj), ])
  normalized_counts <- counts(dds, normalized = TRUE)

} else if (dge_method == "edgeR") {
  design_matrix <- model.matrix(design_formula, data = meta)
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, design_matrix)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design_matrix)
  fit <- glmQLFit(dge, design_matrix)
  
  # contrast 생성
  contrast_str <- paste0(group_var, compare_group, " - ", group_var, base_group)
  contrast_vec <- makeContrasts(contrasts = contrast_str, levels = colnames(design_matrix))
  
  qlf <- glmQLFTest(fit, contrast = contrast_vec)
  res <- topTags(qlf, n = Inf)$table
  res$log2FoldChange <- res$logFC; res$padj <- res$FDR; res$pvalue <- res$PValue; res$baseMean <- 2^res$logCPM
  res_df <- as.data.frame(res)
  normalized_counts <- cpm(dge)

} else if (dge_method == "limma-voom") {
  design_matrix <- model.matrix(design_formula, data = meta)
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, design_matrix)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design_matrix, plot = FALSE)
  fit <- lmFit(v, design_matrix)
  fit <- eBayes(fit)
  
  # contrast 생성
  contrast_str <- paste0(group_var, compare_group, " - ", group_var, base_group)
  contrast_vec <- makeContrasts(contrasts = contrast_str, levels = colnames(design_matrix))
  
  fit_contrast <- contrasts.fit(fit, contrast_vec)
  fit_contrast <- eBayes(fit_contrast)
  
  res <- topTable(fit_contrast, number = Inf, sort.by = "P")
  res$log2FoldChange <- res$logFC; res$padj <- res$adj.P.Val; res$pvalue <- res$P.Value; res$baseMean <- res$AveExpr
  res_df <- as.data.frame(res)
  normalized_counts <- cpm(dge)
  
} else {
  stop("Invalid DGE method.")
}

# --- 5. Annotation (기존 01_...R 스크립트의 로직 재사용) ---
cat("Annotating results...\n")
species_info <- config$databases[[config$species]]
organism_db_name <- species_info$organism_db
if (!require(organism_db_name, character.only = TRUE)) {
  stop(paste("Genome package", organism_db_name, "is not installed."))
}
res_df$symbol <- mapIds(get(organism_db_name),
                        keys = rownames(res_df),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

standard_cols <- c("symbol", "baseMean", "log2FoldChange", "pvalue", "padj")
res_df <- res_df[, c(standard_cols, setdiff(colnames(res_df), standard_cols))]

# --- 6. Save Results (기존 01_...R 스크립트의 로직 재사용) ---
final_results_df <- merge(res_df, as.data.frame(normalized_counts), by = "row.names", sort = FALSE)
rownames(final_results_df) <- final_results_df$Row.names
final_results_df$Row.names <- NULL

# ... (컬럼 순서 재정렬) ...

# 고정된 이름으로 CSV/Excel 저장 (폴더가 이미 쌍별로 분리됨)
output_csv_path <- file.path(output_dir, "final_de_results.csv")
write.csv(final_results_df, file = output_csv_path, row.names = TRUE)
cat(paste("Pairwise results saved to:", output_csv_path, "\n"))

if (isTRUE(config$export$export_to_excel)) {
  output_xlsx_path <- file.path(output_dir, "final_de_results.xlsx")
  write.xlsx(final_results_df, file = output_xlsx_path, rowNames = TRUE, overwrite = TRUE)
  cat(paste("Excel file also saved to:", output_xlsx_path, "\n"))
}

# config 파일 복사
file.copy(from = config_path, to = file.path(output_dir, "config_used.yml"), overwrite = TRUE)
cat("Pairwise analysis finished.\n")