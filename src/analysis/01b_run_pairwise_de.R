# 파일 경로: src/analysis/01b_run_pairwise_de.R
# Pairwise DE analysis
# Usage: Rscript 01b_run_pairwise_de.R [config_path] [compare_group] [base_group] [output_dir]

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(openxlsx)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 01b_run_pairwise_de.R [config_path] [compare_group] [base_group] [output_dir]")
}
config_path <- args[1]
compare_group <- args[2]
base_group <- args[3]
output_dir <- args[4]

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. 설정 및 데이터 로드 ---
config <- yaml.load_file(config_path)
# check.names=FALSE로 설정하여 샘플 이름(컬럼명)이 R에 의해 변형되는 것을 방지 (예: - 가 . 으로 바뀌는 문제)
counts <- read.csv(here(config$count_data_path), row.names = 1, check.names = FALSE)
meta <- read.csv(here(config$metadata_path), row.names = 1)

# 그룹 변수 확인
group_var <- config$de_analysis$group_variable
if (!group_var %in% colnames(meta)) {
    stop(paste("Group variable", group_var, "not found in metadata."))
}

# --- 3. 해당 비교 쌍에 맞는 샘플 필터링 및 정렬 ---
cat(paste("Filtering samples for comparison:", compare_group, "vs", base_group, "\n"))

# Normalization strategy 확인
norm_strategy <- config$de_analysis$advanced_options$pairwise_normalization
if (is.null(norm_strategy)) {
  norm_strategy <- "subset"  # 기본값
}
cat(paste("Normalization strategy:", norm_strategy, "\n"))

# 1) 메타데이터에서 두 그룹에 속하는 샘플만 선택
target_samples <- rownames(meta)[meta[[group_var]] %in% c(compare_group, base_group)]
if (length(target_samples) == 0) {
    stop("No samples found for the specified groups. Check group names in config vs metadata.")
}

# 2) Pre-filtering 및 데이터 준비는 normalization strategy에 따라 다르게 처리
if (norm_strategy == "global") {
  # Global normalization: 전체 샘플로 정규화
  cat("Using GLOBAL normalization: all samples will be used for size factor calculation.\n")
  
  # Ensure sample names match between count data and metadata (전체 데이터)
  # counts may contain more samples than metadata (shared counts file across subgroups)
  count_samples <- colnames(counts)
  meta_samples <- rownames(meta)

  missing_in_counts <- meta_samples[!meta_samples %in% count_samples]
  if (length(missing_in_counts) > 0) {
    stop(paste("Metadata samples not found in count data:", paste(missing_in_counts, collapse=", ")))
  }

  # Subset counts to only the samples present in metadata
  counts <- counts[, meta_samples, drop = FALSE]
  meta <- meta[meta_samples, , drop = FALSE]
  cat(paste("Aligned", ncol(counts), "samples between count data and metadata\n"))
  
  # Pre-filtering on full dataset
  # 기준: 샘플당 count >= min_count 인 샘플이 min_samples개 이상 (샘플 수에 독립적)
  pf_min_count       <- config$de_analysis$advanced_options$prefilter_min_count %||% 1
  pf_min_samples_cfg <- config$de_analysis$advanced_options$prefilter_min_samples %||% "auto"
  if (identical(pf_min_samples_cfg, "auto")) {
    pf_min_samples <- as.integer(min(table(meta[[group_var]])))
  } else {
    pf_min_samples <- as.integer(pf_min_samples_cfg)
  }
  keep_genes <- rowSums(counts >= pf_min_count) >= pf_min_samples
  counts     <- counts[keep_genes, ]
  cat(paste("Pre-filtering (global): kept", sum(keep_genes), "of", length(keep_genes),
            "genes (count >=", pf_min_count, "in >=", pf_min_samples, "samples)\n"))
  
  # Factor 설정 (전체 데이터)
  meta[[group_var]] <- as.factor(meta[[group_var]])
  
  # counts_full과 meta_full 저장 (나중에 사용)
  counts_full <- counts
  meta_full <- meta
  
  # Subset for this comparison (필터링된 데이터에서)
  meta_subset <- meta[target_samples, , drop = FALSE]
  counts_subset <- counts[, target_samples, drop = FALSE]
  counts_subset <- counts_subset[, rownames(meta_subset)]
  
  # Factor 레벨 재설정 (비교 대상 그룹만)
  meta_subset[[group_var]] <- factor(meta_subset[[group_var]], levels = c(base_group, compare_group))
  
} else {
  # Subset normalization: 각 비교마다 해당 두 그룹만으로 정규화
  cat("Using SUBSET normalization: only samples in this comparison will be used.\n")
  
  # 메타데이터와 카운트 데이터를 필터링
  meta_subset <- meta[target_samples, , drop = FALSE]
  counts_subset <- counts[, target_samples, drop = FALSE]
  
  # 카운트 데이터의 컬럼 순서를 메타데이터의 행 순서와 완벽하게 일치시킴
  counts_subset <- counts_subset[, rownames(meta_subset)]
  
  # Pre-filtering (subset data)
  # 기준: 샘플당 count >= min_count 인 샘플이 min_samples개 이상 (샘플 수에 독립적)
  pf_min_count       <- config$de_analysis$advanced_options$prefilter_min_count %||% 1
  pf_min_samples_cfg <- config$de_analysis$advanced_options$prefilter_min_samples %||% "auto"
  if (identical(pf_min_samples_cfg, "auto")) {
    pf_min_samples <- as.integer(min(table(meta_subset[[group_var]])))
  } else {
    pf_min_samples <- as.integer(pf_min_samples_cfg)
  }
  keep_genes    <- rowSums(counts_subset >= pf_min_count) >= pf_min_samples
  counts_subset <- counts_subset[keep_genes, ]
  cat(paste("Pre-filtering (subset): kept", sum(keep_genes), "of", length(keep_genes),
            "genes (count >=", pf_min_count, "in >=", pf_min_samples, "samples)\n"))
  
  # Factor 레벨 재설정
  meta_subset[[group_var]] <- factor(meta_subset[[group_var]], levels = c(base_group, compare_group))
  
  counts_full <- NULL
  meta_full <- NULL
}

# --- 4. 분석 실행 ---
design_formula <- as.formula(config$de_analysis$design_formula)
dge_method <- config$de_analysis$method
cat(paste("Running Pairwise DE:", compare_group, "vs", base_group, "using", dge_method, "\n"))

if (dge_method == "DESeq2") {
  if (norm_strategy == "global") {
    # Global normalization approach
    cat("DESeq2 with global normalization:\n")
    cat("  Step 1: Creating DESeqDataSet with ALL samples...\n")
    
    # 전체 데이터로 DESeqDataSet 생성
    dds_full <- DESeqDataSetFromMatrix(countData = counts_full, 
                                       colData = meta_full, 
                                       design = design_formula)
    
    # 전체 데이터로 size factors 계산 (정규화)
    dds_full <- estimateSizeFactors(dds_full)
    cat(paste("  Step 2: Calculated size factors from", ncol(dds_full), "samples\n"))
    
    # 이제 subset 데이터만 추출하여 분산 추정 및 DE 분석
    cat("  Step 3: Subsetting to comparison samples and running DE...\n")
    dds <- dds_full[, target_samples]
    
    # CRITICAL: subset 후 factor 레벨을 정리하여 사용하지 않는 레벨 제거
    # 이렇게 하지 않으면 model matrix가 full rank가 아니게 됨
    colData(dds)[[group_var]] <- factor(colData(dds)[[group_var]], 
                                         levels = c(base_group, compare_group))
    
    # Design formula를 다시 설정 (정리된 factor 레벨 반영)
    design(dds) <- design_formula
    
    # Size factors는 유지됨! (global normalization의 핵심)
    # 분산 추정 및 DE 분석 실행
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    
    # Results 추출
    res <- results(dds, contrast = c(group_var, compare_group, base_group), 
                   alpha = config$de_analysis$padj_cutoff)
    res_df <- as.data.frame(res[order(res$padj), ])
    
    # Normalized counts (global size factors 사용)
    normalized_counts <- counts(dds, normalized = TRUE)
    
  } else {
    # Subset normalization approach (기존 방식)
    cat("DESeq2 with subset normalization:\n")
    dds <- DESeqDataSetFromMatrix(countData = counts_subset, 
                                  colData = meta_subset, 
                                  design = design_formula)
    dds <- DESeq(dds)
    
    # contrast 명시
    res <- results(dds, contrast = c(group_var, compare_group, base_group), 
                   alpha = config$de_analysis$padj_cutoff)
    res_df <- as.data.frame(res[order(res$padj), ])
    normalized_counts <- counts(dds, normalized = TRUE)
  }

} else if (dge_method == "edgeR") {
  if (norm_strategy == "global") {
    cat("edgeR with global normalization:\n")
    cat("  Step 1: Calculating normalization factors from ALL samples...\n")
    
    # 전체 데이터로 정규화
    design_matrix_full <- model.matrix(design_formula, data = meta_full)
    dge_full <- DGEList(counts = counts_full)
    dge_full <- calcNormFactors(dge_full)
    
    cat("  Step 2: Subsetting to comparison samples...\n")
    # Subset으로 좁히기 (normalization factors 유지)
    dge <- dge_full[, target_samples]
    
    # Subset에 대한 design matrix
    design_matrix <- model.matrix(design_formula, data = meta_subset)
    
    # 분산 추정 및 DE 분석
    dge <- estimateDisp(dge, design_matrix)
    fit <- glmQLFit(dge, design_matrix)
    
  } else {
    cat("edgeR with subset normalization:\n")
    design_matrix <- model.matrix(design_formula, data = meta_subset)
    dge <- DGEList(counts = counts_subset)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design_matrix)
    fit <- glmQLFit(dge, design_matrix)
  }
  
  # Contrast 생성 및 테스트 (공통)
  contrast_str <- paste0(group_var, compare_group, " - ", group_var, base_group)
  contrast_vec <- makeContrasts(contrasts = contrast_str, levels = colnames(design_matrix))
  
  qlf <- glmQLFTest(fit, contrast = contrast_vec)
  res <- topTags(qlf, n = Inf)$table
  res$log2FoldChange <- res$logFC; res$padj <- res$FDR; res$pvalue <- res$PValue; res$baseMean <- 2^res$logCPM
  res_df <- as.data.frame(res)
  normalized_counts <- cpm(dge)

} else if (dge_method == "limma-voom") {
  if (norm_strategy == "global") {
    cat("limma-voom with global normalization:\n")
    cat("  Step 1: Normalizing ALL samples...\n")
    
    # 전체 데이터로 정규화
    design_matrix_full <- model.matrix(design_formula, data = meta_full)
    dge_full <- DGEList(counts = counts_full)
    keep_full <- filterByExpr(dge_full, design_matrix_full)
    dge_full <- dge_full[keep_full, , keep.lib.sizes=FALSE]
    dge_full <- calcNormFactors(dge_full)
    
    cat("  Step 2: Subsetting and running voom transformation...\n")
    # Subset (normalization factors 유지)
    dge <- dge_full[, target_samples]
    design_matrix <- model.matrix(design_formula, data = meta_subset)
    
    v <- voom(dge, design_matrix, plot = FALSE)
    fit <- lmFit(v, design_matrix)
    fit <- eBayes(fit)
    
  } else {
    cat("limma-voom with subset normalization:\n")
    design_matrix <- model.matrix(design_formula, data = meta_subset)
    dge <- DGEList(counts = counts_subset)
    keep <- filterByExpr(dge, design_matrix)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)
    v <- voom(dge, design_matrix, plot = FALSE)
    fit <- lmFit(v, design_matrix)
    fit <- eBayes(fit)
  }
  
  # Contrast 생성 및 테스트 (공통)
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

# --- 5. Annotation (Entrez ID 지원 추가) ---
cat("Annotating results...\n")
species_info <- config$databases[[config$species]]
organism_db_name <- species_info$organism_db
if (!require(organism_db_name, character.only = TRUE)) {
  stop(paste("Genome package", organism_db_name, "is not installed."))
}

valid_ids <- rownames(res_df)
organism_db <- get(organism_db_name)

# ID 타입 자동 감지 로직
first_id <- as.character(valid_ids[1])
is_ensembl <- grepl("^ENS", first_id) # Human: ENSG, Mouse: ENSMUSG
is_entrez <- grepl("^[0-9]+$", first_id) # 숫자로만 구성된 경우 Entrez ID로 간주

# 매핑 실행
if (is_ensembl) {
    cat("Detected ENSEMBL IDs. Attempting to map to SYMBOL...\n")
    tryCatch({
        mapped_symbols <- mapIds(organism_db,
                                keys = valid_ids,
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "first")
        res_df$symbol <- mapped_symbols
    }, error = function(e) {
        cat(paste("Warning: Annotation failed (ENSEMBL):", e$message, "\n"))
        res_df$symbol <- rownames(res_df)
    })
} else if (is_entrez) {
    cat("Detected Entrez IDs (numeric). Attempting to map to SYMBOL...\n")
    tryCatch({
        mapped_symbols <- mapIds(organism_db,
                                keys = valid_ids,
                                column = "SYMBOL",
                                keytype = "ENTREZID", # Keytype을 ENTREZID로 설정
                                multiVals = "first")
        res_df$symbol <- mapped_symbols
    }, error = function(e) {
        cat(paste("Warning: Annotation failed (ENTREZID):", e$message, "\n"))
        res_df$symbol <- rownames(res_df)
    })
} else {
    cat(paste("Warning: Gene IDs do not look like ENSEMBL or Entrez IDs (e.g.,", first_id, "). Skipping annotation assuming they might already be symbols.\n"))
    res_df$symbol <- rownames(res_df)
}

# NA symbol을 gene_id로 대체
if ("symbol" %in% colnames(res_df)) {
  na_count <- sum(is.na(res_df$symbol))
  if (na_count > 0) {
    res_df$symbol[is.na(res_df$symbol)] <- rownames(res_df)[is.na(res_df$symbol)]
    cat(paste("Symbol fallback: filled", na_count, "unmapped gene(s) with gene_id\n"))
  }
}

# 컬럼 정리
if ("symbol" %in% colnames(res_df)) {
    standard_cols <- c("symbol", "baseMean", "log2FoldChange", "pvalue", "padj")
} else {
    standard_cols <- c("baseMean", "log2FoldChange", "pvalue", "padj")
}

cols_to_keep <- intersect(standard_cols, colnames(res_df))
res_df <- res_df[, c(cols_to_keep, setdiff(colnames(res_df), cols_to_keep))]

# --- 6. Save Results ---
# 정규화 카운트와 병합
final_results_df <- merge(res_df, as.data.frame(normalized_counts), by = "row.names", sort = FALSE)
rownames(final_results_df) <- final_results_df$Row.names
final_results_df$Row.names <- NULL

# Reorder columns: DE stats first, then normalized counts grouped by condition and sorted
de_stat_cols <- intersect(c("symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), 
                          colnames(final_results_df))
norm_count_cols <- setdiff(colnames(final_results_df), de_stat_cols)

# Group and sort normalized count columns by condition
# Use the original metadata to get group assignments for target samples
sample_order <- character()

# For each group (base first, then compare), add sorted samples
for (grp in c(base_group, compare_group)) {
  # Get samples belonging to this group from target_samples
  grp_samples <- target_samples[meta[target_samples, group_var] == grp]
  grp_samples_sorted <- sort(grp_samples)  # Alphabetical order within group
  sample_order <- c(sample_order, grp_samples_sorted)
}

# Reorder columns: DE stats + ordered normalized counts
final_results_df <- final_results_df[, c(de_stat_cols, sample_order)]

cat(paste("Column order: DE stats followed by samples grouped by condition (", 
          base_group, "then", compare_group, ") and sorted alphabetically\n"))

# Save main DE results with normalized counts
output_csv_path <- file.path(output_dir, "final_de_results.csv")
write.csv(final_results_df, file = output_csv_path, row.names = TRUE)
cat(paste("Pairwise results saved to:", output_csv_path, "\n"))

# Save normalized counts separately for easier access
# Also reorder columns for normalized counts
normalized_counts_ordered <- normalized_counts[, sample_order]
normalized_counts_path <- file.path(output_dir, "normalized_counts.csv")
write.csv(normalized_counts_ordered, file = normalized_counts_path, row.names = TRUE)
cat(paste("Normalized counts saved to:", normalized_counts_path, "\n"))

if (isTRUE(config$export$export_to_excel)) {
  output_xlsx_path <- file.path(output_dir, "final_de_results.xlsx")
  
  # Create a workbook with multiple sheets
  wb <- createWorkbook()
  
  # Helper: prepend rownames as gene_id column
  with_gene_id <- function(df) cbind(gene_id = rownames(df), df)

  # Sheet 1: DE results with normalized counts
  addWorksheet(wb, "DE_Results")
  writeData(wb, "DE_Results", with_gene_id(final_results_df), rowNames = FALSE)

  # Sheet 2: Normalized counts only (ordered)
  addWorksheet(wb, "Normalized_Counts")
  writeData(wb, "Normalized_Counts", with_gene_id(normalized_counts_ordered), rowNames = FALSE)

  # Sheet 3: Significant genes only
  sig_genes <- final_results_df[!is.na(final_results_df$padj) &
                                 final_results_df$padj < config$de_analysis$padj_cutoff &
                                 abs(final_results_df$log2FoldChange) > config$de_analysis$log2fc_cutoff, ]
  addWorksheet(wb, "Significant_Genes")
  writeData(wb, "Significant_Genes", with_gene_id(sig_genes), rowNames = FALSE)
  
  saveWorkbook(wb, output_xlsx_path, overwrite = TRUE)
  cat(paste("Excel file with multiple sheets saved to:", output_xlsx_path, "\n"))
}

file.copy(from = config_path, to = file.path(output_dir, "config_used.yml"), overwrite = TRUE)
cat("Pairwise analysis finished.\n")