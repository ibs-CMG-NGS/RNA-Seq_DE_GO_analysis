# 파일 경로: src/analysis/02_generate_plots.R
# 사용법:
# 1. (PCA) Rscript 02_generate_plots.R --config config.yml --task pca --output_file output/global_pca_plot.png
# 2. (Volcano) Rscript 02_generate_plots.R --config config.yml --task volcano --input_file output/pairwise/X/final_de_results.csv --output_file output/pairwise/X/volcano_plot.png

# --- 1. Setup: Load libraries and parse arguments ---
suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(optparse)
  library(ggplot2)
  library(ggrepel)
})

# Argument parsing
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yml", 
              help = "Path to the config YAML file", metavar = "character"),
  make_option(c("-t", "--task"), type = "character", 
              help = "Task to perform: 'pca' or 'volcano'", metavar = "character"),
  make_option(c("-i", "--input_file"), type = "character", default = NULL, 
              help = "Path to input DE results file (required for volcano)", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", 
              help = "Path to the output PNG file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$task) || is.null(opt$output_file)) {
  print_help(opt_parser)
  stop("Task and output_file arguments must be supplied.", call. = FALSE)
}

# Load config
config <- yaml.load_file(opt$config)

# --- 2. Execute selected task ---

if (opt$task == "pca") {
  # --- Task 2a: Generate Global PCA Plot ---
  cat(paste("\nGenerating Global PCA plot...\n"))
  
  # Load data
  count_data <- read.csv(here(config$count_data_path), row.names = 1, check.names = FALSE)
  meta_data <- read.csv(here(config$metadata_path), row.names = 1)
  
  # Ensure sample names match between count data and metadata
  count_samples <- colnames(count_data)
  meta_samples <- rownames(meta_data)
  
  if (!all(count_samples %in% meta_samples)) {
    missing_in_meta <- count_samples[!count_samples %in% meta_samples]
    stop(paste("Count data columns not found in metadata:", paste(missing_in_meta, collapse=", ")))
  }
  
  # Reorder metadata to match count data column order
  meta_data <- meta_data[count_samples, , drop = FALSE]
  
  # PCA 플롯의 점 색상을 결정할 그룹 변수 (config에서 읽어옴)
  intgroup <- config$de_analysis$group_variable
  if (!intgroup %in% colnames(meta_data)) {
    stop(paste("PCA intgroup '", intgroup, "' not found in metadata columns."))
  }
  
  # Config에서 고급 옵션 로드 (기본값 설정)
  adv_opts <- config$de_analysis$advanced_options
  blind_mode <- if(!is.null(adv_opts$vst_blind)) adv_opts$vst_blind else TRUE
  pca_ntop <- if(!is.null(adv_opts$pca_ntop)) adv_opts$pca_ntop else 500

  dge_method <- config$de_analysis$method
  
  if (dge_method == "DESeq2") {
    suppressPackageStartupMessages({
      library(DESeq2)
      source(here("src", "utils", "load_data.R"))
    })
    dds <- create_de_object(config_path = opt$config)
    vst_data <- vst(dds, blind = blind_mode)
    pca_data <- plotPCA(vst_data, intgroup = intgroup, returnData = TRUE, ntop = pca_ntop)
    percentVar <- round(100 * attr(pca_data, "percentVar"))

  } else if (dge_method %in% c("edgeR", "limma-voom")) {
    suppressPackageStartupMessages(library(edgeR))
    design_matrix <- model.matrix(as.formula(config$de_analysis$design_formula), data=meta_data)
    dge <- DGEList(counts = count_data)
    keep <- filterByExpr(dge, design_matrix)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)
    logcpm <- cpm(dge, log=TRUE)
    
    pca_res <- prcomp(t(logcpm))
    percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)))
    
    # ggplot을 위해 데이터프레임 생성. intgroup 변수를 동적으로 추가
    pca_data <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2])
    pca_data[[intgroup]] <- meta_data[[intgroup]]
  }
  
  # 공통 ggplot 시각화
  pca_plot <- ggplot(pca_data, aes_string(x = "PC1", y = "PC2", color = intgroup)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle("Global PCA Plot")
  
  ggsave(opt$output_file, plot = pca_plot)
  cat(paste("Global PCA plot saved to:", opt$output_file, "\n"))

} else if (opt$task == "volcano") {
  # --- Task 2b: Generate Pairwise Volcano Plot ---
  cat(paste("\nGenerating Volcano plot for:", opt$input_file, "\n"))
  
  if (is.null(opt$input_file)) {
    stop("Input file (--input_file) is required for volcano task.", call. = FALSE)
  }
  
  res <- read.csv(opt$input_file, row.names = 1)
  vp_aes <- config$plot_aesthetics$volcano
  
  # 데이터 가공
  res$diffexpressed <- "NO"
  res$diffexpressed[res$log2FoldChange > config$de_analysis$log2fc_cutoff & res$padj < config$de_analysis$padj_cutoff] <- "UP"
  res$diffexpressed[res$log2FoldChange < -config$de_analysis$log2fc_cutoff & res$padj < config$de_analysis$padj_cutoff] <- "DOWN"
  
  # 상위 N개 유전자 라벨링 로직
  res$delabel <- NA
  if (vp_aes$label_top_n > 0) {
    significant_genes <- subset(res, diffexpressed != "NO" & !is.na(symbol))
    significant_genes <- significant_genes[order(significant_genes$padj), ]
    top_n_genes <- head(significant_genes, vp_aes$label_top_n)
    res[rownames(top_n_genes), "delabel"] <- top_n_genes$symbol
  }
  
  # ggplot 객체 생성
  volcano_plot <- ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
    geom_point(size = vp_aes$point_size) +
    theme_minimal(base_size = vp_aes$base_font_size) +
    scale_color_manual(values = c("UP" = vp_aes$up_color, "DOWN" = vp_aes$down_color, "NO" = vp_aes$base_color)) +
    geom_text_repel(max.overlaps = Inf, na.rm = TRUE, size = 3.5, box.padding = 0.5) +
    ggtitle(vp_aes$title) +
    geom_vline(xintercept = c(-config$de_analysis$log2fc_cutoff, config$de_analysis$log2fc_cutoff), col = "red", linetype = 'dashed') +
    geom_hline(yintercept = -log10(config$de_analysis$padj_cutoff), col = "red", linetype = 'dashed') +
    theme(legend.position = "top")

  ggsave(opt$output_file, plot = volcano_plot, width = 10, height = 8, dpi = 300)
  cat(paste("Volcano plot saved to:", opt$output_file, "\n"))

} else {
  stop(paste("Invalid task:", opt$task, ". Choose 'pca' or 'volcano'."), call. = FALSE)
}