# 파일 경로: src/analysis/02_generate_plots.R

# --- 1. Setup: Load config and libraries ---

# Suppress startup messages
suppressPackageStartupMessages({
  library(here)
  library(yaml)
})

# Get config file path from command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    cat("No config file provided. Using default 'config.yml'\n")
    config_path <- here("config.yml")
} else {
    config_path <- args[1]
}

# Load the config file
if (!file.exists(config_path)) {
  stop(paste("Config file not found at:", config_path))
}
config <- yaml.load_file(config_path)

# Define output path based on config (needed for saving plots)
output_path <- here(config$output_dir)
# Ensure output directory exists (though Snakemake likely creates it)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Load remaining required libraries for this script
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  # Libraries needed only for DESeq2 PCA path
  if(config$de_analysis$method == "DESeq2") {
      library(DESeq2)
      source(here("src", "utils", "create_de_object.R")) # Source the updated function name
  }
  # Libraries needed only for edgeR/limma PCA path
  if(config$de_analysis$method %in% c("edgeR", "limma-voom")) {
      library(edgeR) # For DGEList, cpm
  }
})

# [수정] 1. config 설정에 따라 PCA 데이터 생성 방식을 분기합니다.
dge_method <- config$de_analysis$method
cat(paste("\nGenerating PCA plot based on data prepared for:", dge_method, "\n"))

# PCA 플롯의 점 색상을 결정할 그룹 변수를 design formula에서 자동으로 추출합니다.
# (예: "~ dex" -> "dex", "~ condition + batch" -> "condition")
intgroup <- all.vars(as.formula(config$de_analysis$design_formula))[1]

if (dge_method == "DESeq2") {
  library(DESeq2)
  # DESeq2 방식: vst 변환 후 plotPCA 함수 사용
  source(here("src", "utils", "load_data.R"))
  
  # [수정] 변경된 함수 이름을 호출합니다.
  dds <- create_de_object()

  vst_data <- vst(dds, blind = FALSE)
  
  # DESeq2의 plotPCA 함수를 사용해 PCA 데이터와 분산(%)을 바로 얻습니다.
  pca_data <- plotPCA(vst_data, intgroup = intgroup, returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

} else if (dge_method %in% c("edgeR", "limma-voom")) {
  library(edgeR)
  # edgeR/limma 방식: log-CPM 변환 후 직접 PCA 수행
  count_data <- read.csv(here(config$count_data_path), row.names = 1)
  meta_data <- read.csv(here(config$metadata_path), row.names = 1)
  design_matrix <- model.matrix(as.formula(config$de_analysis$design_formula), data=meta_data)
  
  dge <- DGEList(counts = count_data)
  keep <- filterByExpr(dge, design_matrix)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)

  # TMM 정규화된 log-CPM 값을 계산합니다.
  logcpm <- cpm(dge, log=TRUE)

  # 기본 R 함수 prcomp()를 사용해 PCA를 직접 계산합니다.
  pca_res <- prcomp(t(logcpm))
  
  # 각 주성분(PC)이 설명하는 분산(%)을 계산합니다.
  percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)))

  # ggplot으로 시각화하기 위해 데이터프레임을 만듭니다.
  pca_data <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], meta_data)

} else {
    stop("PCA plot generation is not defined for the selected DGE method.")
}

# [수정] 2. 공통 ggplot 시각화 코드
# aes_string()을 사용해 intgroup 변수에 담긴 그룹 이름을 동적으로 사용합니다.
pca_plot <- ggplot(pca_data, aes_string(x = "PC1", y = "PC2", color = intgroup)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle(paste("PCA Plot (Data transformed via", dge_method, "method)"))

ggsave(file.path(output_path, "pca_plot.png"), plot = pca_plot)
print("PCA plot saved.")


# --- Volcano Plot ---

# [수정된 부분] 새로운 결과 파일 이름을 사용합니다.
res_path <- file.path(output_path, "final_de_results.csv")
res <- read.csv(res_path, row.names = 1) # row.names=1을 유지하여 유전자 ID를 행 이름으로 사용
vp_aes <- config$plot_aesthetics$volcano

# config에서 plot 속성 가져오기
vp_aes <- config$plot_aesthetics$volcano

# 데이터 가공 (config 경로 수정)
res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > config$de_analysis$log2fc_cutoff & res$padj < config$de_analysis$padj_cutoff] <- "UP" # 수정
res$diffexpressed[res$log2FoldChange < -config$de_analysis$log2fc_cutoff & res$padj < config$de_analysis$padj_cutoff] <- "DOWN" # 수정

# 라벨링할 상위 유전자 데이터 준비
# 1. 라벨을 표시할 컬럼(delabel)을 우선 NA로 초기화합니다.
res$delabel <- NA

# 2. config 파일에 설정된 label_top_n 값이 0보다 클 때만 라벨링을 수행합니다.
if (vp_aes$label_top_n > 0) {
  
  # 3. 유의미한 유전자들만 따로 추출합니다.
  significant_genes <- subset(res, diffexpressed != "NO" & !is.na(symbol))
  
  # 4. padj(adjusted p-value) 기준으로 오름차순 정렬합니다.
  significant_genes <- significant_genes[order(significant_genes$padj), ]
  
  # 5. 정렬된 유전자 중 상위 N개를 선택합니다. (N이 실제 유전자 수보다 많아도 안전)
  top_n_genes <- head(significant_genes, vp_aes$label_top_n)
  
  # 6. 원본 데이터(res)에서 top_n_genes에 해당하는 행의 symbol을 delabel 컬럼에 할당합니다.
  res[rownames(top_n_genes), "delabel"] <- top_n_genes$symbol
}

# ggplot 객체 생성
volcano_plot <- ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_point(size = vp_aes$point_size) +
  theme_minimal(base_size = vp_aes$base_font_size) +
  scale_color_manual(values = c("UP" = vp_aes$up_color, "DOWN" = vp_aes$down_color, "NO" = vp_aes$base_color)) +
  # max.overlaps = Inf 옵션으로 라벨이 최대한 많이 보이도록 설정 (필요시 조절)
  geom_text_repel(max.overlaps = Inf, na.rm = TRUE, size = 3.5, box.padding = 0.5) +
  ggtitle(vp_aes$title) +
  geom_vline(xintercept = c(-config$de_analysis$log2fc_cutoff, config$de_analysis$log2fc_cutoff), col = "red", linetype = 'dashed') + # 수정
geom_hline(yintercept = -log10(config$de_analysis$padj_cutoff), col = "red", linetype = 'dashed') + # 수정
  theme(legend.position = "top")

ggsave(file.path(output_path, "volcano_plot.png"), plot = volcano_plot, width = 10, height = 8, dpi = 300)

print("Volcano plot saved.")
