# 파일 경로: src/utils/load_data.R

# 필요한 라이브러리 로드
library(DESeq2)
library(yaml)
library(here)

create_de_object <- function(config_path = here("config.yml")) {
  
  # 설정 파일 로드
  if (!exists("config")) {
      config <- yaml.load_file(config_path)
  }

  # here()를 사용해 데이터 경로를 프로젝트 최상위 폴더 기준으로 지정
  count_data_path <- here(config$count_data_path)
  metadata_path <- here(config$metadata_path)

  # 데이터 로드
  count_data <- read.csv(count_data_path, row.names = 1, check.names = FALSE)
  meta_data <- read.csv(metadata_path, row.names = 1)

  # Ensure sample names match between count data and metadata
  # Count data columns should match metadata row names
  count_samples <- colnames(count_data)
  meta_samples <- rownames(meta_data)
  
  # Check if samples match
  if (!all(count_samples %in% meta_samples)) {
    missing_in_meta <- count_samples[!count_samples %in% meta_samples]
    stop(paste("Count data columns not found in metadata:", paste(missing_in_meta, collapse=", ")))
  }
  
  # Reorder metadata to match count data column order
  meta_data <- meta_data[count_samples, , drop = FALSE]
  
  cat(paste("Loaded", nrow(count_data), "genes and", ncol(count_data), "samples\n"))

  # config 파일에서 'de_analysis' 섹션의 design_formula를 읽어옵니다.
  if (!"de_analysis" %in% names(config) || !"design_formula" %in% names(config$de_analysis)) {
      stop("ERROR: 'de_analysis' section or 'design_formula' is missing in config.yml.")
  }
  design_formula <- as.formula(config$de_analysis$design_formula)

  # Pre-filtering 적용
  prefilter_thres <- config$de_analysis$advanced_options$prefilter_threshold
  if (!is.null(prefilter_thres) && prefilter_thres > 0) {
      keep_genes <- rowSums(count_data) >= prefilter_thres
      count_data <- count_data[keep_genes, ]
  }

  # DESeqDataSet 객체 생성
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = meta_data,
                                design = design_formula)

  return(dds)
}