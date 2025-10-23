# 파일 경로: src/utils/install_dependencies.R

# --- 1. 핵심 패키지(yaml, here, optparse)를 먼저 확인하고 설치합니다 ---
# 이 스크립트 자체를 실행하는 데 필수적인 패키지들이므로, 가장 먼저 처리합니다.
core_packages <- c("yaml", "here", "optparse")
for (pkg in core_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing core package:", pkg, "\n"))
    install.packages(pkg)
  }
}

# --- 2. 핵심 라이브러리를 로드합니다 ---
# 이제 패키지가 반드시 설치되어 있으므로, 안전하게 라이브러리를 로드할 수 있습니다.
library(yaml)
library(here)


# --- 3. 나머지 패키지 목록을 정의합니다 ---

# CRAN에서 설치할 나머지 패키지 목록
cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "openxlsx",
  "dplyr",
  "forcats"
)

# config 파일을 읽어 분석할 종(species)에 맞는 DB 패키지 이름을 동적으로 결정
config <- yaml.load_file(here("config.yml"))
species_info <- config$databases[[config$species]]
organism_db_pkg <- species_info$organism_db

# Bioconductor에서 설치할 패키지 목록
bioc_packages <- c(
  "DESeq2",
  "edgeR",
  "limma",
  "AnnotationDbi",
  "clusterProfiler",
  "enrichplot",
  organism_db_pkg
)

# --- 4. 나머지 패키지들을 설치합니다 ---

# CRAN 패키지 설치
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing CRAN package:", pkg, "\n"))
    install.packages(pkg)
  }
}

# BiocManager가 없으면 먼저 설치
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor 패키지 설치
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing Bioconductor package:", pkg, "\n"))
    BiocManager::install(pkg, update = FALSE)
  }
}

cat("\nAll required packages are installed and ready. 🚀\n")