#!/usr/bin/env Rscript
# 파일 경로: run_pipeline.R
# RNA-Seq 분석 파이프라인 CLI 진입점

# 필요한 라이브러리 로드
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(here)
})

# 명령줄 인자 정의
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yml",
              help = "Path to config file [default: %default]", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory path [default: timestamp-based directory in output/]", metavar = "DIR"),
  make_option(c("-s", "--step"), type = "character", default = "all",
              help = "Analysis step to run: 'all', 'de', 'plots', 'enrichment', 'go_plots' [default: %default]",
              metavar = "STEP"),
  make_option(c("--install-deps"), action = "store_true", default = FALSE,
              help = "Install required R packages before running analysis"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print verbose output")
)

opt_parser <- OptionParser(
  usage = "Usage: %prog [options]",
  option_list = option_list,
  description = "\nRNA-Seq Differential Expression and GO/KEGG Analysis Pipeline",
  epilogue = "Examples:
  Rscript run_pipeline.R                                # Run full pipeline with default config
  Rscript run_pipeline.R --config my_config.yml         # Use custom config file
  Rscript run_pipeline.R --step de                      # Run only DE analysis
  Rscript run_pipeline.R --output my_output_dir         # Specify output directory
  Rscript run_pipeline.R --install-deps                 # Install dependencies first
"
)

# 인자 파싱
opt <- parse_args(opt_parser)

# 함수: 메시지 출력 (verbose 모드일 때만)
print_verbose <- function(msg) {
  if (opt$verbose) {
    cat(paste0("[", Sys.time(), "] ", msg, "\n"))
  }
}

# 함수: 메시지 출력 (항상)
print_info <- function(msg) {
  cat(paste0("[INFO] ", msg, "\n"))
}

# 함수: 에러 메시지 출력 및 종료
print_error <- function(msg) {
  cat(paste0("[ERROR] ", msg, "\n"), file = stderr())
  quit(status = 1)
}

# 시작 메시지
cat("\n")
cat("============================================================\n")
cat("  RNA-Seq Differential Expression and GO/KEGG Analysis\n")
cat("============================================================\n\n")

# 1. 의존성 설치 (옵션)
if (opt$`install-deps`) {
  print_info("Installing required R packages...")
  tryCatch({
    source(here("src", "utils", "install_dependencies.R"))
    print_info("All dependencies installed successfully.")
  }, error = function(e) {
    print_error(paste("Failed to install dependencies:", e$message))
  })
}

# 2. 설정 파일 로드
config_path <- here(opt$config)
if (!file.exists(config_path)) {
  print_error(paste("Config file not found:", config_path))
}

print_verbose(paste("Loading config from:", config_path))
config <- tryCatch({
  yaml.load_file(config_path)
}, error = function(e) {
  print_error(paste("Failed to load config file:", e$message))
})

# 3. 출력 디렉토리 설정
if (is.null(opt$output)) {
  # 타임스탬프 기반 디렉토리 생성
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  output_path <- here(config$output_dir, timestamp)
} else {
  output_path <- here(opt$output)
}

# 출력 디렉토리 생성
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
  print_info(paste("Created output directory:", output_path))
} else {
  print_info(paste("Using existing output directory:", output_path))
}

# 전역 변수로 설정 (분석 스크립트에서 사용)
assign("config", config, envir = .GlobalEnv)
assign("output_path", output_path, envir = .GlobalEnv)

# 4. 분석 단계 실행
step <- tolower(opt$step)

run_de_analysis <- function() {
  print_info("Step 1: Running differential expression analysis...")
  tryCatch({
    source(here("src", "analysis", "01_run_de_analysis.R"))
    print_info("DE analysis completed successfully.")
  }, error = function(e) {
    print_error(paste("DE analysis failed:", e$message))
  })
}

run_plots <- function() {
  print_info("Step 2: Generating visualization plots...")
  tryCatch({
    source(here("src", "analysis", "02_generate_plots.R"))
    print_info("Plot generation completed successfully.")
  }, error = function(e) {
    print_error(paste("Plot generation failed:", e$message))
  })
}

run_enrichment <- function() {
  print_info("Step 3: Running enrichment analysis...")
  tryCatch({
    source(here("src", "analysis", "03_enrichment_analysis.R"))
    print_info("Enrichment analysis completed successfully.")
  }, error = function(e) {
    print_error(paste("Enrichment analysis failed:", e$message))
  })
}

run_go_plots <- function() {
  print_info("Step 4: Generating GO plots...")
  tryCatch({
    source(here("src", "analysis", "04_generate_go_plots.R"))
    print_info("GO plot generation completed successfully.")
  }, error = function(e) {
    print_error(paste("GO plot generation failed:", e$message))
  })
}

# 분석 단계 실행
if (step == "all") {
  run_de_analysis()
  run_plots()
  run_enrichment()
  run_go_plots()
} else if (step == "de") {
  run_de_analysis()
} else if (step == "plots") {
  run_plots()
} else if (step == "enrichment") {
  run_enrichment()
} else if (step == "go_plots") {
  run_go_plots()
} else {
  print_error(paste("Invalid step:", step, "\nValid options: all, de, plots, enrichment, go_plots"))
}

# 완료 메시지
cat("\n")
cat("============================================================\n")
cat("  Analysis completed successfully!\n")
cat(paste0("  Results saved to: ", output_path, "\n"))
cat("============================================================\n\n")
