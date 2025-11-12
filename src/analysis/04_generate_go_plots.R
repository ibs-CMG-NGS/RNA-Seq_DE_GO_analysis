# 파일 경로: src/analysis/04_generate_go_plots.R
# 사용법: Rscript 04_generate_go_plots.R --config config.yml --output_dir [path/to/output_pair_folder]

# --- 1. Setup: Load config and libraries ---
suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(forcats)
})

# Argument parsing
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yml", 
              help = "Path to the config YAML file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", 
              help = "Path to the output directory containing GO CSVs", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("output_dir argument must be supplied.", call. = FALSE)
}

# Load config
config <- yaml.load_file(opt$config)

# [수정] output_path 변수를 Snakemake 인자로부터 설정
output_path <- opt$output_dir

# --- 2. 설정값 로드 ---
cat("\n--- Running Step 4: Generating GO Bar Plots ---\n")
bar_aes <- config$go_barplot
gene_sets_to_plot <- config$enrichment$gene_lists

# --- 3. 유전자 그룹별로 반복 실행 ---
for (gene_set in gene_sets_to_plot) {
  
  cat(paste("\nProcessing gene set:", gene_set, "\n"))
  
  # --- 4. 모든 GO 결과 파일(BP, CC, MF)을 하나로 합치기 ---
  all_go_results <- list()
  for (ont in c("BP", "CC", "MF")) {
    # [수정] input_csv 경로는 인자로 받은 output_path를 기준으로 합니다.
    input_csv <- file.path(output_path, paste0("go_enrichment_", gene_set, "_", ont, ".csv"))
    
    if (file.exists(input_csv)) {
      go_df <- read.csv(input_csv)
      if (nrow(go_df) > 0) { 
        if (!"ONTOLOGY" %in% names(go_df)) {
          go_df$ONTOLOGY <- ont
        }
        all_go_results[[ont]] <- go_df
      }
    }
  }
  
  if (length(all_go_results) == 0) {
    cat(paste("No GO result files found for gene set:", gene_set, ". Skipping plot generation.\n"))
    next
  }
  combined_go <- bind_rows(all_go_results)
  
  # --- 5. 각 Ontology 별로 상위 N개 데이터 선택 ---
  plot_data <- combined_go %>%
    filter(ONTOLOGY %in% bar_aes$namespaces) %>%
    group_by(ONTOLOGY) %>%
    arrange(p.adjust) %>%
    slice_head(n = bar_aes$top_n) %>%
    ungroup()

  if (nrow(plot_data) == 0) {
    cat(paste("No significant GO terms to plot for gene set:", gene_set, "after filtering. Skipping.\n"))
    next
  }

  # --- 6. ggplot Facet으로 Subplot 생성 ---
  plot_data <- plot_data %>%
    mutate(Description = fct_reorder(Description, -log10(p.adjust)))
  
  go_bar <- ggplot(plot_data, aes(x = -log10(p.adjust), y = Description, fill = ONTOLOGY)) +
    geom_col() +
    facet_wrap(~ ONTOLOGY, nrow = 3, scales = "free_y") +
    scale_fill_manual(values = bar_aes$colors) +
    labs(
      title = paste("Top", bar_aes$top_n, "GO Terms for", gene_set, "regulated genes"),
      x = "-log10(Adjusted P-value)",
      y = "GO Term"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # --- 7. 플롯 저장 ---
  output_plot_path <- file.path(output_path, paste0("go_barplot_", gene_set, ".png"))
  ggsave(output_plot_path, plot = go_bar, width = 10, height = 15, dpi = 300) 
  
  cat(paste("Successfully generated and saved:", output_plot_path, "\n"))
}

cat("\n--- GO Bar Plot generation finished. ---\n")