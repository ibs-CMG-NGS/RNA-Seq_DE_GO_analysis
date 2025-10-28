# 파일 경로: src/analysis/04_generate_go_plots.R
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
config <- yaml.load_file(config_path) # Now 'config' is available

# Define output path based on config
output_path <- here(config$output_dir)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Load remaining required libraries for this script
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats) # Make sure forcats is loaded for fct_reorder
})

cat("\n--- Running Step 4: Generating GO Bar Plots ---\n")

# --- 1. 설정값 로드 ---
bar_aes <- config$go_barplot
gene_sets_to_plot <- config$enrichment$gene_lists

# --- 2. 유전자 그룹별로 반복 실행 ---
for (gene_set in gene_sets_to_plot) {
  
  cat(paste("\nProcessing gene set:", gene_set, "\n"))
  
  # --- 3. 모든 GO 결과 파일(BP, CC, MF)을 하나로 합치기 ---
  all_go_results <- list()
  for (ont in c("BP", "CC", "MF")) {
    input_csv <- file.path(output_path, paste0("go_enrichment_", gene_set, "_", ont, ".csv"))
    if (file.exists(input_csv)) {
      go_df <- read.csv(input_csv)
      if (nrow(go_df) > 0) { # 결과가 있는 파일만 추가
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
  
  # --- 4. [수정] 각 Ontology 별로 상위 N개 데이터 선택 ---
  
  plot_data <- combined_go %>%
    filter(ONTOLOGY %in% bar_aes$namespaces) %>%
    # ONTOLOGY(BP, CC, MF)로 그룹을 나눕니다.
    group_by(ONTOLOGY) %>%
    # 각 그룹 내에서 p.adjust 기준으로 정렬합니다.
    arrange(p.adjust) %>%
    # 각 그룹의 상위 N개(top_n)를 선택합니다.
    slice_head(n = bar_aes$top_n) %>%
    ungroup()

  if (nrow(plot_data) == 0) {
    cat(paste("No significant GO terms to plot for gene set:", gene_set, "after filtering. Skipping.\n"))
    next
  }

  # --- 5. [수정] ggplot Facet으로 Subplot 생성 ---
  
  # y축(Description)을 p.adjust 값에 따라 정렬하기 위해 factor 레벨 재정렬
  plot_data <- plot_data %>%
    mutate(Description = reorder(Description, -log10(p.adjust)))
  
  go_bar <- ggplot(plot_data, aes(x = -log10(p.adjust), y = Description, fill = ONTOLOGY)) +
    geom_col() +
    # [핵심] facet_wrap을 사용하여 ONTOLOGY별로 1x3 Subplot을 생성합니다.
    # scales = "free_y" 옵션은 각 subplot의 y축 눈금을 독립적으로 만듭니다.
    facet_wrap(~ ONTOLOGY, nrow = 3, scales = "free_y") + # ncol =3
    scale_fill_manual(values = bar_aes$colors) +
    labs(
      title = paste("Top", bar_aes$top_n, "GO Terms for", gene_set, "regulated genes"),
      x = "-log10(Adjusted P-value)",
      y = "GO Term"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 10),
      # facet 제목이 각 subplot 위에 표시되므로, 범례(legend)는 더 이상 필요 없습니다.
      legend.position = "none",
      # facet 제목 텍스트 스타일링
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # --- 6. 플롯 저장 ---
  output_plot_path <- file.path(output_path, paste0("go_barplot_", gene_set, ".png"))
  ggsave(output_plot_path, plot = go_bar, width = 10, height = 15, dpi = 300) # 가로 길이를 늘려 subplot이 잘 보이도록 조정
  
  cat(paste("Successfully generated and saved:", output_plot_path, "\n"))
}


cat("\n--- GO Bar Plot generation finished. ---\n")
