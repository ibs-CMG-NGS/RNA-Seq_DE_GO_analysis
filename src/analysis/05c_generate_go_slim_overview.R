# File: src/analysis/05c_generate_go_slim_overview.R
#
# Purpose: 03_enrichment_analysis.R의 GO Slim rollup 단계(go_slim_{geneset}_{ontology}.csv)를
# ontology별로 up/down 합쳐서 "30초 방향 파악용" 대칭 bar chart로 만든다.
#
# Usage: Rscript 05c_generate_go_slim_overview.R [config_path] [compare_group] [base_group] [output_dir]

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 05c_generate_go_slim_overview.R [config_path] [compare_group] [base_group] [output_dir]")
}
config_path   <- args[1]
compare_group <- args[2]
base_group    <- args[3]
output_dir    <- args[4]

config <- yaml.load_file(config_path)

slim_cfg <- config$enrichment$go_slim
slim_enabled <- if (is.null(slim_cfg) || is.null(slim_cfg$enabled)) TRUE else isTRUE(slim_cfg$enabled)

if (!slim_enabled) {
  cat("[05c_generate_go_slim_overview] go_slim.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

slim_level  <- ifelse(is.null(slim_cfg$level), 3, slim_cfg$level)
ontologies  <- config$enrichment$go_ontologies %||% c("BP", "CC", "MF")

cat(paste("[05c_generate_go_slim_overview]", compare_group, "vs", base_group, "\n"))

save_placeholder <- function(ont, msg) {
  empty_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg, size = 5, hjust = 0.5) +
    theme_void()
  ggsave(file.path(output_dir, paste0("go_slim_overview_", ont, ".png")), plot = empty_plot,
         width = 10, height = 6, bg = "white")
}

# BP는 항상 파일이 생성됨을 보장한다(Snakemake output으로 이 파일 하나만 추적하므로).
# CC/MF는 결과가 있을 때만 생성되는 보너스 산출물.
for (ont in ontologies) {
  up_csv   <- file.path(output_dir, paste0("go_slim_up_", ont, ".csv"))
  down_csv <- file.path(output_dir, paste0("go_slim_down_", ont, ".csv"))

  up_df   <- if (file.exists(up_csv))   read.csv(up_csv,   stringsAsFactors = FALSE) else NULL
  down_df <- if (file.exists(down_csv)) read.csv(down_csv, stringsAsFactors = FALSE) else NULL

  if ((is.null(up_df) || nrow(up_df) == 0) && (is.null(down_df) || nrow(down_df) == 0)) {
    cat(paste("  [", ont, "] go_slim 결과 없음 — 건너뜀\n"))
    if (ont == "BP") save_placeholder(ont, sprintf("No significant GO Slim (level %d) terms for BP.", slim_level))
    next
  }

  combined <- bind_rows(
    if (!is.null(up_df)   && nrow(up_df)   > 0) mutate(up_df,   direction = "Up")   else NULL,
    if (!is.null(down_df) && nrow(down_df) > 0) mutate(down_df, direction = "Down") else NULL
  )

  out_csv <- file.path(output_dir, paste0("go_slim_overview_", ont, ".csv"))
  write.csv(combined, out_csv, row.names = FALSE)

  plot_df <- combined %>%
    mutate(signed_count = ifelse(direction == "Up", Count, -Count)) %>%
    arrange(desc(abs(signed_count)))

  p <- ggplot(plot_df, aes(x = reorder(Description, abs(signed_count)),
                            y = signed_count, fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#E74C3C", "Down" = "#3498DB")) +
    scale_y_continuous(labels = abs) +
    labs(title = sprintf("GO Slim Overview - %s (level %d)", ont, slim_level),
         subtitle = paste(compare_group, "vs", base_group),
         x = NULL, y = "Gene Count", fill = "Direction") +
    theme_minimal(base_size = 12)

  out_png <- file.path(output_dir, paste0("go_slim_overview_", ont, ".png"))
  plot_height <- max(6, nrow(plot_df) * 0.3)
  ggsave(out_png, plot = p, width = 10, height = plot_height, bg = "white")

  cat(sprintf("  [%s] %d Slim 범주(Up %d / Down %d) 저장: %s, %s\n",
              ont, nrow(plot_df),
              sum(plot_df$direction == "Up"), sum(plot_df$direction == "Down"),
              basename(out_csv), basename(out_png)))
}

# BP가 go_ontologies 설정에서 아예 빠진 경우를 위한 안전장치 — Snakemake가 이 파일 하나만
# output으로 추적하므로 항상 존재를 보장해야 한다.
if (!"BP" %in% ontologies && !file.exists(file.path(output_dir, "go_slim_overview_BP.png"))) {
  save_placeholder("BP", "BP not included in enrichment.go_ontologies config.")
}

cat("[05c_generate_go_slim_overview] Done.\n")
