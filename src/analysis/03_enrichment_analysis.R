# 파일 경로: src/analysis/03_enrichment_analysis.R
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

# Define output path based on config
output_path <- here(config$output_dir)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Load remaining required libraries for this script
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(AnnotationDbi) # Needed for mapIds
  # [수정] Ensure the correct organism DB package is loaded AND assigned
  species_info <- config$databases[[config$species]]
  organism_db_name <- species_info$organism_db # e.g., "org.Hs.eg.db"

  # Load the library using the string name
  if (!require(organism_db_name, character.only = TRUE)) {
      stop(paste("Required organism DB package", organism_db_name, "is not installed."))
  }
  # [핵심!] Convert the string name into the actual R object and assign it
  organism_db <- get(organism_db_name)

})

# --- 2. DE 분석 결과 로드 ---
res_path <- file.path(output_path, "final_de_results.csv")
res <- read.csv(res_path, row.names = 1)

# --- 3. 유전자 목록 및 Ontology 조합에 따른 반복 분석 ---
cat("Starting enrichment analysis based on config settings...\n")

for (gene_set in config$enrichment$gene_lists) {
  
  cat(paste("\n--- Preparing gene list for:", gene_set, "regulated genes ---\n"))
  
  # ... (유전자 목록 선택 및 ID 변환 부분은 이전과 동일) ...
  significant_genes <- subset(res, !is.na(padj) & padj < config$de_analysis$padj_cutoff)
  gene_list <- if (gene_set == "up") {
    subset(significant_genes, log2FoldChange > config$de_analysis$log2fc_cutoff)
  } else if (gene_set == "down") {
    subset(significant_genes, log2FoldChange < -config$de_analysis$log2fc_cutoff)
  } else { # "total"
    significant_genes
  }
  if (nrow(gene_list) == 0) { cat("No significant genes. Skipping.\n"); next }
  entrez_ids <- mapIds(organism_db, keys = rownames(gene_list), column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  if (length(entrez_ids) == 0) { cat("No Entrez IDs mapped. Skipping.\n"); next }

  # --- GO Enrichment Analysis (내부 반복문) ---
  for (ont in config$enrichment$go_ontologies) {
    cat(paste("Running GO analysis for", gene_set, "genes, Ontology:", ont, "\n"))
    
    go_results <- enrichGO(gene = entrez_ids, OrgDb = organism_db, keyType = 'ENTREZID', ont = ont,
                           pAdjustMethod = "BH", pvalueCutoff = config$enrichment$pvalue_cutoff,
                           qvalueCutoff = config$enrichment$qvalue_cutoff)
    
    if (!is.null(go_results) && nrow(go_results) > 0) {
      
      plot_df <- as.data.frame(go_results)
      
      # Define x_var and check if it exists in plot_df
      x_var <- dp_aes$x_axis_variable
      if (! x_var %in% colnames(plot_df)) {
        message(sprintf("Warning: x_axis_variable '%s' not found in plot_df columns. Available columns: %s", x_var, paste(colnames(plot_df), collapse = ", ")))
      }
      
      # [수정] dplyr 파이프라인으로 데이터 가공을 하나로 합칩니다.
      plot_df <- as.data.frame(go_results) %>%
        # 1. GeneRatio를 숫자로 변환
        mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
        # 2. p.adjust 기준으로 정렬
        arrange(p.adjust) %>%
        # 3. 상위 N개 선택
        head(dp_aes$show_n_categories) %>%
        # 4. y축 정렬
        mutate(Description = fct_reorder(Description, .data[[dp_aes$x_axis_variable]]))
      
      # Ensure numeric and handle p.adjust == 0
      eps <- 1e-300
      plot_df <- plot_df %>%
        mutate(p.adjust = as.numeric(p.adjust),
               Count = as.numeric(Count),
               p.adjust = ifelse(is.na(p.adjust), NA_real_, ifelse(p.adjust == 0, eps, p.adjust)),
               log10padj = -log10(p.adjust))
      
      bad_count <- plot_df %>% 
        filter(is.na(.data[[x_var]]) | is.na(log10padj) | is.na(Count) | !is.finite(.data[[x_var]]) | !is.finite(log10padj)) %>% 
        nrow()
      message("Rows that will be removed by ggplot due to NA/Inf in mapped aesthetics: ", bad_count)
      
      go_dotplot <- ggplot(plot_df, aes(x = .data[[x_var]], y = Description, 
                                               color = log10padj, size = Count)) +
        geom_point() +
        scale_color_gradient(low = "blue", high = "red") +
        labs(
          title = paste("GO Enrichment -", ont, "(", gene_set, "regulated)"),
          x = x_var, y = "GO Term", color = "-log10(p.adjust)", size = "Gene Count"
        ) +
        theme_minimal(base_size = 14)

      out_plot <- paste0("go_dotplot_", gene_set, "_", ont, ".png")
      ggsave(file.path(output_path, out_plot), plot = go_dotplot, width = 10, height = 8)
    }
    # 결과 파일 저장은 그대로 유지
    out_csv <- paste0("go_enrichment_", gene_set, "_", ont, ".csv")
    write.csv(as.data.frame(go_results), file.path(output_path, out_csv))
  }

  # --- KEGG Pathway Analysis (동일한 로직 적용) ---
  cat(paste("Running KEGG analysis for", gene_set, "genes\n"))
  kegg_results <- enrichKEGG(gene = entrez_ids, organism = kegg_organism, pvalueCutoff = config$enrichment$pvalue_cutoff)
  
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    
    # Define x_var_kegg and check if it exists in plot_df_kegg
    x_var_kegg <- dp_aes$x_axis_variable
    plot_df_kegg <- as.data.frame(kegg_results)
    if (! x_var_kegg %in% colnames(plot_df_kegg)) {
      message(sprintf("Warning: x_axis_variable '%s' not found in plot_df_kegg columns. Available columns: %s", x_var_kegg, paste(colnames(plot_df_kegg), collapse = ", ")))
    }
    
    # [수정] KEGG 부분도 동일하게 dplyr 파이프라인으로 개선합니다.
    plot_df_kegg <- as.data.frame(kegg_results) %>%
      mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
      arrange(p.adjust) %>%
      head(dp_aes$show_n_categories) %>%
      mutate(Description = fct_reorder(Description, .data[[dp_aes$x_axis_variable]]))
    
    # Ensure numeric and handle p.adjust == 0
    eps <- 1e-300
    plot_df_kegg <- plot_df_kegg %>%
      mutate(p.adjust = as.numeric(p.adjust),
             Count = as.numeric(Count),
             p.adjust = ifelse(is.na(p.adjust), NA_real_, ifelse(p.adjust == 0, eps, p.adjust)),
             log10padj_kegg = -log10(p.adjust))
    
    bad_count_kegg <- plot_df_kegg %>% 
      filter(is.na(.data[[x_var_kegg]]) | is.na(log10padj_kegg) | is.na(Count) | !is.finite(.data[[x_var_kegg]]) | !is.finite(log10padj_kegg)) %>% 
      nrow()
    message("Rows that will be removed by ggplot due to NA/Inf in mapped aesthetics: ", bad_count_kegg)
    
    kegg_dotplot <- ggplot(plot_df_kegg, aes(x = .data[[x_var_kegg]], y = Description, 
                                                 color = log10padj_kegg, size = Count)) +
      geom_point() +
      scale_color_gradient(low = "blue", high = "red") +
      labs(
        title = paste("KEGG Pathways (", gene_set, "regulated)"),
        x = x_var_kegg, y = "KEGG Pathway", color = "-log10(p.adjust)", size = "Gene Count"
      ) +
      theme_minimal(base_size = 14)
    
    out_plot_kegg <- paste0("kegg_dotplot_", gene_set, ".png")
    ggsave(file.path(output_path, out_plot_kegg), plot = kegg_dotplot, width = 10, height = 8)
  }
  # 결과 파일 저장은 그대로 유지
  out_csv_kegg <- paste0("kegg_enrichment_", gene_set, ".csv")
  write.csv(as.data.frame(kegg_results), file.path(output_path, out_csv_kegg))
}


cat("\nEnrichment analysis pipeline finished successfully! 🚀\n")

