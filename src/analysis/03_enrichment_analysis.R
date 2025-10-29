# 파일 경로: src/analysis/03_enrichment_analysis.R
# --- 1. Setup: Load config and libraries ---

# Suppress startup messages
suppressPackageStartupMessages({
  library(here)
  library.dynam('yaml', 'yaml', '/home/ygkim/program/anaconda3/envs/rna-seq-de-go-analysis/lib/R/library') # Trying explicit load just in case
  library(yaml)
})

# Get config file path
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    cat("[DEBUG] No config path from args, using default.\n")
    config_path <- here("config.yml")
} else {
    cat(paste("[DEBUG] Config path from args:", args[1], "\n"))
    config_path <- args[1]
}

# Check if config file exists
if (!file.exists(config_path)) {
  stop(paste("[FATAL] Config file not found at:", config_path))
}
cat(paste("[DEBUG] Config file found at:", config_path, "\n"))

# Load the config file
config <- NULL # Initialize to NULL
tryCatch({
    config <- yaml.load_file(config_path)
    cat("[DEBUG] Config file loaded successfully.\n")
}, error = function(e){
    stop(paste("[FATAL] Failed to load or parse config YAML:", e$message))
})

# Check if config object is valid
if (is.null(config)) {
    stop("[FATAL] Config object is NULL after loading.")
}
cat("[DEBUG] Config object is not NULL.\n")

# Define output path
output_path <- here(config$output_dir)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
cat(paste("[DEBUG] Output path set to:", output_path, "\n"))

# Load remaining libraries AND DEFINE dp_aes
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(AnnotationDbi)

  # --- Check and define organism_db AND kegg_organism ---
  if (!"databases" %in% names(config) || !config$species %in% names(config$databases)) {
      stop("[FATAL] 'databases' section or species entry missing in config.")
  }
  species_info <- config$databases[[config$species]] # <-- species_info 정의

  # Organism DB
  if (!"organism_db" %in% names(species_info)){ # <-- species_info 사용
      stop("[FATAL] 'organism_db' key missing under species entry in config.")
  }
  organism_db_name <- species_info$organism_db
  cat(paste("[DEBUG] Organism DB name:", organism_db_name, "\n"))
  if (!require(organism_db_name, character.only = TRUE)) {
      stop(paste("[FATAL] Required organism DB package", organism_db_name, "is not installed."))
  }
  organism_db <- get(organism_db_name)
  cat("[DEBUG] organism_db object created successfully.\n")

  # [수정] KEGG Organism Code - species_info 변수와 "kegg_code" 키 이름 확인
  if (!"kegg_code" %in% names(species_info)){ # <-- species_info 사용, "kegg_code" 확인
      stop("[FATAL] 'kegg_code' key missing under species entry in config.")
  }
  kegg_organism <- species_info$kegg_code # <-- species_info 사용
  cat(paste("[DEBUG] KEGG organism code:", kegg_organism, "\n"))
  if (is.null(kegg_organism) || kegg_organism == ""){
      stop("[FATAL] kegg_organism code is empty or NULL.")
  }

  # --- Check and define dp_aes ---
  if (!"plot_aesthetics" %in% names(config)) {
      stop("[FATAL] 'plot_aesthetics' section missing in config.")
  }
  cat("[DEBUG] 'plot_aesthetics' section found in config.\n")

  if (!"dotplot" %in% names(config$plot_aesthetics)) {
      stop("[FATAL] 'dotplot' subsection missing under 'plot_aesthetics' in config.")
  }
  cat("[DEBUG] 'dotplot' subsection found under 'plot_aesthetics'.\n")

  # Define dp_aes
  dp_aes <- config$plot_aesthetics$dotplot

  # Final check if dp_aes was assigned
  if (is.null(dp_aes)) {
       stop("[FATAL] dp_aes is NULL after assignment!")
  }
  cat("[DEBUG] dp_aes object created successfully. Contains:\n")
  print(dp_aes) # Print the content of dp_aes

})

# --- [수정] 2. DE 분석 결과 로드 ---
# 이 부분이 Setup 블록 다음, for 반복문 이전에 와야 합니다!
cat("[INFO] Loading DE analysis results...\n")
res_path <- file.path(output_path, "final_de_results.csv")
if (!file.exists(res_path)) {
    stop(paste("[FATAL] DE results file not found at:", res_path))
}
res <- read.csv(res_path, row.names = 1)
cat("[INFO] DE analysis results loaded successfully.\n")

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
    
    # Check if results are not empty before proceeding to plotting
      if (!is.null(go_results) && nrow(go_results) > 0) {

        # --- [NEW DEBUGGING STEP] ---
        cat("[DEBUG] Inside GO results block. Checking for dp_aes...\n")
        if (!exists("dp_aes")) {
            stop("[FATAL] dp_aes object does NOT exist right before plotting!")
        } else {
            cat("[DEBUG] dp_aes object FOUND right before plotting. Content:\n")
            print(dp_aes)
        }
      
      # [수정] dplyr 파이프라인으로 데이터 가공을 하나로 합칩니다.
      plot_df <- as.data.frame(go_results) %>%
        mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
        arrange(p.adjust) %>%
        head(dp_aes$show_n_categories) %>% # <--- First use
        mutate(Description = fct_reorder(Description, .data[[dp_aes$x_axis_variable]])) # <--- Second use

      x_var <- dp_aes$x_axis_variable # <--- Third use
      
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

      # --- [NEW DEBUGGING STEP] ---
      cat("[DEBUG] Inside KEGG results block. Checking for dp_aes...\n")
      if (!exists("dp_aes")) {
           stop("[FATAL] dp_aes object does NOT exist right before plotting!")
      } else {
           cat("[DEBUG] dp_aes object FOUND right before plotting. Content:\n")
           print(dp_aes)
      }
    
    # [수정] KEGG 부분도 동일하게 dplyr 파이프라인으로 개선합니다.
    plot_df_kegg <- as.data.frame(kegg_results) %>%
        mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
        arrange(p.adjust) %>%
        head(dp_aes$show_n_categories) %>% # <--- Use
        mutate(Description = fct_reorder(Description, .data[[dp_aes$x_axis_variable]])) # <--- Use

      x_var_kegg <- dp_aes$x_axis_variable # <--- Use
    
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
