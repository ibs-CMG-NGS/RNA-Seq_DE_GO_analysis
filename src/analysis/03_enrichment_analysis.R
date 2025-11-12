# 파일 경로: src/analysis/03_enrichment_analysis.R
# 사용법: Rscript 03_enrichment_analysis.R --config config.yml --input_csv [path/to/de_results.csv] --output_dir [path/to/output_pair_folder]

# --- 1. Setup: Load config and libraries ---
suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(AnnotationDbi)
})

# Argument parsing
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yml", 
              help = "Path to the config YAML file", metavar = "character"),
  make_option(c("-i", "--input_csv"), type = "character", 
              help = "Path to input DE results file (e.g., final_de_results.csv)", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", 
              help = "Path to the output directory for this pair", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_csv) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("input_csv and output_dir arguments must be supplied.", call. = FALSE)
}

# Load config
config <- yaml.load_file(opt$config)

# [수정] output_path 변수를 Snakemake 인자로부터 설정
output_path <- opt$output_dir
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# --- 1b. Load remaining libraries and define variables ---
species_info <- config$databases[[config$species]]
organism_db_name <- species_info$organism_db
kegg_organism <- species_info$kegg_code
if (!require(organism_db_name, character.only = TRUE)) {
  stop(paste("Genome package", organism_db_name, "is not installed."))
}
organism_db <- get(organism_db_name)
dp_aes <- config$plot_aesthetics$dotplot

# --- 2. DE 분석 결과 로드 ---
# [수정] res_path를 Snakemake 인자로부터 설정
res_path <- opt$input_csv
res <- read.csv(res_path, row.names = 1)

# --- 3. 유전자 목록 및 Ontology 조합에 따른 반복 분석 ---
cat("Starting enrichment analysis based on config settings...\n")

for (gene_set in config$enrichment$gene_lists) {
  
  cat(paste("\n--- Preparing gene list for:", gene_set, "regulated genes ---\n"))
  
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
    
    # CSV 저장 (output_path는 스크립트 상단에서 정의됨)
    out_csv <- paste0("go_enrichment_", gene_set, "_", ont, ".csv")
    write.csv(as.data.frame(go_results), file.path(output_path, out_csv))
    
    if (!is.null(go_results) && nrow(go_results) > 0) {
      
      plot_df <- as.data.frame(go_results) %>%
        mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
        arrange(p.adjust) %>%
        head(dp_aes$show_n_categories) %>%
        mutate(Description = fct_reorder(Description, .data[[dp_aes$x_axis_variable]]))
      
      x_var <- dp_aes$x_axis_variable
      
      go_dotplot <- ggplot(plot_df, aes_string(x = x_var, y = "Description", 
                                               color = "-log10(p.adjust)", size = "Count")) +
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
  }

  # --- KEGG Pathway Analysis ---
  cat(paste("Running KEGG analysis for", gene_set, "genes\n"))
  kegg_results <- enrichKEGG(gene = entrez_ids, organism = kegg_organism, pvalueCutoff = config$enrichment$pvalue_cutoff)
  
  out_csv_kegg <- paste0("kegg_enrichment_", gene_set, ".csv")
  write.csv(as.data.frame(kegg_results), file.path(output_path, out_csv_kegg))
  
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {

    plot_df_kegg <- as.data.frame(kegg_results) %>%
      mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
      arrange(p.adjust) %>%
      head(dp_aes$show_n_categories)
      
    x_var_kegg <- dp_aes$x_axis_variable
    plot_df_kegg <- plot_df_kegg %>%
      mutate(Description = fct_reorder(Description, .data[[x_var_kegg]]))
    
    kegg_dotplot <- ggplot(plot_df_kegg, aes_string(x = x_var_kegg, y = "Description", 
                                                 color = "-log10(p.adjust)", size = "Count")) +
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
} # end of for (gene_set ...)

cat("\nEnrichment analysis pipeline finished successfully! 🚀\n")