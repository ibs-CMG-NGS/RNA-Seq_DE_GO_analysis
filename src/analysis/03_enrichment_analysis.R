# 파일 경로: src/analysis/03_enrichment_analysis.R
# 사용법: Rscript 03_enrichment_analysis.R --config config.yml --input_csv [path/to/de_results.csv] --output_dir [path/to/output_pair_folder] --task [go|kegg] --geneset [up|down|total] --ontology [BP|CC|MF]

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
              help = "Path to the output directory for this pair", metavar = "character"),
  make_option(c("-t", "--task"), type = "character", 
              help = "Task to run: 'go' or 'kegg'", metavar = "character"),
  make_option(c("-g", "--geneset"), type = "character", 
              help = "Gene set to use: 'up', 'down', or 'total'", metavar = "character"),
  make_option(c("-n", "--ontology"), type = "character", default = NULL,
              help = "Ontology for GO: 'BP', 'CC', or 'MF' (required for 'go' task)", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_csv) || is.null(opt$output_dir) || is.null(opt$task) || is.null(opt$geneset)) {
  print_help(opt_parser)
  stop("input_csv, output_dir, task, and geneset arguments must be supplied.", call. = FALSE)
}

if (opt$task == "go" && is.null(opt$ontology)) {
  print_help(opt_parser)
  stop("ontology argument must be supplied for GO task.", call. = FALSE)
}

# Load config
if (!file.exists(opt$config)) {
  stop(paste("Config file not found at:", opt$config))
}
config <- yaml.load_file(opt$config)

# Debug: Print config info
cat(paste("Using config file:", opt$config, "\n"))
cat(paste("Species:", config$species, "\n"))

# Define output path based on arguments
output_path <- opt$output_dir
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# --- 1b. Load remaining libraries and define variables ---
if (!"databases" %in% names(config) || !config$species %in% names(config$databases)) {
    stop("[FATAL] 'databases' section or species entry missing in config.")
}
species_info <- config$databases[[config$species]]

# Organism DB
if (!"organism_db" %in% names(species_info)){
    stop("[FATAL] 'organism_db' key missing under species entry in config.")
}
organism_db_name <- species_info$organism_db
cat(paste("Loading organism database:", organism_db_name, "\n"))
if (!require(organism_db_name, character.only = TRUE)) {
    stop(paste("[FATAL] Required organism DB package", organism_db_name, "is not installed."))
}
organism_db <- get(organism_db_name)

# KEGG Organism Code
if (!"kegg_code" %in% names(species_info)){
    stop("[FATAL] 'kegg_code' key missing under species entry in config.")
}
kegg_organism <- species_info$kegg_code
if (is.null(kegg_organism) || kegg_organism == ""){
    stop("[FATAL] kegg_organism code is empty or NULL.")
}

# Dotplot aesthetics
if (!"plot_aesthetics" %in% names(config) || !"dotplot" %in% names(config$plot_aesthetics)) {
    stop("[FATAL] 'plot_aesthetics' section or 'dotplot' subsection missing.")
}
dp_aes <- config$plot_aesthetics$dotplot


# --- 2. DE 분석 결과 로드 ---
res_path <- opt$input_csv
if (!file.exists(res_path)) {
    stop(paste("[FATAL] DE results file not found at:", res_path))
}

# Read CSV without setting row names first
res <- read.csv(res_path, stringsAsFactors = FALSE, check.names = FALSE)

# The first column should contain gene IDs (Ensembl IDs)
# Set it as rownames explicitly
if (ncol(res) > 0) {
  gene_id_col <- res[, 1]
  rownames(res) <- as.character(gene_id_col)  # Force to character
  res <- res[, -1]  # Remove the first column after setting rownames
  cat(paste("Loaded", nrow(res), "genes from DE results\n"))
  cat(paste("Example gene IDs:", paste(head(rownames(res), 3), collapse=", "), "\n"))
}


# --- 3. 유전자 목록 준비 ---
gene_set <- opt$geneset
cat(paste("\n--- Preparing gene list for:", gene_set, "regulated genes ---\n"))

# Debug: Check if rownames are preserved
cat(paste("First few rownames in res:", paste(head(rownames(res), 3), collapse=", "), "\n"))

significant_genes <- res[!is.na(res$padj) & res$padj < config$de_analysis$padj_cutoff, ]

gene_list <- if (gene_set == "up") {
  significant_genes[significant_genes$log2FoldChange > config$de_analysis$log2fc_cutoff, ]
} else if (gene_set == "down") {
  significant_genes[significant_genes$log2FoldChange < -config$de_analysis$log2fc_cutoff, ]
} else { # "total"
  significant_genes
}

if (nrow(gene_list) == 0) { 
  cat("No significant genes. Skipping.\n")
  quit(status=0) 
}

# Get gene IDs from rownames (these are already Entrez IDs in this pipeline)
gene_ids <- rownames(gene_list)

cat(paste("Total genes in list:", length(gene_ids), "\n"))
cat(paste("First few gene IDs:", paste(head(gene_ids, 3), collapse=", "), "\n"))

# The gene IDs are already Entrez IDs, so use them directly
entrez_ids <- as.character(gene_ids)

# Remove any NA values
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

cat(paste("Using", length(entrez_ids), "Entrez IDs for enrichment analysis\n"))

if (length(entrez_ids) == 0) { 
  cat("No valid Entrez IDs. Skipping.\n")
  quit(status=0) 
}

# --- 4. Task 실행 ---

if (opt$task == "go") {
  # --- GO Enrichment Analysis ---
  ont <- opt$ontology
  cat(paste("Running GO analysis for", gene_set, "genes, Ontology:", ont, "\n"))
  
  # Force garbage collection to free memory before enrichGO
  gc()
  
  # Wrap enrichGO in tryCatch to handle potential errors gracefully
  go_results <- tryCatch({
    # Use simpler universe setting to avoid memory issues
    enrichGO(gene = entrez_ids, 
             OrgDb = organism_db, 
             keyType = 'ENTREZID', 
             ont = ont,
             pAdjustMethod = "BH", 
             pvalueCutoff = config$enrichment$pvalue_cutoff,
             qvalueCutoff = config$enrichment$qvalue_cutoff,
             readable = FALSE,  # Don't convert IDs to symbols (can cause issues)
             pool = FALSE)      # Don't pool gene sets (more stable)
  }, error = function(e) {
    cat(paste("Error in enrichGO:", e$message, "\n"))
    cat("Returning NULL result.\n")
    return(NULL)
  })
  
  # CSV 저장
  out_csv <- paste0("go_enrichment_", gene_set, "_", ont, ".csv")
  write.csv(as.data.frame(go_results), file.path(output_path, out_csv))
  
  out_plot <- paste0("go_dotplot_", gene_set, "_", ont, ".png")
  
  if (!is.null(go_results) && nrow(go_results) > 0) {
    plot_df <- as.data.frame(go_results) %>%
      mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
      arrange(p.adjust) %>%
      head(dp_aes$show_n_categories) %>%
      mutate(Description = fct_reorder(Description, .data[[dp_aes$x_axis_variable]]))
    
    go_dotplot <- ggplot(plot_df, aes_string(x = dp_aes$x_axis_variable, y = "Description", 
                                             color = "-log10(p.adjust)", size = "Count")) +
      geom_point() +
      scale_color_gradient(low = dp_aes$high_color, high = dp_aes$low_color) +
      labs(
        title = paste("GO Enrichment -", ont, "(", gene_set, "regulated)"),
        x = dp_aes$x_axis_variable, y = "GO Term", color = "-log10(p.adjust)", size = "Gene Count"
      ) +
      theme_minimal(base_size = dp_aes$font_size)

    ggsave(file.path(output_path, out_plot), plot = go_dotplot, width = 10, height = 8)
  } else {
    # No enrichment results - create an empty/placeholder plot
    cat("No GO enrichment results. Creating placeholder plot.\n")
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = paste("No significant GO enrichment found\nfor", gene_set, "regulated genes -", ont), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(output_path, out_plot), plot = empty_plot, width = 10, height = 8)
  }

} else if (opt$task == "kegg") {
  # --- KEGG Pathway Analysis ---
  cat(paste("Running KEGG analysis for", gene_set, "genes\n"))
  kegg_results <- enrichKEGG(gene = entrez_ids, organism = kegg_organism, pvalueCutoff = config$enrichment$pvalue_cutoff)
  
  out_csv_kegg <- paste0("kegg_enrichment_", gene_set, ".csv")
  write.csv(as.data.frame(kegg_results), file.path(output_path, out_csv_kegg))
  
  out_plot_kegg <- paste0("kegg_dotplot_", gene_set, ".png")
  
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    plot_df_kegg <- as.data.frame(kegg_results) %>%
      mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
      arrange(p.adjust) %>%
      head(dp_aes$show_n_categories) %>%
      mutate(Description = fct_reorder(Description, .data[[dp_aes$x_axis_variable]]))
      
    kegg_dotplot <- ggplot(plot_df_kegg, aes_string(x = dp_aes$x_axis_variable, y = "Description", 
                                                 color = "-log10(p.adjust)", size = "Count")) +
      geom_point() +
      scale_color_gradient(low = dp_aes$high_color, high = dp_aes$low_color) +
      labs(
        title = paste("KEGG Pathways (", gene_set, "regulated)"),
        x = dp_aes$x_axis_variable, y = "KEGG Pathway", color = "-log10(p.adjust)", size = "Gene Count"
      ) +
      theme_minimal(base_size = dp_aes$font_size)
    
    ggsave(file.path(output_path, out_plot_kegg), plot = kegg_dotplot, width = 10, height = 8)
  } else {
    # No enrichment results - create an empty/placeholder plot
    cat("No KEGG enrichment results. Creating placeholder plot.\n")
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = paste("No significant KEGG enrichment found\nfor", gene_set, "regulated genes"), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(output_path, out_plot_kegg), plot = empty_plot, width = 10, height = 8)
  }
} else {
  stop(paste("Invalid task:", opt$task))
}

cat("\nEnrichment analysis step finished successfully! 🚀\n")