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

# Handle case with no significant genes
if (nrow(gene_list) == 0) { 
  cat("No significant genes. Creating empty output files.\n")
  
  # Create empty CSV file
  if (opt$task == "go") {
    out_csv <- paste0("go_enrichment_", gene_set, "_", opt$ontology, ".csv")
    out_plot <- paste0("go_dotplot_", gene_set, "_", opt$ontology, ".png")
    
    # Empty CSV with header
    empty_df <- data.frame(
      ID = character(),
      Description = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer()
    )
    write.csv(empty_df, file.path(output_path, out_csv), row.names = FALSE)
    
    # Placeholder plot
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0("No significant ", gene_set, " regulated genes found\n",
                            "(padj < ", config$de_analysis$padj_cutoff, 
                            ", |log2FC| > ", config$de_analysis$log2fc_cutoff, ")"), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(output_path, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
    
  } else if (opt$task == "kegg") {
    out_csv <- paste0("kegg_enrichment_", gene_set, ".csv")
    out_plot <- paste0("kegg_dotplot_", gene_set, ".png")
    
    # Empty CSV with header
    empty_df <- data.frame(
      ID = character(),
      Description = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer()
    )
    write.csv(empty_df, file.path(output_path, out_csv), row.names = FALSE)
    
    # Placeholder plot
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0("No significant ", gene_set, " regulated genes found\n",
                            "(padj < ", config$de_analysis$padj_cutoff, 
                            ", |log2FC| > ", config$de_analysis$log2fc_cutoff, ")"), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(output_path, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
  }
  
  cat("Empty output files created successfully.\n")
  quit(status=0) 
}

# Get gene IDs from rownames
gene_ids <- rownames(gene_list)

cat(paste("Total genes in list:", length(gene_ids), "\n"))
cat(paste("First few gene IDs:", paste(head(gene_ids, 3), collapse=", "), "\n"))

# --- Determine gene ID type and convert to Entrez if needed ---
# Step 1: Get gene_id_type from config (if specified)
gene_id_type <- if ("gene_id_type" %in% names(config)) {
  toupper(config$gene_id_type)
} else {
  # Auto-detect based on ID pattern
  sample_id <- head(gene_ids[!is.na(gene_ids)], 1)
  if (length(sample_id) == 0) {
    stop("No valid gene IDs found for type detection")
  }
  
  detected_type <- if (grepl("^ENSMUSG[0-9]+", sample_id)) {
    "ENSEMBL"  # Mouse Ensembl
  } else if (grepl("^ENSG[0-9]+", sample_id)) {
    "ENSEMBL"  # Human Ensembl
  } else if (grepl("^[0-9]+$", sample_id)) {
    "ENTREZID"  # Numeric = Entrez
  } else if (grepl("^[A-Z][A-Z0-9]+$", sample_id)) {
    "SYMBOL"  # Gene symbol
  } else {
    warning(paste("Could not auto-detect gene ID type. Sample ID:", sample_id))
    "UNKNOWN"
  }
  
  cat(paste("Auto-detected gene_id_type:", detected_type, "(based on sample ID:", sample_id, ")\n"))
  detected_type
}

cat(paste("Gene ID type:", gene_id_type, "\n"))

# Step 2: Convert to Entrez IDs if needed
if (gene_id_type == "ENTREZID") {
  # Already Entrez IDs - use directly
  cat("Gene IDs are already in Entrez format. No conversion needed.\n")
  entrez_ids <- as.character(gene_ids)
  
} else if (gene_id_type %in% c("ENSEMBL", "ENSEMBLID", "SYMBOL")) {
  # Need to convert to Entrez
  cat(paste("Converting", gene_id_type, "IDs to Entrez IDs...\n"))
  
  # Determine keytype for mapIds
  keytype <- if (gene_id_type == "SYMBOL") {
    "SYMBOL"
  } else {
    # Try both ENSEMBL and ENSEMBLID
    "ENSEMBL"
  }
  
  entrez_ids <- tryCatch({
    mapIds(organism_db, 
           keys = gene_ids,
           column = "ENTREZID",
           keytype = keytype,
           multiVals = "first")
  }, error = function(e) {
    if (gene_id_type %in% c("ENSEMBL", "ENSEMBLID")) {
      # Try alternative keytype
      alt_keytype <- if (keytype == "ENSEMBL") "ENSEMBLID" else "ENSEMBL"
      cat(paste("Retrying with keytype:", alt_keytype, "\n"))
      tryCatch({
        mapIds(organism_db, 
               keys = gene_ids,
               column = "ENTREZID",
               keytype = alt_keytype,
               multiVals = "first")
      }, error = function(e2) {
        stop(paste("Failed to convert gene IDs. Error:", e2$message))
      })
    } else {
      stop(paste("Failed to convert gene IDs. Error:", e$message))
    }
  })
  
  # Remove NA values
  n_before <- length(entrez_ids)
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  n_after <- length(entrez_ids)
  
  cat(paste("Successfully converted", n_after, "out of", n_before, "genes to Entrez IDs\n"))
  cat(paste("Conversion rate:", round(100 * n_after / n_before, 1), "%\n"))
  
} else {
  stop(paste("Unsupported gene_id_type:", gene_id_type, 
             "\nSupported types: ENSEMBL, ENTREZID, SYMBOL"))
}

# Final validation
if (length(entrez_ids) == 0) { 
  cat("No valid Entrez IDs after conversion. Creating empty output files.\n")
  
  # Create empty CSV file
  if (opt$task == "go") {
    out_csv <- paste0("go_enrichment_", gene_set, "_", opt$ontology, ".csv")
    out_plot <- paste0("go_dotplot_", gene_set, "_", opt$ontology, ".png")
    
    # Empty CSV with header
    empty_df <- data.frame(
      ID = character(),
      Description = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer()
    )
    write.csv(empty_df, file.path(output_path, out_csv), row.names = FALSE)
    
    # Placeholder plot
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0("No valid Entrez IDs after conversion\n",
                            "Gene ID type: ", gene_id_type), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(output_path, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
    
  } else if (opt$task == "kegg") {
    out_csv <- paste0("kegg_enrichment_", gene_set, ".csv")
    out_plot <- paste0("kegg_dotplot_", gene_set, ".png")
    
    # Empty CSV with header
    empty_df <- data.frame(
      ID = character(),
      Description = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer()
    )
    write.csv(empty_df, file.path(output_path, out_csv), row.names = FALSE)
    
    # Placeholder plot
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0("No valid Entrez IDs after conversion\n",
                            "Gene ID type: ", gene_id_type), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(output_path, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
  }
  
  cat("Empty output files created successfully.\n")
  quit(status=0) 
}

# --- 4. Task 실행 ---

if (opt$task == "go") {
  # --- GO Enrichment Analysis ---
  ont <- opt$ontology
  cat(paste("Running GO analysis for", gene_set, "genes, Ontology:", ont, "\n"))
  
  # Force garbage collection to free memory before enrichGO
  gc()
  
  # Limit the number of genes to avoid memory issues
  # If too many genes, take top genes by significance
  max_genes_for_go <- 2000
  if (length(entrez_ids) > max_genes_for_go) {
    cat(paste("Warning: Too many genes (", length(entrez_ids), "). Limiting to top", max_genes_for_go, "by p-value.\n"))
    # Sort by adjusted p-value and take top genes
    top_indices <- order(res_sig$padj)[1:min(max_genes_for_go, nrow(res_sig))]
    entrez_ids <- entrez_ids[top_indices]
    cat(paste("Using", length(entrez_ids), "genes for GO enrichment.\n"))
  }
  
  # Wrap enrichGO in tryCatch to handle potential errors gracefully
  # Note: Use pvalueCutoff=1 and qvalueCutoff=1 to get ALL results
  #       Filtering will be applied later for visualization
  go_results <- tryCatch({
    # Use simpler universe setting to avoid memory issues
    enrichGO(gene = entrez_ids, 
             OrgDb = organism_db, 
             keyType = 'ENTREZID', 
             ont = ont,
             pAdjustMethod = "BH", 
             pvalueCutoff = 1.0,  # Get all results
             qvalueCutoff = 1.0,  # Get all results
             readable = FALSE,  # Don't convert IDs to symbols (can cause issues)
             pool = FALSE,      # Don't pool gene sets (more stable)
             minGSSize = ifelse(is.null(config$enrichment$min_gs_size), 10, config$enrichment$min_gs_size),
             maxGSSize = ifelse(is.null(config$enrichment$max_gs_size), 500, config$enrichment$max_gs_size))
  }, error = function(e) {
    cat(paste("Error in enrichGO:", e$message, "\n"))
    cat("Returning NULL result.\n")
    return(NULL)
  })
  
  # Apply gene count filter if specified
  min_gene_count <- ifelse(is.null(config$enrichment$min_gene_count), 1, config$enrichment$min_gene_count)
  if (!is.null(go_results) && nrow(go_results) > 0) {
    go_df <- as.data.frame(go_results)
    original_count <- nrow(go_df)
    go_df <- go_df %>% filter(Count >= min_gene_count)
    if (nrow(go_df) < original_count) {
      cat(sprintf("Filtered out %d terms with Count < %d\n", original_count - nrow(go_df), min_gene_count))
    }
    go_results@result <- go_df
  }
  
  # CSV 저장 (모든 결과 저장)
  out_csv <- paste0("go_enrichment_", gene_set, "_", ont, ".csv")
  write.csv(as.data.frame(go_results), file.path(output_path, out_csv))
  cat(sprintf("Saved %d GO terms to %s\n", nrow(as.data.frame(go_results)), out_csv))
  
  # Dotplot 생성 (필터링 적용)
  out_plot <- paste0("go_dotplot_", gene_set, "_", ont, ".png")
  plot_top_n <- ifelse(is.null(config$enrichment$plot_top_n), dp_aes$show_n_categories, config$enrichment$plot_top_n)
  
  if (!is.null(go_results) && nrow(go_results) > 0) {
    # Filter for visualization: apply p-value cutoff and top N
    plot_df <- as.data.frame(go_results) %>%
      filter(p.adjust < config$enrichment$pvalue_cutoff) %>%
      mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
      arrange(p.adjust) %>%
      head(plot_top_n) %>%
      mutate(Description = fct_reorder(Description, .data[[dp_aes$x_axis_variable]]))
    
    if (nrow(plot_df) > 0) {
      go_dotplot <- ggplot(plot_df, aes_string(x = dp_aes$x_axis_variable, y = "Description", 
                                               color = "-log10(p.adjust)", size = "Count")) +
        geom_point() +
        scale_color_gradient(low = dp_aes$high_color, high = dp_aes$low_color) +
        labs(
          title = paste("GO Enrichment -", ont, "(", gene_set, "regulated)"),
          x = dp_aes$x_axis_variable, y = "GO Term", color = "-log10(p.adjust)", size = "Gene Count"
        ) +
        theme_minimal(base_size = dp_aes$font_size)

      ggsave(file.path(output_path, out_plot), plot = go_dotplot, width = 10, height = 8, bg = "white")
      cat(sprintf("Saved dotplot with %d terms to %s\n", nrow(plot_df), out_plot))
    } else {
      cat(sprintf("No terms pass p.adjust < %.3f threshold for plotting\n", config$enrichment$pvalue_cutoff))
      # Create placeholder plot
      empty_plot <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = sprintf("No significant GO enrichment\n(p.adjust < %.3f)\nfor %s regulated genes - %s", 
                              config$enrichment$pvalue_cutoff, gene_set, ont), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(output_path, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
  }

} else if (opt$task == "kegg") {
  # --- KEGG Pathway Analysis ---
  cat(paste("Running KEGG analysis for", gene_set, "genes\n"))
  
  # Use pvalueCutoff=1 to get ALL results
  kegg_results <- enrichKEGG(gene = entrez_ids, organism = kegg_organism, 
                            pvalueCutoff = 1.0,  # Get all results
                            minGSSize = ifelse(is.null(config$enrichment$min_gs_size), 10, config$enrichment$min_gs_size),
                            maxGSSize = ifelse(is.null(config$enrichment$max_gs_size), 500, config$enrichment$max_gs_size))
  
  # Apply gene count filter if specified
  min_gene_count <- ifelse(is.null(config$enrichment$min_gene_count), 1, config$enrichment$min_gene_count)
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    kegg_df <- as.data.frame(kegg_results)
    original_count <- nrow(kegg_df)
    kegg_df <- kegg_df %>% filter(Count >= min_gene_count)
    if (nrow(kegg_df) < original_count) {
      cat(sprintf("Filtered out %d pathways with Count < %d\n", original_count - nrow(kegg_df), min_gene_count))
    }
    kegg_results@result <- kegg_df
  }
  
  # CSV 저장 (모든 결과 저장)
  out_csv_kegg <- paste0("kegg_enrichment_", gene_set, ".csv")
  write.csv(as.data.frame(kegg_results), file.path(output_path, out_csv_kegg))
  cat(sprintf("Saved %d KEGG pathways to %s\n", nrow(as.data.frame(kegg_results)), out_csv_kegg))
  
  # Dotplot 생성 (필터링 적용)
  out_plot_kegg <- paste0("kegg_dotplot_", gene_set, ".png")
  plot_top_n <- ifelse(is.null(config$enrichment$plot_top_n), dp_aes$show_n_categories, config$enrichment$plot_top_n)
  
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    # Filter for visualization
    plot_df_kegg <- as.data.frame(kegg_results) %>%
      filter(p.adjust < config$enrichment$pvalue_cutoff) %>%
      mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
      arrange(p.adjust) %>%
      head(plot_top_n) %>%
      mutate(Description = fct_reorder(Description, .data[[dp_aes$x_axis_variable]]))
    
    if (nrow(plot_df_kegg) > 0) {
      kegg_dotplot <- ggplot(plot_df_kegg, aes_string(x = dp_aes$x_axis_variable, y = "Description", 
                                                   color = "-log10(p.adjust)", size = "Count")) +
        geom_point() +
        scale_color_gradient(low = dp_aes$high_color, high = dp_aes$low_color) +
        labs(
          title = paste("KEGG Pathways (", gene_set, "regulated)"),
          x = dp_aes$x_axis_variable, y = "KEGG Pathway", color = "-log10(p.adjust)", size = "Gene Count"
        ) +
        theme_minimal(base_size = dp_aes$font_size)
      
      ggsave(file.path(output_path, out_plot_kegg), plot = kegg_dotplot, width = 10, height = 8, bg = "white")
      cat(sprintf("Saved KEGG dotplot with %d pathways to %s\n", nrow(plot_df_kegg), out_plot_kegg))
    } else {
      cat(sprintf("No pathways pass p.adjust < %.3f threshold for plotting\n", config$enrichment$pvalue_cutoff))
      # Create placeholder plot
      empty_plot <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = sprintf("No significant KEGG enrichment\n(p.adjust < %.3f)\nfor %s regulated genes", 
                              config$enrichment$pvalue_cutoff, gene_set), 
                 size = 6, hjust = 0.5) +
        theme_void()
      ggsave(file.path(output_path, out_plot_kegg), plot = empty_plot, width = 10, height = 8, bg = "white")
    }
  } else {
    # No results at all
    cat("No KEGG enrichment results. Creating placeholder plot.\n")
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
              label = paste("No KEGG enrichment found\nfor", gene_set, "regulated genes"), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(output_path, out_plot_kegg), plot = empty_plot, width = 10, height = 8, bg = "white")
  }
} else {
  stop(paste("Invalid task:", opt$task))
}

cat("\nEnrichment analysis step finished successfully! 🚀\n")