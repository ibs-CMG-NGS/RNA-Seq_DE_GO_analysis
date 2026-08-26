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

# 중간 산출물(geneset x ontology 조합으로 개수가 곱해지는 CSV/PNG)은 pair 루트가 아니라
# enrichment/(CSV)·plots/(PNG) 서브폴더에 저장해 루트를 final_*.xlsx 등 요약본 위주로 정리한다.
enrichment_dir <- file.path(output_path, "enrichment")
plots_dir      <- file.path(output_path, "plots")
dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,      showWarnings = FALSE, recursive = TRUE)

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
    write.csv(empty_df, file.path(enrichment_dir, out_csv), row.names = FALSE)
    
    # Placeholder plot
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0("No significant ", gene_set, " regulated genes found\n",
                            "(padj < ", config$de_analysis$padj_cutoff, 
                            ", |log2FC| > ", config$de_analysis$log2fc_cutoff, ")"), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(plots_dir, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
    
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
    write.csv(empty_df, file.path(enrichment_dir, out_csv), row.names = FALSE)
    
    # Placeholder plot
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0("No significant ", gene_set, " regulated genes found\n",
                            "(padj < ", config$de_analysis$padj_cutoff, 
                            ", |log2FC| > ", config$de_analysis$log2fc_cutoff, ")"), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(plots_dir, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
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
    cat(paste("Warning: mapIds failed with keytype '", keytype, "':", e$message, "\n"))
    cat("No Entrez ID mapping possible — empty enrichment results will be written.\n")
    setNames(rep(NA_character_, length(gene_ids)), gene_ids)
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

# CRITICAL: Remove duplicate Entrez IDs (keep unique only)
n_before_dedup <- length(entrez_ids)
entrez_ids <- unique(as.character(entrez_ids))
n_after_dedup <- length(entrez_ids)

if (n_before_dedup != n_after_dedup) {
  cat(sprintf("WARNING: Removed %d duplicate Entrez IDs (%d unique genes remain)\n", 
              n_before_dedup - n_after_dedup, n_after_dedup))
}

cat(paste("Final number of unique Entrez IDs for enrichment:", length(entrez_ids), "\n"))
cat(paste("Sample Entrez IDs:", paste(head(entrez_ids, 5), collapse=", "), "\n"))

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
    write.csv(empty_df, file.path(enrichment_dir, out_csv), row.names = FALSE)
    
    # Placeholder plot
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0("No valid Entrez IDs after conversion\n",
                            "Gene ID type: ", gene_id_type), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(plots_dir, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
    
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
    write.csv(empty_df, file.path(enrichment_dir, out_csv), row.names = FALSE)
    
    # Placeholder plot
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0("No valid Entrez IDs after conversion\n",
                            "Gene ID type: ", gene_id_type), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(plots_dir, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
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
  
  # Limit the number of genes to avoid memory issues and crashes
  # enrichGO can crash with too many genes due to AnnotationDbi memory issues
  max_genes_for_go <- 5000
  if (length(entrez_ids) > max_genes_for_go) {
    cat(paste("WARNING: Too many genes (", length(entrez_ids), ") for GO enrichment.\n"))
    cat(paste("Limiting to top", max_genes_for_go, "genes by adjusted p-value to prevent memory issues.\n"))
    
    # Get the order of genes by adjusted p-value (from original DE results)
    # Match entrez_ids back to original gene_list to get p-values
    original_genes <- rownames(gene_list)
    
    # Create a lookup for p-values
    padj_lookup <- setNames(gene_list$padj, original_genes)
    
    # For converted IDs, we need to map back
    if (gene_id_type != "ENTREZID") {
      # Map Entrez IDs back to original IDs to get p-values
      # This is approximate, but sufficient for limiting genes
      entrez_to_padj <- gene_list$padj[match(names(entrez_ids), original_genes)]
      entrez_to_padj[is.na(entrez_to_padj)] <- 1  # Set NA p-values to 1
      
      # Sort by p-value and take top genes
      top_indices <- order(entrez_to_padj)[1:min(max_genes_for_go, length(entrez_ids))]
      entrez_ids <- entrez_ids[top_indices]
    } else {
      # For ENTREZID, can match directly
      entrez_to_padj <- gene_list$padj[match(entrez_ids, original_genes)]
      entrez_to_padj[is.na(entrez_to_padj)] <- 1
      
      top_indices <- order(entrez_to_padj)[1:min(max_genes_for_go, length(entrez_ids))]
      entrez_ids <- entrez_ids[top_indices]
    }
    
    cat(paste("Using", length(entrez_ids), "genes for GO enrichment.\n"))
  }
  
  cat(paste("Running enrichGO with", length(entrez_ids), "unique genes...\n"))
  
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
  write.csv(as.data.frame(go_results), file.path(enrichment_dir, out_csv))
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

      ggsave(file.path(plots_dir, out_plot), plot = go_dotplot, width = 10, height = 8, bg = "white")
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
      ggsave(file.path(plots_dir, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
    }
  } else {
    # No results at all
    cat("No GO enrichment results. Creating placeholder plot.\n")
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
              label = paste("No GO enrichment found\nfor", gene_set, "regulated genes -", ont),
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(plots_dir, out_plot), plot = empty_plot, width = 10, height = 8, bg = "white")
  }

  # GO term clustering(term_cluster)과 GO Slim rollup(go_slim) 둘 다에서 쓰는 공용 헬퍼/값
  parse_ratio <- function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    parts[1] / parts[2]
  }
  compute_fold_enrichment <- function(go_df) {
    mapply(function(gr, br) parse_ratio(gr) / parse_ratio(br), go_df$GeneRatio, go_df$BgRatio)
  }

  # --- GO term clustering (up/down 세트에만 적용, "total"은 제외) ---
  # FDR<fdr_cutoff & FoldEnrichment>fold_enrichment_cutoff로 유의 term만 추린 뒤,
  # term간 유전자 중복도(Jaccard)로 pairwise_termsim() + treeplot() 클러스터링.
  tc_cfg <- config$enrichment$term_cluster
  tc_enabled <- if (is.null(tc_cfg) || is.null(tc_cfg$enabled)) TRUE else isTRUE(tc_cfg$enabled)

  if (tc_enabled && gene_set %in% c("up", "down") && !is.null(go_results) && nrow(go_results) > 0) {
    fdr_cutoff  <- ifelse(is.null(tc_cfg$fdr_cutoff), 0.05, tc_cfg$fdr_cutoff)
    fe_cutoff   <- ifelse(is.null(tc_cfg$fold_enrichment_cutoff), 2.0, tc_cfg$fold_enrichment_cutoff)
    # similarity_cutoff가 기본 모드(CMG-SeqViewer 자체 클러스터링 기본값과 동일: Jaccard 0.7).
    # 하위호환용으로 similarity_cutoff를 명시적으로 null로 두면 n_clusters(고정 k) 모드로 대체.
    similarity_cutoff <- if ("similarity_cutoff" %in% names(tc_cfg)) tc_cfg$similarity_cutoff else 0.7
    n_clusters  <- ifelse(is.null(tc_cfg$n_clusters), 5, tc_cfg$n_clusters)
    hclust_method <- ifelse(is.null(tc_cfg$hclust_method), "average", tc_cfg$hclust_method)
    # show_top_n은 PNG treeplot 시각화 전용 상한이다(treeplot()이 term이 많아질수록 렌더링이
    # 급격히 느려지고 200~500개 구간에서 내부 에러도 나므로). CSV/CMG용 엑셀 내보내기는 이
    # 제한과 무관하게 유의 term 전체(모든 클러스터+싱글톤)를 대상으로 한다.
    show_top_n  <- ifelse(is.null(tc_cfg$show_top_n), 30, tc_cfg$show_top_n)

    go_df_tc <- as.data.frame(go_results)
    go_df_tc$FoldEnrichment <- compute_fold_enrichment(go_df_tc)
    sig_ids <- go_df_tc$ID[!is.na(go_df_tc$p.adjust) & go_df_tc$p.adjust < fdr_cutoff &
                            !is.na(go_df_tc$FoldEnrichment) & go_df_tc$FoldEnrichment > fe_cutoff]

    out_termcluster     <- file.path(plots_dir, paste0("go_termcluster_", gene_set, "_", ont, ".png"))
    out_termcluster_csv <- file.path(enrichment_dir, paste0("go_termcluster_", gene_set, "_", ont, ".csv"))
    min_terms_needed <- 3

    if (length(sig_ids) >= min_terms_needed) {
      go_sig <- go_results
      go_sig@result <- go_df_tc[go_df_tc$ID %in% sig_ids, ]

      # --- 유의 term 전체를 클러스터링(show_top_n 제한 없음) — CSV/CMG 엑셀 내보내기용 ---
      # pairwise_termsim() 자체는 term이 수백 개여도 빠르다(느려지는/깨지는 건 treeplot()
      # 렌더링 쪽이라 아래에서 따로 표시 개수만 제한한다).
      clus_result <- tryCatch({
        # pairwise_termsim()은 showCategory 기본값이 200이라, 명시하지 않으면 유의 term이
        # 200개를 넘을 때 조용히 상위 200개로만 유사도 행렬을 잘라버린다(실제 확인된 실패
        # 사례: "subscript out of bounds"). 유의 term 전체를 클러스터링해야 하므로 명시한다.
        go_sig_termsim <- pairwise_termsim(go_sig, method = "JC", showCategory = length(sig_ids))

        keep <- seq_len(length(sig_ids))
        termsim2 <- go_sig_termsim@termsim[keep, keep]
        termsim2[is.na(termsim2)] <- 0
        termsim2 <- termsim2 + t(termsim2)
        diag(termsim2) <- 1
        hc_manual <- stats::hclust(stats::as.dist(1 - termsim2), method = hclust_method)

        if (!is.null(similarity_cutoff)) {
          clus <- stats::cutree(hc_manual, h = 1 - similarity_cutoff)
        } else {
          clus <- stats::cutree(hc_manual, k = min(n_clusters, length(sig_ids)))
        }
        list(termsim = go_sig_termsim, clus = clus)
      }, error = function(e) {
        cat(paste("[term_cluster] Failed:", e$message, "\n"))
        NULL
      })

      if (!is.null(clus_result)) {
        clus <- clus_result$clus
        effective_n <- length(unique(clus))
        n_singletons <- sum(table(clus) == 1)
        if (!is.null(similarity_cutoff)) {
          cat(sprintf("[term_cluster] similarity_cutoff=%.2f -> %d terms, %d clusters, %d singletons\n",
                       similarity_cutoff, length(sig_ids), effective_n, n_singletons))
        } else {
          cat(sprintf("[term_cluster] n_clusters=%d -> %d terms, %d clusters, %d singletons\n",
                       n_clusters, length(sig_ids), effective_n, n_singletons))
        }

        # --- CSV 저장: 유의 term 전체(모든 클러스터 + 싱글톤)의 클러스터 배정 ---
        # 주의: pairwise_termsim() 결과의 @termsim 행/열 이름은 GO ID가 아니라 Description이다.
        cluster_df <- go_df_tc[match(names(clus), go_df_tc$Description),
                                c("ID", "Description", "GeneRatio", "BgRatio", "FoldEnrichment",
                                  "pvalue", "p.adjust", "qvalue", "Count", "geneID")]
        cluster_df$cluster <- clus[cluster_df$Description]
        cluster_df <- cluster_df[order(cluster_df$cluster, cluster_df$p.adjust), ]
        write.csv(cluster_df, out_termcluster_csv, row.names = FALSE)
        cat(sprintf("Saved GO term cluster assignments (%d terms, %d clusters) to %s\n",
                     nrow(cluster_df), effective_n, basename(out_termcluster_csv)))

        # --- PNG 시각화 (best-effort, show_top_n으로 상위 term만 대표 표시 — 실패해도 위 CSV는 이미 저장됨) ---
        show_n <- min(as.integer(show_top_n), length(sig_ids))
        tp <- tryCatch({
          treeplot(clus_result$termsim, showCategory = show_n,
                   cluster.params = list(method = hclust_method, n = min(effective_n, show_n)))
        }, error = function(e) {
          cat(paste("[term_cluster] treeplot failed:", e$message, "\n"))
          NULL
        })

        if (!is.null(tp)) {
          # term/클러스터 수가 많으면 라벨이 겹치지 않도록 세로 길이를 동적으로 늘림
          plot_height <- max(8, show_n * 0.22, min(effective_n, show_n) * 0.4)
          ggsave(out_termcluster, plot = tp, width = 12, height = plot_height, bg = "white", limitsize = FALSE)
          cat(sprintf("Saved GO term cluster plot (%d of %d terms shown) to %s\n",
                       show_n, length(sig_ids), basename(out_termcluster)))
        }
      }
    } else {
      cat(sprintf("[term_cluster] Only %d terms pass FDR<%.3f & FoldEnrichment>%.1f (need >= %d) — skipping clustering, creating placeholder.\n",
                   length(sig_ids), fdr_cutoff, fe_cutoff, min_terms_needed))
      empty_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = sprintf("Not enough significant GO terms to cluster\n(FDR<%.3f & FoldEnrichment>%.1f)\nfor %s regulated genes - %s\n(%d term(s) found, need >= %d)",
                                 fdr_cutoff, fe_cutoff, gene_set, ont, length(sig_ids), min_terms_needed),
                 size = 5, hjust = 0.5) +
        theme_void()
      ggsave(out_termcluster, plot = empty_plot, width = 12, height = 8, bg = "white")
    }
  }

  # --- GO Slim rollup (up/down 세트에만 적용, "total"은 제외) ---
  # FDR<fdr_cutoff & FoldEnrichment>fold_enrichment_cutoff로 유의 term만 추린 뒤,
  # clusterProfiler::gofilter()로 GO DAG level 기준 상위 범주만 남긴다. 이 CSV들은
  # 05c_generate_go_slim_overview.R이 pair 단위로 up/down을 합쳐 대칭 bar chart를 그리는 데 쓰인다.
  slim_cfg <- config$enrichment$go_slim
  slim_enabled <- if (is.null(slim_cfg) || is.null(slim_cfg$enabled)) TRUE else isTRUE(slim_cfg$enabled)

  if (slim_enabled && gene_set %in% c("up", "down") && !is.null(go_results) && nrow(go_results) > 0) {
    slim_level <- ifelse(is.null(slim_cfg$level), 3, slim_cfg$level)
    slim_fdr   <- ifelse(is.null(slim_cfg$fdr_cutoff), 0.05, slim_cfg$fdr_cutoff)
    slim_fe    <- ifelse(is.null(slim_cfg$fold_enrichment_cutoff), 2.0, slim_cfg$fold_enrichment_cutoff)

    go_df_slim <- as.data.frame(go_results)
    go_df_slim$FoldEnrichment <- compute_fold_enrichment(go_df_slim)
    slim_sig_ids <- go_df_slim$ID[!is.na(go_df_slim$p.adjust) & go_df_slim$p.adjust < slim_fdr &
                                   !is.na(go_df_slim$FoldEnrichment) & go_df_slim$FoldEnrichment > slim_fe]

    out_slim_csv <- file.path(enrichment_dir, paste0("go_slim_", gene_set, "_", ont, ".csv"))
    if (length(slim_sig_ids) > 0) {
      go_sig_slim <- go_results
      go_sig_slim@result <- go_df_slim[go_df_slim$ID %in% slim_sig_ids, ]
      slim_result <- tryCatch(gofilter(go_sig_slim, level = slim_level), error = function(e) {
        cat(paste("[go_slim] gofilter failed:", e$message, "\n"))
        NULL
      })
      if (!is.null(slim_result) && nrow(slim_result) > 0) {
        write.csv(as.data.frame(slim_result), out_slim_csv, row.names = FALSE)
        cat(sprintf("[go_slim] level=%d -> %d/%d terms retained, saved to %s\n",
                     slim_level, nrow(slim_result), length(slim_sig_ids), basename(out_slim_csv)))
      } else {
        cat(sprintf("[go_slim] level=%d -> 0 terms retained (해당 level에 남는 term 없음) — CSV 생략.\n", slim_level))
      }
    } else {
      cat("[go_slim] No significant terms to roll up — CSV 생략.\n")
    }
  }

  # --- rrvgo 의미론적(semantic) 축약 (up/down 세트에만 적용, "total"은 제외) ---
  # FDR<fdr_cutoff & FoldEnrichment>fold_enrichment_cutoff로 유의 term 전체를 대상으로
  # GO DAG 의미 거리 기반 축약(parentTerm)을 계산한다. Jaccard 클러스터링/CMG 엑셀
  # 내보내기(term_cluster)의 입력(유의 term 전체 목록)은 그대로 유지되고, rrvgo 결과는
  # 별도 CSV/시각화로만 추가된다(term 목록을 사전에 줄이는 필터로 쓰지 않음).
  rrvgo_cfg <- config$enrichment$rrvgo
  rrvgo_enabled <- if (is.null(rrvgo_cfg) || is.null(rrvgo_cfg$enabled)) TRUE else isTRUE(rrvgo_cfg$enabled)

  if (rrvgo_enabled && gene_set %in% c("up", "down") && !is.null(go_results) && nrow(go_results) > 0) {
    suppressPackageStartupMessages(library(rrvgo))

    rrvgo_fdr       <- ifelse(is.null(rrvgo_cfg$fdr_cutoff), 0.05, rrvgo_cfg$fdr_cutoff)
    rrvgo_fe        <- ifelse(is.null(rrvgo_cfg$fold_enrichment_cutoff), 2.0, rrvgo_cfg$fold_enrichment_cutoff)
    rrvgo_method    <- ifelse(is.null(rrvgo_cfg$method), "Rel", rrvgo_cfg$method)
    rrvgo_threshold <- ifelse(is.null(rrvgo_cfg$threshold), 0.7, rrvgo_cfg$threshold)
    rrvgo_score_by  <- ifelse(is.null(rrvgo_cfg$score_by), "fdr", rrvgo_cfg$score_by)  # "fdr" | "count"

    go_df_rr <- as.data.frame(go_results)
    go_df_rr$FoldEnrichment <- compute_fold_enrichment(go_df_rr)
    rr_sig_ids <- go_df_rr$ID[!is.na(go_df_rr$p.adjust) & go_df_rr$p.adjust < rrvgo_fdr &
                               !is.na(go_df_rr$FoldEnrichment) & go_df_rr$FoldEnrichment > rrvgo_fe]

    out_rrvgo_csv     <- file.path(enrichment_dir, paste0("go_rrvgo_", gene_set, "_", ont, ".csv"))
    out_rrvgo_treemap <- file.path(plots_dir, paste0("go_rrvgo_treemap_", gene_set, "_", ont, ".png"))
    out_rrvgo_scatter <- file.path(plots_dir, paste0("go_rrvgo_scatter_", gene_set, "_", ont, ".png"))

    if (length(rr_sig_ids) >= 2) {
      go_df_rr_sig <- go_df_rr[go_df_rr$ID %in% rr_sig_ids, ]
      rr_scores <- if (rrvgo_score_by == "count") {
        setNames(go_df_rr_sig$Count, go_df_rr_sig$ID)
      } else {
        setNames(-log10(go_df_rr_sig$p.adjust), go_df_rr_sig$ID)
      }

      rr_result <- tryCatch({
        simMatrix <- calculateSimMatrix(rr_sig_ids, orgdb = organism_db_name, ont = ont, method = rrvgo_method)
        reducedTerms <- reduceSimMatrix(simMatrix, scores = rr_scores, threshold = rrvgo_threshold, orgdb = organism_db_name)
        list(simMatrix = simMatrix, reducedTerms = reducedTerms)
      }, error = function(e) {
        cat(paste("[rrvgo] Failed:", e$message, "\n"))
        NULL
      })

      if (!is.null(rr_result)) {
        reducedTerms <- rr_result$reducedTerms
        write.csv(reducedTerms, out_rrvgo_csv, row.names = FALSE)
        cat(sprintf("[rrvgo] %d terms -> %d parent groups (threshold=%.2f), saved to %s\n",
                     length(rr_sig_ids), length(unique(reducedTerms$parent)), rrvgo_threshold, basename(out_rrvgo_csv)))

        sp <- tryCatch(scatterPlot(rr_result$simMatrix, reducedTerms), error = function(e) {
          cat(paste("[rrvgo] scatterPlot failed:", e$message, "\n"))
          NULL
        })
        if (!is.null(sp)) {
          ggsave(out_rrvgo_scatter, plot = sp, width = 10, height = 8, bg = "white")
          cat(sprintf("[rrvgo] Saved scatter plot to %s\n", basename(out_rrvgo_scatter)))
        }

        tm_ok <- tryCatch({
          png(out_rrvgo_treemap, width = 12, height = 8, units = "in", res = 300, bg = "white")
          treemapPlot(reducedTerms)
          dev.off()
          TRUE
        }, error = function(e) {
          cat(paste("[rrvgo] treemapPlot failed:", e$message, "\n"))
          if (dev.cur() > 1) dev.off()
          FALSE
        })
        if (isTRUE(tm_ok)) cat(sprintf("[rrvgo] Saved treemap to %s\n", basename(out_rrvgo_treemap)))
      }
    } else {
      cat("[rrvgo] Fewer than 2 significant terms — skipping semantic reduction.\n")
    }
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
  write.csv(as.data.frame(kegg_results), file.path(enrichment_dir, out_csv_kegg))
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
      
      ggsave(file.path(plots_dir, out_plot_kegg), plot = kegg_dotplot, width = 10, height = 8, bg = "white")
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
      ggsave(file.path(plots_dir, out_plot_kegg), plot = empty_plot, width = 10, height = 8, bg = "white")
    }
  } else {
    # No results at all
    cat("No KEGG enrichment results. Creating placeholder plot.\n")
    empty_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
              label = paste("No KEGG enrichment found\nfor", gene_set, "regulated genes"), 
               size = 6, hjust = 0.5) +
      theme_void()
    ggsave(file.path(plots_dir, out_plot_kegg), plot = empty_plot, width = 10, height = 8, bg = "white")
  }
} else {
  stop(paste("Invalid task:", opt$task))
}

cat("\nEnrichment analysis step finished successfully! 🚀\n")