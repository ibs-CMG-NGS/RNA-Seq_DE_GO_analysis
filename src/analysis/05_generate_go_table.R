# File: src/analysis/05_generate_go_table.R
# 
# Purpose: Generate GO enrichment results as publication-ready Excel file
#
# Usage: Rscript 05_generate_go_table.R [config_path] [compare_group] [base_group] [output_dir]

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
  library(openxlsx)
  library(stringr)
  library(AnnotationDbi)
})

# --- 1. Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 05_generate_go_table.R [config_path] [compare_group] [base_group] [output_dir]")
}
config_path <- args[1]
compare_group <- args[2]
base_group <- args[3]
output_dir <- args[4]

# --- 2. Load config ---
if (!file.exists(config_path)) {
  stop(paste("Config file not found at:", config_path))
}
config <- yaml.load_file(config_path)

cat("\n==============================================\n")
cat("  Generating GO Summary Table for Publication\n")
cat("==============================================\n")
cat(paste("Comparison:", compare_group, "vs", base_group, "\n"))
cat(paste("Output directory:", output_dir, "\n\n"))

# --- 2b. Load organism database for gene symbol conversion ---
cat(paste("Species:", config$species, "\n"))

if (!"databases" %in% names(config) || !config$species %in% names(config$databases)) {
  stop("[FATAL] 'databases' section or species entry missing in config.")
}
species_info <- config$databases[[config$species]]

if (!"organism_db" %in% names(species_info)) {
  stop("[FATAL] 'organism_db' key missing under species entry in config.")
}
organism_db_name <- species_info$organism_db
cat(paste("Loading organism database:", organism_db_name, "\n"))
if (!require(organism_db_name, character.only = TRUE, quietly = TRUE)) {
  stop(paste("[FATAL] Required organism DB package", organism_db_name, "is not installed."))
}
organism_db <- get(organism_db_name)
cat("  ✓ Organism database loaded successfully\n\n")

# --- 2c. Function to convert Entrez IDs to Gene Symbols ---
convert_entrez_to_symbols <- function(entrez_ids_string, organism_db) {
  # entrez_ids_string: e.g., "1234/5678/9012"
  # Returns: e.g., "GENEA/GENEB/GENEC"
  
  if (is.na(entrez_ids_string) || entrez_ids_string == "" || is.null(entrez_ids_string)) {
    return("")
  }
  
  # Split Entrez IDs
  entrez_ids <- strsplit(as.character(entrez_ids_string), "/")[[1]]
  
  # Convert to symbols
  symbols <- tryCatch({
    mapIds(organism_db, 
           keys = entrez_ids,
           column = "SYMBOL",
           keytype = "ENTREZID",
           multiVals = "first")
  }, error = function(e) {
    # If conversion fails, return original IDs
    return(entrez_ids)
  })
  
  # Replace NAs with original Entrez IDs
  symbols[is.na(symbols)] <- entrez_ids[is.na(symbols)]
  
  # Combine with "/"
  return(paste(symbols, collapse = "/"))
}

# --- 3a. GO file collection function ---
collect_go_results <- function(output_dir, gene_set, ontology) {
  # Reads GO enrichment results for specified gene_set (up/down/total)
  # and ontology (BP/CC/MF)
  
  file_path <- file.path(output_dir, 
                         paste0("go_enrichment_", gene_set, "_", ontology, ".csv"))
  
  if (!file.exists(file_path)) {
    cat(paste("  [WARNING] File not found:", basename(file_path), "\n"))
    return(NULL)
  }
  
  # Skip if file contains only empty/quoted string (no real data)
  raw_content <- trimws(paste(readLines(file_path, warn = FALSE), collapse = ""))
  if (raw_content == "" || raw_content == '""') {
    cat(paste("  [INFO] Empty results for", gene_set, ontology, "\n"))
    return(NULL)
  }

  # Read GO enrichment results
  go_res <- read.csv(file_path, stringsAsFactors = FALSE)

  # Skip if empty
  if (nrow(go_res) == 0) {
    cat(paste("  [INFO] Empty results for", gene_set, ontology, "\n"))
    return(NULL)
  }
  
  # Add metadata columns
  go_res$GeneSet <- toupper(gene_set)  # UP, DOWN, TOTAL
  go_res$Ontology <- ontology           # BP, CC, MF
  
  cat(paste("  ✓ Loaded", nrow(go_res), "terms from", gene_set, ontology, "\n"))
  
  return(go_res)
}

# --- 3b. KEGG file collection function ---
collect_kegg_results <- function(output_dir, gene_set) {
  # Reads KEGG enrichment results for specified gene_set (up/down/total)
  
  file_path <- file.path(output_dir, 
                         paste0("kegg_enrichment_", gene_set, ".csv"))
  
  if (!file.exists(file_path)) {
    cat(paste("  [WARNING] File not found:", basename(file_path), "\n"))
    return(NULL)
  }

  # Skip if file contains only empty/quoted string (no real data)
  raw_content <- trimws(paste(readLines(file_path, warn = FALSE), collapse = ""))
  if (raw_content == "" || raw_content == '""') {
    cat(paste("  [INFO] Empty KEGG results for", gene_set, "\n"))
    return(NULL)
  }

  # Read KEGG enrichment results
  kegg_res <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Skip if empty
  if (nrow(kegg_res) == 0) {
    cat(paste("  [INFO] Empty KEGG results for", gene_set, "\n"))
    return(NULL)
  }
  
  # Add metadata column
  kegg_res$GeneSet <- toupper(gene_set)  # UP, DOWN, TOTAL
  
  cat(paste("  ✓ Loaded", nrow(kegg_res), "pathways from", gene_set, "\n"))
  
  return(kegg_res)
}

# --- 4. Collect all GO and KEGG results ---
cat("Collecting GO enrichment results...\n")

gene_sets <- c("up", "down", "total")
ontologies <- c("BP", "CC", "MF")

all_go_results <- list()

for (geneset in gene_sets) {
  for (ont in ontologies) {
    go_data <- collect_go_results(output_dir, geneset, ont)
    if (!is.null(go_data)) {
      all_go_results[[paste(geneset, ont, sep="_")]] <- go_data
    }
  }
}

# Collect KEGG results
cat("\nCollecting KEGG enrichment results...\n")
all_kegg_results <- list()

for (geneset in gene_sets) {
  kegg_data <- collect_kegg_results(output_dir, geneset)
  if (!is.null(kegg_data)) {
    all_kegg_results[[geneset]] <- kegg_data
  }
}

# 결과가 하나도 없으면 종료
if (length(all_go_results) == 0 && length(all_kegg_results) == 0) {
  cat("\n[WARNING] No GO or KEGG enrichment results found. Skipping table generation.\n")
  quit(save = "no", status = 0)
}

# --- 5. Data formatting ---
cat("\nFormatting results for publication...\n")

format_go_table <- function(go_df, organism_db) {
  # Formats GO enrichment results for publication
  
  # Convert Entrez IDs to Gene Symbols
  cat("  Converting Entrez IDs to Gene Symbols for GO...\n")
  go_df$geneSymbol <- sapply(go_df$geneID, function(x) {
    convert_entrez_to_symbols(x, organism_db)
  })
  cat("  ✓ Conversion complete\n")
  
  # Select and order columns (use dplyr::select explicitly to avoid conflicts with AnnotationDbi)
  formatted <- go_df %>%
    dplyr::select(
      GeneSet,
      Ontology,
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count,
      geneSymbol
    ) %>%
    # 컬럼명을 논문 친화적으로 변경
    dplyr::rename(
      `Gene Set` = GeneSet,
      `GO ID` = ID,
      `GO Term` = Description,
      `Gene Ratio` = GeneRatio,
      `Background Ratio` = BgRatio,
      `P-value` = pvalue,
      `Adjusted P-value` = p.adjust,
      `Q-value` = qvalue,
      `Gene Count` = Count,
      `Gene Symbols` = geneSymbol
    ) %>%
    # p-value로 정렬 (가장 유의한 것부터)
    dplyr::arrange(Ontology, `Gene Set`, `Adjusted P-value`)
  
  return(formatted)
}

format_kegg_table <- function(kegg_df, organism_db) {
  # Formats KEGG enrichment results for publication
  
  # Convert Entrez IDs to Gene Symbols
  cat("  Converting Entrez IDs to Gene Symbols for KEGG...\n")
  kegg_df$geneSymbol <- sapply(kegg_df$geneID, function(x) {
    convert_entrez_to_symbols(x, organism_db)
  })
  cat("  ✓ Conversion complete\n")
  
  # Select and order columns
  formatted <- kegg_df %>%
    dplyr::select(
      GeneSet,
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count,
      geneSymbol
    ) %>%
    # 컬럼명을 논문 친화적으로 변경
    dplyr::rename(
      `Gene Set` = GeneSet,
      `KEGG ID` = ID,
      `KEGG Pathway` = Description,
      `Gene Ratio` = GeneRatio,
      `Background Ratio` = BgRatio,
      `P-value` = pvalue,
      `Adjusted P-value` = p.adjust,
      `Q-value` = qvalue,
      `Gene Count` = Count,
      `Gene Symbols` = geneSymbol
    ) %>%
    # p-value로 정렬 (가장 유의한 것부터)
    dplyr::arrange(`Gene Set`, `Adjusted P-value`)
  
  return(formatted)
}

# Format GO results
formatted_go <- NULL
if (length(all_go_results) > 0) {
  combined_go <- bind_rows(all_go_results)
  formatted_go <- format_go_table(combined_go, organism_db)
  cat(paste("  Total GO terms collected:", nrow(formatted_go), "\n"))
}

# Format KEGG results
formatted_kegg <- NULL
if (length(all_kegg_results) > 0) {
  combined_kegg <- bind_rows(all_kegg_results)
  formatted_kegg <- format_kegg_table(combined_kegg, organism_db)
  cat(paste("  Total KEGG pathways collected:", nrow(formatted_kegg), "\n"))
}

# --- 6. Create Excel file (organize by worksheet) ---
cat("\nCreating Excel workbook...\n")

# Create Excel workbook
wb <- createWorkbook()

# Define styles
header_style <- createStyle(
  fontSize = 11,
  fontName = "Arial",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  fgFill = "#4472C4",
  fontColour = "#FFFFFF",
  border = "TopBottomLeftRight",
  borderColour = "#000000",
  wrapText = FALSE
)

text_style <- createStyle(
  fontSize = 10,
  fontName = "Arial",
  halign = "left",
  valign = "center",
  border = "TopBottomLeftRight",
  borderColour = "#D3D3D3",
  wrapText = FALSE
)

number_style <- createStyle(
  fontSize = 10,
  fontName = "Arial",
  halign = "right",
  valign = "center",
  border = "TopBottomLeftRight",
  borderColour = "#D3D3D3",
  numFmt = "0.00E+00"  # Scientific notation
)

pvalue_style <- createStyle(
  fontSize = 10,
  fontName = "Arial",
  halign = "right",
  valign = "center",
  border = "TopBottomLeftRight",
  borderColour = "#D3D3D3",
  numFmt = "0.000"
)

# --- 6a. Create sheets for GO and KEGG results ---

# Sheet 1: GO Results by Gene Set and Ontology
if (!is.null(formatted_go)) {
  cat("  Creating GO enrichment sheets...\n")
  
  for (gs in c("UP", "DOWN", "TOTAL")) {
    for (ont in c("BP", "CC", "MF")) {
      sheet_data <- formatted_go %>% filter(`Gene Set` == gs, Ontology == ont)
      
      if (nrow(sheet_data) > 0) {
        sheet_name <- paste0(gs, "_", ont)
        cat(paste("    Adding sheet:", sheet_name, "\n"))
        
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, sheet_data, startRow = 1, startCol = 1, headerStyle = header_style)
        
        # Apply header style
        addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(sheet_data), gridExpand = TRUE)
        
        # Apply column-specific styles
        text_cols <- c(1, 2, 3, 4, 5, 6, 11)  # GeneSet, Ontology, GO ID, GO Term, Ratios, Gene Symbols
        addStyle(wb, sheet_name, text_style, rows = 2:(nrow(sheet_data) + 1), cols = text_cols, gridExpand = TRUE)
        
        # P-value columns
        pvalue_cols <- c(7, 8, 9)  # P-value, Adjusted P-value, Q-value
        addStyle(wb, sheet_name, pvalue_style, rows = 2:(nrow(sheet_data) + 1), cols = pvalue_cols, gridExpand = TRUE)
        
        # Gene Count column
        count_col <- 10
        addStyle(wb, sheet_name, text_style, rows = 2:(nrow(sheet_data) + 1), cols = count_col, gridExpand = TRUE)
        
        # Auto-size columns
        setColWidths(wb, sheet_name, cols = 1:ncol(sheet_data), widths = "auto")
      }
    }
  }
}

# Sheet 2: KEGG Results by Gene Set
if (!is.null(formatted_kegg)) {
  cat("  Creating KEGG enrichment sheets...\n")
  
  for (gs in c("UP", "DOWN", "TOTAL")) {
    sheet_data <- formatted_kegg %>% filter(`Gene Set` == gs)
    
    if (nrow(sheet_data) > 0) {
      sheet_name <- paste0("KEGG_", gs)
      cat(paste("    Adding sheet:", sheet_name, "\n"))
      
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, sheet_data, startRow = 1, startCol = 1, headerStyle = header_style)
      
      # Apply header style
      addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(sheet_data), gridExpand = TRUE)
      
      # Apply column-specific styles
      text_cols <- c(1, 2, 3, 4, 5, 10)  # GeneSet, KEGG ID, Pathway, Ratios, Gene Symbols
      addStyle(wb, sheet_name, text_style, rows = 2:(nrow(sheet_data) + 1), cols = text_cols, gridExpand = TRUE)
      
      # P-value columns
      pvalue_cols <- c(6, 7, 8)  # P-value, Adjusted P-value, Q-value
      addStyle(wb, sheet_name, pvalue_style, rows = 2:(nrow(sheet_data) + 1), cols = pvalue_cols, gridExpand = TRUE)
      
      # Gene Count column
      count_col <- 9
      addStyle(wb, sheet_name, text_style, rows = 2:(nrow(sheet_data) + 1), cols = count_col, gridExpand = TRUE)
      
      # Auto-size columns
      setColWidths(wb, sheet_name, cols = 1:ncol(sheet_data), widths = "auto")
    }
  }
}

# --- 6b. Old code removal: No longer creating "All_Results" sheet ---
# The old "All_Results" sheet has been removed to avoid confusion
# Users can refer to individual sheets for each gene set and ontology combination

# --- 6c. Create Analysis Info sheet ---
cat("  Creating Analysis Info sheet...\n")

# Count GO terms
go_total <- 0
go_up <- 0
go_down <- 0
go_total_all <- 0

if (!is.null(formatted_go)) {
  go_up <- formatted_go %>% filter(`Gene Set` == "UP") %>% nrow()
  go_down <- formatted_go %>% filter(`Gene Set` == "DOWN") %>% nrow()
  go_total_all <- formatted_go %>% filter(`Gene Set` == "TOTAL") %>% nrow()
  go_total <- nrow(formatted_go)
}

# Count KEGG pathways
kegg_total <- 0
kegg_up <- 0
kegg_down <- 0
kegg_total_all <- 0

if (!is.null(formatted_kegg)) {
  kegg_up <- formatted_kegg %>% filter(`Gene Set` == "UP") %>% nrow()
  kegg_down <- formatted_kegg %>% filter(`Gene Set` == "DOWN") %>% nrow()
  kegg_total_all <- formatted_kegg %>% filter(`Gene Set` == "TOTAL") %>% nrow()
  kegg_total <- nrow(formatted_kegg)
}

# Create analysis info data frame
analysis_info <- data.frame(
  Parameter = c(
    "Comparison",
    "Analysis Date",
    "DE Method",
    "P-value Cutoff (padj)",
    "Log2FC Cutoff",
    "GO P-value Cutoff",
    "GO Q-value Cutoff",
    "Species",
    "Organism Database",
    "Total GO Terms Found",
    "UP-regulated GO Terms",
    "DOWN-regulated GO Terms",
    "TOTAL GO Terms",
    "Total KEGG Pathways Found",
    "UP-regulated KEGG Pathways",
    "DOWN-regulated KEGG Pathways",
    "TOTAL KEGG Pathways"
  ),
  Value = c(
    paste(compare_group, "vs", base_group),
    format(Sys.Date(), "%Y-%m-%d"),
    ifelse(is.null(config$de_analysis$method), "DESeq2", config$de_analysis$method),
    ifelse(is.null(config$de_analysis$padj_cutoff), 0.05, config$de_analysis$padj_cutoff),
    ifelse(is.null(config$de_analysis$log2fc_cutoff), 0.5, config$de_analysis$log2fc_cutoff),
    ifelse(is.null(config$enrichment$pvalue_cutoff), 0.05, config$enrichment$pvalue_cutoff),
    ifelse(is.null(config$enrichment$qvalue_cutoff), 0.25, config$enrichment$qvalue_cutoff),
    config$species,
    organism_db_name,
    go_total,
    go_up,
    go_down,
    go_total_all,
    kegg_total,
    kegg_up,
    kegg_down,
    kegg_total_all
  ),
  stringsAsFactors = FALSE
)

# Add Analysis Info sheet at the beginning
addWorksheet(wb, "Analysis_Info", gridLines = TRUE)
writeData(wb, "Analysis_Info", analysis_info, startRow = 1, startCol = 1)

# Style the Analysis Info sheet
info_header_style <- createStyle(
  fontSize = 12,
  fontColour = "#FFFFFF",
  halign = "center",
  fgFill = "#4472C4",
  border = "TopBottomLeftRight",
  borderColour = "#000000",
  textDecoration = "bold"
)

info_text_style <- createStyle(
  fontSize = 11,
  halign = "left",
  valign = "center",
  border = "TopBottomLeftRight",
  borderColour = "#CCCCCC",
  wrapText = FALSE
)

# Apply styles
addStyle(wb, "Analysis_Info", info_header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
addStyle(wb, "Analysis_Info", info_text_style, rows = 2:(nrow(analysis_info) + 1), cols = 1:2, gridExpand = TRUE)

# Set column widths
setColWidths(wb, "Analysis_Info", cols = 1, widths = 30)
setColWidths(wb, "Analysis_Info", cols = 2, widths = 25)

cat("  ✓ Analysis Info sheet created\n")

# --- 7. Save Excel file ---
output_file <- file.path(output_dir, "final_go_results.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\n==============================================\n")
cat("✓ Enrichment Summary Table successfully generated!\n")
cat("==============================================\n")
cat(paste("Output file:", output_file, "\n"))
cat(paste("Total sheets:", length(names(wb)), "\n"))
cat("\nSheet organization:\n")
cat("  Analysis Information:\n")
cat("    • Analysis_Info (parameters and summary statistics)\n\n")
if (!is.null(formatted_go)) {
  cat("  GO Enrichment sheets:\n")
  for (gs in c("UP", "DOWN", "TOTAL")) {
    for (ont in c("BP", "CC", "MF")) {
      cat(paste("    •", paste0(gs, "_", ont), "\n"))
    }
  }
}
if (!is.null(formatted_kegg)) {
  cat("  KEGG Pathway sheets:\n")
  for (gs in c("UP", "DOWN", "TOTAL")) {
    cat(paste("    • KEGG_", gs, "\n", sep=""))
  }
}
cat("\nThis file is ready for supplementary material!\n\n")
