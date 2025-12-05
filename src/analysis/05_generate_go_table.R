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

# --- 3. GO file collection function ---
collect_go_results <- function(output_dir, gene_set, ontology) {
  # Reads GO enrichment results for specified gene_set (up/down/total)
  # and ontology (BP/CC/MF)
  
  file_path <- file.path(output_dir, 
                         paste0("go_enrichment_", gene_set, "_", ontology, ".csv"))
  
  if (!file.exists(file_path)) {
    cat(paste("  [WARNING] File not found:", basename(file_path), "\n"))
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

# --- 4. Collect all GO results ---
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

# 결과가 하나도 없으면 종료
if (length(all_go_results) == 0) {
  cat("\n[WARNING] No GO enrichment results found. Skipping GO table generation.\n")
  quit(save = "no", status = 0)
}

# --- 5. Data formatting ---
cat("\nFormatting results for publication...\n")

format_go_table <- function(go_df, organism_db) {
  # Formats GO enrichment results for publication
  
  # Convert Entrez IDs to Gene Symbols
  cat("  Converting Entrez IDs to Gene Symbols...\n")
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

# Combine all results into one dataframe
combined_go <- bind_rows(all_go_results)
formatted_go <- format_go_table(combined_go, organism_db)

cat(paste("  Total GO terms collected:", nrow(formatted_go), "\n"))

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

# --- 6a. Sheet 1: Summary (All results combined) ---
cat("  Creating 'All_Results' sheet...\n")
addWorksheet(wb, "All_Results")
writeData(wb, "All_Results", formatted_go, startRow = 1, startCol = 1, headerStyle = header_style)

# Apply header style
addStyle(wb, "All_Results", header_style, rows = 1, cols = 1:ncol(formatted_go), gridExpand = TRUE)

# Apply column-specific styles
if (nrow(formatted_go) > 0) {
  # Text columns
  text_cols <- c(1, 2, 3, 4, 5, 6, 11)  # GeneSet, Ontology, GO ID, GO Term, Ratios, Gene Symbols
  addStyle(wb, "All_Results", text_style, rows = 2:(nrow(formatted_go) + 1), cols = text_cols, gridExpand = TRUE)
  
  # P-value columns (scientific notation)
  pval_cols <- c(7, 8, 9)  # P-value, Adjusted P-value, Q-value
  addStyle(wb, "All_Results", number_style, rows = 2:(nrow(formatted_go) + 1), cols = pval_cols, gridExpand = TRUE)
  
  # Count column
  addStyle(wb, "All_Results", text_style, rows = 2:(nrow(formatted_go) + 1), cols = 10, gridExpand = TRUE)
}

# Auto-adjust column widths
setColWidths(wb, "All_Results", cols = 1:ncol(formatted_go), widths = "auto")
setColWidths(wb, "All_Results", cols = 4, widths = 50)    # GO Term
setColWidths(wb, "All_Results", cols = 11, widths = 60)   # Gene Symbols
freezePane(wb, "All_Results", firstRow = TRUE)

# --- 6b. Gene Set-specific sheets (UP, DOWN, TOTAL) ---
for (geneset in c("UP", "DOWN", "TOTAL")) {
  subset_data <- formatted_go %>% dplyr::filter(`Gene Set` == geneset)
  
  if (nrow(subset_data) > 0) {
    sheet_name <- paste0(geneset, "_regulated")
    cat(paste("  Creating", paste0("'", sheet_name, "'"), "sheet...\n"))
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, subset_data, startRow = 1, startCol = 1, headerStyle = header_style)
    
    # 스타일 적용
    addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(subset_data), gridExpand = TRUE)
    addStyle(wb, sheet_name, text_style, rows = 2:(nrow(subset_data) + 1), cols = text_cols, gridExpand = TRUE)
    addStyle(wb, sheet_name, number_style, rows = 2:(nrow(subset_data) + 1), cols = pval_cols, gridExpand = TRUE)
    addStyle(wb, sheet_name, text_style, rows = 2:(nrow(subset_data) + 1), cols = 10, gridExpand = TRUE)
    
    # 컬럼 너비 조정
    setColWidths(wb, sheet_name, cols = 1:ncol(subset_data), widths = "auto")
    setColWidths(wb, sheet_name, cols = 4, widths = 50)
    setColWidths(wb, sheet_name, cols = 11, widths = 60)
    freezePane(wb, sheet_name, firstRow = TRUE)
  }
}

# --- 6c. Ontology-specific sheets (BP, CC, MF) ---
for (ont in c("BP", "CC", "MF")) {
  subset_data <- formatted_go %>% dplyr::filter(Ontology == ont)
  
  if (nrow(subset_data) > 0) {
    ont_fullname <- switch(ont,
                           "BP" = "Biological_Process",
                           "CC" = "Cellular_Component",
                           "MF" = "Molecular_Function")
    
    cat(paste("  Creating", paste0("'", ont_fullname, "'"), "sheet...\n"))
    
    addWorksheet(wb, ont_fullname)
    writeData(wb, ont_fullname, subset_data, startRow = 1, startCol = 1, headerStyle = header_style)
    
    # Apply styles
    addStyle(wb, ont_fullname, header_style, rows = 1, cols = 1:ncol(subset_data), gridExpand = TRUE)
    addStyle(wb, ont_fullname, text_style, rows = 2:(nrow(subset_data) + 1), cols = text_cols, gridExpand = TRUE)
    addStyle(wb, ont_fullname, number_style, rows = 2:(nrow(subset_data) + 1), cols = pval_cols, gridExpand = TRUE)
    addStyle(wb, ont_fullname, text_style, rows = 2:(nrow(subset_data) + 1), cols = 10, gridExpand = TRUE)
    
    # Adjust column widths
    setColWidths(wb, ont_fullname, cols = 1:ncol(subset_data), widths = "auto")
    setColWidths(wb, ont_fullname, cols = 4, widths = 50)
    setColWidths(wb, ont_fullname, cols = 11, widths = 60)
    freezePane(wb, ont_fullname, firstRow = TRUE)
  }
}

# --- 6d. Top terms summary sheet (Most significant results) ---
cat("  Creating 'Top_Terms' summary sheet...\n")

top_n <- 20  # Top 20 from each category

top_terms <- formatted_go %>%
  group_by(`Gene Set`, Ontology) %>%
  dplyr::arrange(`Adjusted P-value`) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  dplyr::arrange(Ontology, `Gene Set`, `Adjusted P-value`)

if (nrow(top_terms) > 0) {
  addWorksheet(wb, "Top_Terms")
  writeData(wb, "Top_Terms", top_terms, startRow = 1, startCol = 1, headerStyle = header_style)
  
  # 스타일 적용
  addStyle(wb, "Top_Terms", header_style, rows = 1, cols = 1:ncol(top_terms), gridExpand = TRUE)
  addStyle(wb, "Top_Terms", text_style, rows = 2:(nrow(top_terms) + 1), cols = text_cols, gridExpand = TRUE)
  addStyle(wb, "Top_Terms", number_style, rows = 2:(nrow(top_terms) + 1), cols = pval_cols, gridExpand = TRUE)
  addStyle(wb, "Top_Terms", text_style, rows = 2:(nrow(top_terms) + 1), cols = 10, gridExpand = TRUE)
  
  setColWidths(wb, "Top_Terms", cols = 1:ncol(top_terms), widths = "auto")
  setColWidths(wb, "Top_Terms", cols = 4, widths = 50)
  setColWidths(wb, "Top_Terms", cols = 11, widths = 60)
  freezePane(wb, "Top_Terms", firstRow = TRUE)
}

# --- 6e. Metadata sheet (Analysis information) ---
cat("  Creating 'Analysis_Info' sheet...\n")

metadata <- data.frame(
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
    "UP-regulated Terms",
    "DOWN-regulated Terms",
    "TOTAL Terms"
  ),
  Value = c(
    paste(compare_group, "vs", base_group),
    format(Sys.Date(), "%Y-%m-%d"),
    config$de_analysis$method,
    config$de_analysis$padj_cutoff,
    config$de_analysis$log2fc_cutoff,
    config$enrichment$pvalue_cutoff,
    config$enrichment$qvalue_cutoff,
    config$species,
    config$databases[[config$species]]$organism_db,
    nrow(formatted_go),
    nrow(formatted_go %>% dplyr::filter(`Gene Set` == "UP")),
    nrow(formatted_go %>% dplyr::filter(`Gene Set` == "DOWN")),
    nrow(formatted_go %>% dplyr::filter(`Gene Set` == "TOTAL"))
  ),
  stringsAsFactors = FALSE
)

addWorksheet(wb, "Analysis_Info")
writeData(wb, "Analysis_Info", metadata, startRow = 1, startCol = 1, headerStyle = header_style)

# 스타일 적용
addStyle(wb, "Analysis_Info", header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
addStyle(wb, "Analysis_Info", text_style, rows = 2:(nrow(metadata) + 1), cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "Analysis_Info", cols = 1:2, widths = "auto")
setColWidths(wb, "Analysis_Info", cols = 1, widths = 30)

# --- 7. Save Excel file ---
output_file <- file.path(output_dir, "final_go_results.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\n==============================================\n")
cat("✓ GO Summary Table successfully generated!\n")
cat("==============================================\n")
cat(paste("Output file:", output_file, "\n"))
cat(paste("Total sheets:", length(names(wb)), "\n"))
cat("\nSheet contents:\n")
cat("  • All_Results: Complete GO enrichment results\n")
cat("  • UP_regulated: Up-regulated gene GO terms\n")
cat("  • DOWN_regulated: Down-regulated gene GO terms\n")
cat("  • TOTAL_regulated: All significant gene GO terms\n")
cat("  • Biological_Process: BP ontology terms\n")
cat("  • Cellular_Component: CC ontology terms\n")
cat("  • Molecular_Function: MF ontology terms\n")
cat("  • Top_Terms: Top 20 most significant terms per category\n")
cat("  • Analysis_Info: Analysis parameters and metadata\n")
cat("\nThis file is ready for supplementary material!\n\n")
