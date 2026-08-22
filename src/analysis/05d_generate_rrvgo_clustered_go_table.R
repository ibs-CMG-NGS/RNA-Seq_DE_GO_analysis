# File: src/analysis/05d_generate_rrvgo_clustered_go_table.R
#
# Purpose: Aggregate the per-(geneset, ontology) go_rrvgo_*.csv files
# (produced by 03_enrichment_analysis.R's rrvgo semantic-reduction step) into
# a single CMG-SeqViewer "clustered GO" Excel file, following the same
# column contract as 05b_generate_clustered_go_table.R
# (docs/RRVGO_CLUSTERED_GO_FORMAT.md), plus two rrvgo-specific columns:
#   - one sheet per Gene Set x Ontology group
#   - cluster_id as a zero-padded digit string, globally unique across the
#     whole file (not just within one group)
#   - singleton terms (final cluster size 1) get cluster_id = "Singleton"
#   - adjusted-p column must be named "Adjusted P-value" (never "padj")
#   - "Representative Term" = rrvgo's parentTerm (semantic cluster label)
#   - "Algorithm" = fixed "Semantic (rrvgo)" (distinguishes this file's
#     clustering method from 05b's Jaccard/gene-overlap clustering)
#
# go_rrvgo_*.csv (rrvgo::reduceSimMatrix() output) does not carry
# GeneRatio/BgRatio/pvalue/Count/geneID, so this script joins it back to
# go_enrichment_*.csv (the full enrichResult dump, written for every gene
# set/ontology by 03_enrichment_analysis.R) on GO ID.
#
# Usage: Rscript 05d_generate_rrvgo_clustered_go_table.R [config_path] [compare_group] [base_group] [output_dir]

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
  library(openxlsx)
  library(AnnotationDbi)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 05d_generate_rrvgo_clustered_go_table.R [config_path] [compare_group] [base_group] [output_dir]")
}
config_path   <- args[1]
compare_group <- args[2]
base_group    <- args[3]
output_dir    <- args[4]

config <- yaml.load_file(config_path)

rr_cfg <- config$enrichment$rrvgo
rr_enabled <- if (is.null(rr_cfg) || is.null(rr_cfg$enabled)) TRUE else isTRUE(rr_cfg$enabled)

if (!rr_enabled) {
  cat("[05d_generate_rrvgo_clustered_go_table] rrvgo.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

cat(paste("[05d_generate_rrvgo_clustered_go_table]", compare_group, "vs", base_group, "\n"))

# --- 2. Organism DB 로드 (Entrez -> Symbol 변환용) ---
species_info <- config$databases[[config$species]]
organism_db_name <- species_info$organism_db
if (!require(organism_db_name, character.only = TRUE, quietly = TRUE)) {
  stop(paste("[FATAL] Required organism DB package", organism_db_name, "is not installed."))
}
organism_db <- get(organism_db_name)

convert_entrez_column_to_symbols <- function(entrez_ids_strings, organism_db) {
  non_empty <- entrez_ids_strings[!is.na(entrez_ids_strings) & entrez_ids_strings != ""]
  all_ids <- unique(unlist(strsplit(as.character(non_empty), "/")))
  lookup <- tryCatch({
    mapIds(organism_db, keys = all_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  }, error = function(e) setNames(all_ids, all_ids))
  vapply(entrez_ids_strings, function(x) {
    if (is.na(x) || x == "") return("")
    ids <- strsplit(as.character(x), "/")[[1]]
    symbols <- lookup[ids]
    symbols[is.na(symbols)] <- ids[is.na(symbols)]
    paste(symbols, collapse = "/")
  }, character(1), USE.NAMES = FALSE)
}

# --- 3. 그룹별 go_rrvgo_*.csv 수집 + cluster_id 전역 유일화 ---
gene_sets  <- c("up", "down")
ontologies <- config$enrichment$go_ontologies %||% c("BP", "CC", "MF")

# 그룹(Gene Set x Ontology) 내부의 로컬 cluster 번호를 전역 zero-padded 문자열로 재부여.
# 최종 클러스터 크기가 1인 term은 "Singleton"으로 표시(05b와 동일 규칙).
relabel_clusters <- function(cluster_vec, start_counter) {
  sizes <- table(cluster_vec)
  multi_ids <- names(sizes)[sizes > 1]
  multi_ids <- multi_ids[order(suppressWarnings(as.numeric(multi_ids)))]
  mapping <- setNames(sprintf("%03d", start_counter + seq_along(multi_ids) - 1L), multi_ids)
  new_id <- ifelse(as.character(cluster_vec) %in% multi_ids,
                    mapping[as.character(cluster_vec)],
                    "Singleton")
  list(ids = unname(new_id), next_counter = start_counter + length(multi_ids))
}

wb <- createWorkbook()
header_style <- createStyle(fontSize = 11, fontName = "Arial", textDecoration = "bold",
                             halign = "center", valign = "center", fgFill = "#4472C4",
                             fontColour = "#FFFFFF", border = "TopBottomLeftRight", borderColour = "#000000")

global_counter <- 1L
sheet_summary <- list()

for (gs in gene_sets) {
  for (ont in ontologies) {
    in_csv <- file.path(output_dir, paste0("go_rrvgo_", gs, "_", ont, ".csv"))
    if (!file.exists(in_csv)) next
    d <- read.csv(in_csv, stringsAsFactors = FALSE)
    if (nrow(d) == 0) next

    # go_rrvgo_*.csv(reduceSimMatrix 출력)에는 GeneRatio/BgRatio/pvalue/Count/geneID가
    # 없으므로, 같은 gene set/ontology의 go_enrichment_*.csv(원본 enrichResult 전체)와
    # GO ID 기준으로 join해서 채운다.
    enrich_csv <- file.path(output_dir, paste0("go_enrichment_", gs, "_", ont, ".csv"))
    if (!file.exists(enrich_csv)) {
      cat(sprintf("[05d_generate_rrvgo_clustered_go_table] %s missing — skipping %s/%s.\n",
                   basename(enrich_csv), gs, ont))
      next
    }
    e <- read.csv(enrich_csv, stringsAsFactors = FALSE)
    d <- merge(d, e[, c("ID", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "geneID")],
               by.x = "go", by.y = "ID", all.x = TRUE)
    if (any(is.na(d$Count))) {
      cat(sprintf("[05d_generate_rrvgo_clustered_go_table] %d/%d rrvgo term(s) not found in %s — dropped.\n",
                   sum(is.na(d$Count)), nrow(d), basename(enrich_csv)))
      d <- d[!is.na(d$Count), ]
    }
    if (nrow(d) == 0) next

    relabeled <- relabel_clusters(d$cluster, global_counter)
    global_counter <- relabeled$next_counter

    d$geneSymbol <- convert_entrez_column_to_symbols(d$geneID, organism_db)

    d_out <- data.frame(
      `Gene Set`            = toupper(gs),
      `Ontology`            = ont,
      `GO ID`               = d$go,
      `GO Term`             = d$term,
      `Gene Ratio`          = d$GeneRatio,
      `Background Ratio`    = d$BgRatio,
      `P-value`             = d$pvalue,
      `Adjusted P-value`    = d$p.adjust,
      `Gene Count`          = d$Count,
      `Gene Symbols`        = d$geneSymbol,
      `cluster_id`          = relabeled$ids,
      `Representative Term` = d$parentTerm,
      `Algorithm`           = "Semantic (rrvgo)",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    d_out <- d_out[order(d_out$cluster_id == "Singleton", d_out$cluster_id, d_out$`Adjusted P-value`), ]

    sheet_name <- substr(paste0(toupper(gs), "_", ont), 1, 31)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, d_out, headerStyle = header_style)
    addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(d_out), gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = 1:ncol(d_out), widths = "auto")
    sheet_summary[[sheet_name]] <- nrow(d_out)

    n_clusters_final <- length(unique(d_out$cluster_id[d_out$cluster_id != "Singleton"]))
    n_singletons <- sum(d_out$cluster_id == "Singleton")
    cat(sprintf("  %s: %d terms, %d clusters, %d singletons\n", sheet_name, nrow(d_out), n_clusters_final, n_singletons))
  }
}

output_file <- file.path(output_dir, "final_go_rrvgo_clustered_results.xlsx")

if (length(sheet_summary) == 0) {
  # Snakemake는 이 rule이 output_file을 반드시 만들 것으로 기대하므로(선언된 output),
  # 유의 term이 하나도 없어 클러스터링할 게 없는 경우에도 빈 placeholder를 써야 한다.
  cat("[05d_generate_rrvgo_clustered_go_table] No go_rrvgo_*.csv files found — writing empty placeholder.\n")
  addWorksheet(wb, "No Results")
  writeData(wb, "No Results",
            data.frame(Message = "No significant GO terms available to cluster for this comparison."))
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat(paste("[05d_generate_rrvgo_clustered_go_table] Saved empty placeholder:", output_file, "\n"))
  quit(save = "no", status = 0)
}

saveWorkbook(wb, output_file, overwrite = TRUE)
cat(paste("[05d_generate_rrvgo_clustered_go_table] Saved:", output_file,
          "(", length(sheet_summary), "sheets,", global_counter - 1L, "clusters total )\n"))
