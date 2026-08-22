# 파일 경로: src/analysis/11_run_group_enrichment.R
# time_series 클러스터 / coexpression 모듈처럼 "그룹 컬럼(cluster_id, module_id 등)으로
# 나뉜 유의 유전자 리스트"에 대해 그룹별(+전체) GO/KEGG enrichment를 수행하고,
# 개별 CSV + final_go_results.xlsx(그룹×ontology 시트 구성)로 저장한다.
#
# 추가로 03_enrichment_analysis.R의 term_cluster(Jaccard)/rrvgo(의미론적 축약) 단계를
# 그룹이 가변적인 이 스크립트에 맞게 이식해서, 그룹별로도 too-many-significant-terms
# 문제를 요약해서 볼 수 있게 한다. pairwise와 달리 "all"(TOTAL, 전체 클러스터/모듈을
# 합친 유전자 집합)도 요약 대상에 포함한다 — pairwise의 "total"은 up+down 방향이
# 섞여서 제외하지만, 여기 TOTAL은 방향 충돌이 없는 "실험 전체 그림" 요약이라 유의미함.
#
# 03_enrichment_analysis.R(개별 GO/KEGG 실행 + term_cluster/rrvgo)과
# 05_generate_go_table.R / 05b_generate_clustered_go_table.R / 05d_generate_rrvgo_clustered_go_table.R
# (CSV -> xlsx)의 핵심 로직을 그룹이 가변적인 이 두 모듈(time_series/coexpression_modules)에
# 맞게 하나의 스크립트 안에서 재구성한 것.
#
# 사용법:
#   Rscript 11_run_group_enrichment.R <config_path> <input_csv> <group_col> <group_prefix> <output_dir>
#   예) Rscript 11_run_group_enrichment.R config.yml \
#         output/proj/time_series/time_series_significant_genes.csv \
#         cluster_id cluster output/proj/time_series

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(AnnotationDbi)
  library(dplyr)
  library(openxlsx)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript 11_run_group_enrichment.R <config_path> <input_csv> <group_col> <group_prefix> <output_dir>")
}
config_path <- args[1]
input_csv   <- args[2]
group_col   <- args[3]
group_prefix <- args[4]
output_dir  <- args[5]

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

config <- yaml.load_file(config_path)
enr_cfg <- config$enrichment %||% list()
go_ontologies  <- enr_cfg$go_ontologies %||% c("BP", "CC", "MF")
min_gs_size    <- enr_cfg$min_gs_size    %||% 10
max_gs_size    <- enr_cfg$max_gs_size    %||% 500
min_gene_count <- enr_cfg$min_gene_count %||% 1
max_genes_for_go <- 5000

species_info <- config$databases[[config$species]]
organism_db_name <- species_info$organism_db
cat(paste("[11_run_group_enrichment] Loading organism database:", organism_db_name, "\n"))
if (!require(organism_db_name, character.only = TRUE, quietly = TRUE)) {
  stop(paste("[FATAL] Required organism DB package", organism_db_name, "is not installed."))
}
organism_db <- get(organism_db_name)
kegg_organism <- species_info$kegg_code

# --- 2. 입력 CSV 로드 & 그룹 목록 구성 ---
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene_id" %in% colnames(df)) stop(paste("[FATAL]", input_csv, "does not have a 'gene_id' column."))
if (!group_col %in% colnames(df)) stop(paste("[FATAL]", input_csv, "does not have a", group_col, "column."))

group_values <- sort(unique(df[[group_col]]))
cat(paste("[11_run_group_enrichment] Found", length(group_values), "groups in", group_col, ":",
          paste(group_values, collapse = ", "), "\n"))

# 그룹 라벨 -> gene_id 벡터 (그룹별 + 전체("all"))
gene_lists <- setNames(
  lapply(group_values, function(v) df$gene_id[df[[group_col]] == v]),
  paste0(group_prefix, group_values)
)
gene_lists[["all"]] <- df$gene_id

# 그룹별 p-value 유사 지표(정렬/상위 컷용) — p_value 또는 padj 중 존재하는 것 사용
sort_col <- intersect(c("p_value", "padj"), colnames(df))[1]

# --- 3. Gene ID 타입 감지 (03_enrichment_analysis.R과 동일 로직) ---
detect_gene_id_type <- function(gene_ids) {
  if ("gene_id_type" %in% names(config)) return(toupper(config$gene_id_type))
  sample_id <- head(gene_ids[!is.na(gene_ids)], 1)
  if (grepl("^ENSMUSG[0-9]+", sample_id) || grepl("^ENSG[0-9]+", sample_id)) return("ENSEMBL")
  if (grepl("^[0-9]+$", sample_id)) return("ENTREZID")
  "SYMBOL"
}

convert_to_entrez <- function(gene_ids, gene_id_type) {
  if (gene_id_type == "ENTREZID") return(unique(as.character(gene_ids)))
  keytype <- if (gene_id_type == "SYMBOL") "SYMBOL" else "ENSEMBL"
  entrez <- tryCatch({
    mapIds(organism_db, keys = gene_ids, column = "ENTREZID", keytype = keytype, multiVals = "first")
  }, error = function(e) setNames(rep(NA_character_, length(gene_ids)), gene_ids))
  unique(as.character(entrez[!is.na(entrez)]))
}

write_empty_result <- function(path) {
  empty_df <- data.frame(ID = character(), Description = character(), GeneRatio = character(),
                          BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
                          qvalue = numeric(), geneID = character(), Count = integer())
  write.csv(empty_df, path, row.names = FALSE)
}

# --- 4. 그룹별 GO/KEGG enrichment 실행 ---
for (label in names(gene_lists)) {
  gene_ids <- gene_lists[[label]]
  cat(paste("\n[11_run_group_enrichment] Group:", label, "-", length(gene_ids), "genes\n"))

  gene_id_type <- detect_gene_id_type(gene_ids)
  entrez_ids <- convert_to_entrez(gene_ids, gene_id_type)

  if (length(entrez_ids) > max_genes_for_go && !is.na(sort_col)) {
    ord <- order(df[[sort_col]][match(gene_ids, df$gene_id)])
    top_ids <- gene_ids[ord][seq_len(min(max_genes_for_go, length(gene_ids)))]
    entrez_ids <- convert_to_entrez(top_ids, gene_id_type)
    cat(paste("  Limited to top", length(entrez_ids), "genes by", sort_col, "for enrichment\n"))
  }

  # GO
  for (ont in go_ontologies) {
    out_csv <- file.path(output_dir, paste0("go_enrichment_", label, "_", ont, ".csv"))
    if (length(entrez_ids) == 0) {
      write_empty_result(out_csv)
      next
    }
    go_res <- tryCatch({
      enrichGO(gene = entrez_ids, OrgDb = organism_db, keyType = "ENTREZID", ont = ont,
               pAdjustMethod = "BH", pvalueCutoff = 1.0, qvalueCutoff = 1.0,
               readable = FALSE, pool = FALSE,
               minGSSize = min_gs_size, maxGSSize = max_gs_size)
    }, error = function(e) {
      cat(paste("  [WARN] enrichGO failed for", label, ont, ":", conditionMessage(e), "\n"))
      NULL
    })
    if (is.null(go_res) || nrow(go_res) == 0) {
      write_empty_result(out_csv)
    } else {
      go_df <- as.data.frame(go_res) %>% filter(Count >= min_gene_count)
      write.csv(go_df, out_csv, row.names = FALSE)
      cat(paste("  Saved", nrow(go_df), "GO(", ont, ") terms\n"))
    }
  }

  # KEGG
  out_kegg_csv <- file.path(output_dir, paste0("kegg_enrichment_", label, ".csv"))
  if (length(entrez_ids) == 0) {
    write_empty_result(out_kegg_csv)
  } else {
    kegg_res <- tryCatch({
      enrichKEGG(gene = entrez_ids, organism = kegg_organism, pvalueCutoff = 1.0,
                 minGSSize = min_gs_size, maxGSSize = max_gs_size)
    }, error = function(e) {
      cat(paste("  [WARN] enrichKEGG failed for", label, ":", conditionMessage(e), "\n"))
      NULL
    })
    if (is.null(kegg_res) || nrow(kegg_res) == 0) {
      write_empty_result(out_kegg_csv)
    } else {
      kegg_df <- as.data.frame(kegg_res) %>% filter(Count >= min_gene_count)
      write.csv(kegg_df, out_kegg_csv, row.names = FALSE)
      cat(paste("  Saved", nrow(kegg_df), "KEGG pathways\n"))
    }
  }
}

# --- 4.5 그룹별 term_cluster(Jaccard) + rrvgo(의미론적 축약) ---
# 03_enrichment_analysis.R과 동일한 방식으로, 방금 저장한 go_enrichment_{label}_{ont}.csv를
# 다시 읽어 FDR/FoldEnrichment로 유의 term만 추리고 두 방식으로 각각 클러스터링한다.
# 이 단계에서는 (03과 동일하게) "로컬" cluster 번호만 CSV에 저장하고, 전역 유일 cluster_id
# 부여는 아래 6/7단계(xlsx 집계 시점)에서 처리한다.
parse_ratio <- function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1] / parts[2]
}
compute_fold_enrichment <- function(go_df) {
  mapply(function(gr, br) parse_ratio(gr) / parse_ratio(br), go_df$GeneRatio, go_df$BgRatio)
}

tc_cfg <- enr_cfg$term_cluster
tc_enabled <- if (is.null(tc_cfg) || is.null(tc_cfg$enabled)) TRUE else isTRUE(tc_cfg$enabled)
rr_cfg <- enr_cfg$rrvgo
rr_enabled <- if (is.null(rr_cfg) || is.null(rr_cfg$enabled)) TRUE else isTRUE(rr_cfg$enabled)

if (tc_enabled || rr_enabled) {
  cat("\n[11_run_group_enrichment] Running term_cluster/rrvgo per group...\n")

  tc_fdr_cutoff  <- ifelse(is.null(tc_cfg$fdr_cutoff), 0.05, tc_cfg$fdr_cutoff)
  tc_fe_cutoff   <- ifelse(is.null(tc_cfg$fold_enrichment_cutoff), 2.0, tc_cfg$fold_enrichment_cutoff)
  tc_similarity_cutoff <- if (!is.null(tc_cfg) && "similarity_cutoff" %in% names(tc_cfg)) tc_cfg$similarity_cutoff else 0.7
  tc_n_clusters  <- ifelse(is.null(tc_cfg$n_clusters), 5, tc_cfg$n_clusters)
  tc_hclust_method <- ifelse(is.null(tc_cfg$hclust_method), "average", tc_cfg$hclust_method)
  tc_show_top_n  <- ifelse(is.null(tc_cfg$show_top_n), 30, tc_cfg$show_top_n)
  tc_min_terms   <- 3

  rr_fdr_cutoff  <- ifelse(is.null(rr_cfg$fdr_cutoff), 0.05, rr_cfg$fdr_cutoff)
  rr_fe_cutoff   <- ifelse(is.null(rr_cfg$fold_enrichment_cutoff), 2.0, rr_cfg$fold_enrichment_cutoff)
  rr_method      <- ifelse(is.null(rr_cfg$method), "Rel", rr_cfg$method)
  rr_threshold   <- ifelse(is.null(rr_cfg$threshold), 0.7, rr_cfg$threshold)
  rr_score_by    <- ifelse(is.null(rr_cfg$score_by), "fdr", rr_cfg$score_by)
  rr_min_terms   <- 2

  for (label in names(gene_lists)) {
    for (ont in go_ontologies) {
      enrich_csv <- file.path(output_dir, paste0("go_enrichment_", label, "_", ont, ".csv"))
      if (!file.exists(enrich_csv)) next
      go_df <- tryCatch(read.csv(enrich_csv, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(go_df) || nrow(go_df) == 0) next
      go_df$FoldEnrichment <- compute_fold_enrichment(go_df)

      # --- term_cluster (Jaccard) ---
      if (tc_enabled) {
        sig_ids <- go_df$ID[!is.na(go_df$p.adjust) & go_df$p.adjust < tc_fdr_cutoff &
                             !is.na(go_df$FoldEnrichment) & go_df$FoldEnrichment > tc_fe_cutoff]
        out_tc_csv <- file.path(output_dir, paste0("go_termcluster_", label, "_", ont, ".csv"))
        out_tc_png <- file.path(output_dir, paste0("go_termcluster_", label, "_", ont, ".png"))

        if (length(sig_ids) >= tc_min_terms) {
          go_sig_result <- go_df[go_df$ID %in% sig_ids, ]
          # geneInCategory()/get_similarity_matrix()는 rownames(result)로 geneSets를 인덱싱하는
          # 게 아니라 result$ID 값 자체로 매칭하지만, geneInCategory()는 리스트 이름을
          # rownames(x@result)에서 가져온다 — read.csv()로 되읽은 data.frame은 기본 정수
          # rownames이라 명시적으로 ID와 맞춰주지 않으면 매칭이 깨져 유사도가 전부 0(=전부
          # singleton)이 되는 조용한 실패가 난다.
          rownames(go_sig_result) <- go_sig_result$ID
          go_sig <- new("enrichResult", result = go_sig_result,
                         pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1,
                         organism = "unknown", ontology = ont, gene = as.character(sig_ids),
                         keytype = "ENTREZID", universe = character(),
                         gene2Symbol = character(), geneSets = list())
          clus_result <- tryCatch({
            go_sig_termsim <- pairwise_termsim(go_sig, method = "JC", showCategory = length(sig_ids))
            keep <- seq_len(length(sig_ids))
            termsim2 <- go_sig_termsim@termsim[keep, keep]
            termsim2[is.na(termsim2)] <- 0
            termsim2 <- termsim2 + t(termsim2)
            diag(termsim2) <- 1
            hc_manual <- stats::hclust(stats::as.dist(1 - termsim2), method = tc_hclust_method)
            clus <- if (!is.null(tc_similarity_cutoff)) {
              stats::cutree(hc_manual, h = 1 - tc_similarity_cutoff)
            } else {
              stats::cutree(hc_manual, k = min(tc_n_clusters, length(sig_ids)))
            }
            list(termsim = go_sig_termsim, clus = clus)
          }, error = function(e) {
            cat(paste("[term_cluster]", label, ont, "failed:", e$message, "\n"))
            NULL
          })

          if (!is.null(clus_result)) {
            clus <- clus_result$clus
            effective_n <- length(unique(clus))
            cluster_df <- go_df[match(names(clus), go_df$Description),
                                 c("ID", "Description", "GeneRatio", "BgRatio", "FoldEnrichment",
                                   "pvalue", "p.adjust", "qvalue", "Count", "geneID")]
            cluster_df$cluster <- clus[cluster_df$Description]
            cluster_df <- cluster_df[order(cluster_df$cluster, cluster_df$p.adjust), ]
            write.csv(cluster_df, out_tc_csv, row.names = FALSE)
            cat(sprintf("[term_cluster] %s/%s: %d terms -> %d clusters\n", label, ont, length(sig_ids), effective_n))

            show_n <- min(as.integer(tc_show_top_n), length(sig_ids))
            tp <- tryCatch({
              treeplot(clus_result$termsim, showCategory = show_n,
                       cluster.params = list(method = tc_hclust_method, n = min(effective_n, show_n)))
            }, error = function(e) {
              cat(paste("[term_cluster] treeplot failed:", e$message, "\n"))
              NULL
            })
            if (!is.null(tp)) {
              plot_height <- max(8, show_n * 0.22, min(effective_n, show_n) * 0.4)
              ggsave(out_tc_png, plot = tp, width = 12, height = plot_height, bg = "white", limitsize = FALSE)
            }
          }
        }
      }

      # --- rrvgo (의미론적 축약) ---
      if (rr_enabled) {
        rr_sig_ids <- go_df$ID[!is.na(go_df$p.adjust) & go_df$p.adjust < rr_fdr_cutoff &
                                !is.na(go_df$FoldEnrichment) & go_df$FoldEnrichment > rr_fe_cutoff]
        out_rr_csv     <- file.path(output_dir, paste0("go_rrvgo_", label, "_", ont, ".csv"))
        out_rr_treemap <- file.path(output_dir, paste0("go_rrvgo_treemap_", label, "_", ont, ".png"))
        out_rr_scatter <- file.path(output_dir, paste0("go_rrvgo_scatter_", label, "_", ont, ".png"))

        if (length(rr_sig_ids) >= rr_min_terms) {
          suppressPackageStartupMessages(library(rrvgo))
          go_df_rr_sig <- go_df[go_df$ID %in% rr_sig_ids, ]
          rr_scores <- if (rr_score_by == "count") {
            setNames(go_df_rr_sig$Count, go_df_rr_sig$ID)
          } else {
            setNames(-log10(go_df_rr_sig$p.adjust), go_df_rr_sig$ID)
          }

          rr_result <- tryCatch({
            simMatrix <- calculateSimMatrix(rr_sig_ids, orgdb = organism_db_name, ont = ont, method = rr_method)
            reducedTerms <- reduceSimMatrix(simMatrix, scores = rr_scores, threshold = rr_threshold, orgdb = organism_db_name)
            list(simMatrix = simMatrix, reducedTerms = reducedTerms)
          }, error = function(e) {
            cat(paste("[rrvgo]", label, ont, "failed:", e$message, "\n"))
            NULL
          })

          if (!is.null(rr_result)) {
            reducedTerms <- rr_result$reducedTerms
            write.csv(reducedTerms, out_rr_csv, row.names = FALSE)
            cat(sprintf("[rrvgo] %s/%s: %d terms -> %d parent groups\n",
                        label, ont, length(rr_sig_ids), length(unique(reducedTerms$parent))))

            sp <- tryCatch(scatterPlot(rr_result$simMatrix, reducedTerms), error = function(e) NULL)
            if (!is.null(sp)) ggsave(out_rr_scatter, plot = sp, width = 10, height = 8, bg = "white")

            tm_ok <- tryCatch({
              png(out_rr_treemap, width = 12, height = 8, units = "in", res = 300, bg = "white")
              treemapPlot(reducedTerms)
              dev.off()
              TRUE
            }, error = function(e) {
              if (dev.cur() > 1) dev.off()
              FALSE
            })
          }
        }
      }
    }
  }
}

# --- 5. CSV 취합 -> final_go_results.xlsx (05_generate_go_table.R 패턴 재사용) ---
cat("\n[11_run_group_enrichment] Building final_go_results.xlsx...\n")

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

read_csv_if_nonempty <- function(path) {
  if (!file.exists(path)) return(NULL)
  raw <- trimws(paste(readLines(path, warn = FALSE), collapse = ""))
  if (raw == "" || raw == '""') return(NULL)
  d <- read.csv(path, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(NULL)
  d
}

# 05_generate_go_table.R(DEG 기반 final_go_results.xlsx)과 동일한 형식으로 맞추기 위한
# Gene Set 라벨 변환: "all" -> TOTAL, "cluster1"/"module2" -> "Cluster01"/"Module02"
format_gene_set_label <- function(label, group_prefix) {
  if (label == "all") return("TOTAL")
  num <- suppressWarnings(as.integer(sub(paste0("^", group_prefix), "", label)))
  prefix_title <- paste0(toupper(substr(group_prefix, 1, 1)), substr(group_prefix, 2, nchar(group_prefix)))
  if (!is.na(num)) return(sprintf("%s%02d", prefix_title, num))
  label
}

wb <- createWorkbook()
header_style <- createStyle(fontSize = 11, fontName = "Arial", textDecoration = "bold",
                             halign = "center", valign = "center", fgFill = "#4472C4",
                             fontColour = "#FFFFFF", border = "TopBottomLeftRight", borderColour = "#000000")
text_style <- createStyle(fontSize = 10, fontName = "Arial", halign = "left", valign = "center",
                           border = "TopBottomLeftRight", borderColour = "#D3D3D3")
pvalue_style <- createStyle(fontSize = 10, fontName = "Arial", halign = "right", valign = "center",
                             border = "TopBottomLeftRight", borderColour = "#D3D3D3", numFmt = "0.000")

group_labels <- names(gene_lists)
sheet_summary <- list()

for (label in group_labels) {
  gene_set_label <- format_gene_set_label(label, group_prefix)

  for (ont in go_ontologies) {
    path <- file.path(output_dir, paste0("go_enrichment_", label, "_", ont, ".csv"))
    d <- read_csv_if_nonempty(path)
    if (is.null(d)) next
    d$geneSymbol <- convert_entrez_column_to_symbols(d$geneID, organism_db)
    d$GeneSet <- gene_set_label
    d$Ontology <- ont
    d <- d %>%
      dplyr::select(GeneSet, Ontology, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, geneSymbol) %>%
      dplyr::rename(`Gene Set` = GeneSet, `GO ID` = ID, `GO Term` = Description, `Gene Ratio` = GeneRatio,
                    `Background Ratio` = BgRatio, `P-value` = pvalue, `Adjusted P-value` = p.adjust,
                    `Q-value` = qvalue, `Gene Count` = Count, `Gene Symbols` = geneSymbol) %>%
      dplyr::arrange(`Adjusted P-value`)

    sheet_name <- substr(paste0(label, "_", ont), 1, 31)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, d, headerStyle = header_style)
    addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(d), gridExpand = TRUE)
    addStyle(wb, sheet_name, text_style, rows = 2:(nrow(d) + 1), cols = c(1, 2, 3, 4, 5, 6, 11), gridExpand = TRUE)
    addStyle(wb, sheet_name, pvalue_style, rows = 2:(nrow(d) + 1), cols = c(7, 8, 9), gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = 1:ncol(d), widths = "auto")
    sheet_summary[[sheet_name]] <- nrow(d)
  }

  kegg_path <- file.path(output_dir, paste0("kegg_enrichment_", label, ".csv"))
  kd <- read_csv_if_nonempty(kegg_path)
  if (!is.null(kd)) {
    kd$geneSymbol <- convert_entrez_column_to_symbols(kd$geneID, organism_db)
    kd$GeneSet <- gene_set_label
    kd <- kd %>%
      dplyr::select(GeneSet, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, geneSymbol) %>%
      dplyr::rename(`Gene Set` = GeneSet, `KEGG ID` = ID, `KEGG Pathway` = Description, `Gene Ratio` = GeneRatio,
                    `Background Ratio` = BgRatio, `P-value` = pvalue, `Adjusted P-value` = p.adjust,
                    `Q-value` = qvalue, `Gene Count` = Count, `Gene Symbols` = geneSymbol) %>%
      dplyr::arrange(`Adjusted P-value`)

    sheet_name <- substr(paste0("KEGG_", label), 1, 31)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, kd, headerStyle = header_style)
    addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(kd), gridExpand = TRUE)
    addStyle(wb, sheet_name, text_style, rows = 2:(nrow(kd) + 1), cols = c(1, 2, 3, 4, 5, 10), gridExpand = TRUE)
    addStyle(wb, sheet_name, pvalue_style, rows = 2:(nrow(kd) + 1), cols = c(6, 7, 8), gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = 1:ncol(kd), widths = "auto")
    sheet_summary[[sheet_name]] <- nrow(kd)
  }
}

if (length(sheet_summary) == 0) {
  addWorksheet(wb, "No Results")
  writeData(wb, "No Results", data.frame(Message = "No significant GO or KEGG enrichment results for any group."))
}

# Analysis Info 시트
info_df <- data.frame(
  Parameter = c("Group column", "Groups (excl. all)", "GO ontologies", "Organism Database",
                "Min gene set size", "Max gene set size", "Min gene count", "Sheets with results"),
  Value = c(group_col, paste(group_values, collapse = ", "), paste(go_ontologies, collapse = ", "),
            organism_db_name, min_gs_size, max_gs_size, min_gene_count, length(sheet_summary)),
  stringsAsFactors = FALSE
)
addWorksheet(wb, "Analysis_Info")
writeData(wb, "Analysis_Info", info_df)
addStyle(wb, "Analysis_Info", header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "Analysis_Info", cols = 1:2, widths = c(25, 40))

output_xlsx <- file.path(output_dir, "final_go_results.xlsx")
saveWorkbook(wb, output_xlsx, overwrite = TRUE)
cat(paste("[11_run_group_enrichment] final_go_results.xlsx saved:", output_xlsx,
          "(", length(sheet_summary), "sheets with data)\n"))

# --- 6. go_termcluster_*.csv 취합 -> final_go_clustered_results.xlsx ---
# (05b_generate_clustered_go_table.R과 동일한 컬럼 계약. group_labels 전체(TOTAL 포함)를
# 대상으로 한다는 점만 pairwise용 05b와 다름.)
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

wb_tc <- createWorkbook()
tc_global_counter <- 1L
tc_sheet_summary <- list()

if (!tc_enabled) {
  cat("[11_run_group_enrichment] term_cluster.enabled is not true — writing placeholder final_go_clustered_results.xlsx.\n")
  addWorksheet(wb_tc, "No Results")
  writeData(wb_tc, "No Results", data.frame(Message = "term_cluster is disabled (enrichment.term_cluster.enabled: false)."))
} else {
  for (label in group_labels) {
    gene_set_label <- format_gene_set_label(label, group_prefix)
    for (ont in go_ontologies) {
      in_csv <- file.path(output_dir, paste0("go_termcluster_", label, "_", ont, ".csv"))
      if (!file.exists(in_csv)) next
      d <- read.csv(in_csv, stringsAsFactors = FALSE)
      if (nrow(d) == 0) next

      relabeled <- relabel_clusters(d$cluster, tc_global_counter)
      tc_global_counter <- relabeled$next_counter
      d$geneSymbol <- convert_entrez_column_to_symbols(d$geneID, organism_db)

      d_out <- data.frame(
        `Gene Set`         = gene_set_label,
        `Ontology`         = ont,
        `GO ID`            = d$ID,
        `GO Term`          = d$Description,
        `Gene Ratio`       = d$GeneRatio,
        `Background Ratio` = d$BgRatio,
        `P-value`          = d$pvalue,
        `Adjusted P-value` = d$p.adjust,
        `Gene Count`       = d$Count,
        `Gene Symbols`     = d$geneSymbol,
        `cluster_id`       = relabeled$ids,
        check.names = FALSE, stringsAsFactors = FALSE
      )
      d_out <- d_out[order(d_out$cluster_id == "Singleton", d_out$cluster_id, d_out$`Adjusted P-value`), ]

      sheet_name <- substr(paste0(label, "_", ont), 1, 31)
      addWorksheet(wb_tc, sheet_name)
      writeData(wb_tc, sheet_name, d_out, headerStyle = header_style)
      addStyle(wb_tc, sheet_name, header_style, rows = 1, cols = 1:ncol(d_out), gridExpand = TRUE)
      setColWidths(wb_tc, sheet_name, cols = 1:ncol(d_out), widths = "auto")
      tc_sheet_summary[[sheet_name]] <- nrow(d_out)
    }
  }

  if (length(tc_sheet_summary) == 0) {
    addWorksheet(wb_tc, "No Results")
    writeData(wb_tc, "No Results", data.frame(Message = "No significant GO terms available to cluster for any group."))
  }
}

out_tc_xlsx <- file.path(output_dir, "final_go_clustered_results.xlsx")
saveWorkbook(wb_tc, out_tc_xlsx, overwrite = TRUE)
cat(paste("[11_run_group_enrichment] final_go_clustered_results.xlsx saved:", out_tc_xlsx,
          "(", length(tc_sheet_summary), "sheets,", tc_global_counter - 1L, "clusters total )\n"))

# --- 7. go_rrvgo_*.csv 취합 -> final_go_rrvgo_clustered_results.xlsx ---
# (05d_generate_rrvgo_clustered_go_table.R과 동일한 컬럼 계약. TOTAL 포함.)
wb_rr <- createWorkbook()
rr_global_counter <- 1L
rr_sheet_summary <- list()

if (!rr_enabled) {
  cat("[11_run_group_enrichment] rrvgo.enabled is not true — writing placeholder final_go_rrvgo_clustered_results.xlsx.\n")
  addWorksheet(wb_rr, "No Results")
  writeData(wb_rr, "No Results", data.frame(Message = "rrvgo is disabled (enrichment.rrvgo.enabled: false)."))
} else {
  for (label in group_labels) {
    gene_set_label <- format_gene_set_label(label, group_prefix)
    for (ont in go_ontologies) {
      in_csv <- file.path(output_dir, paste0("go_rrvgo_", label, "_", ont, ".csv"))
      if (!file.exists(in_csv)) next
      d <- read.csv(in_csv, stringsAsFactors = FALSE)
      if (nrow(d) == 0) next

      enrich_csv <- file.path(output_dir, paste0("go_enrichment_", label, "_", ont, ".csv"))
      if (!file.exists(enrich_csv)) next
      e <- read.csv(enrich_csv, stringsAsFactors = FALSE)
      d <- merge(d, e[, c("ID", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "geneID")],
                 by.x = "go", by.y = "ID", all.x = TRUE)
      d <- d[!is.na(d$Count), ]
      if (nrow(d) == 0) next

      relabeled <- relabel_clusters(d$cluster, rr_global_counter)
      rr_global_counter <- relabeled$next_counter
      d$geneSymbol <- convert_entrez_column_to_symbols(d$geneID, organism_db)

      d_out <- data.frame(
        `Gene Set`            = gene_set_label,
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
        check.names = FALSE, stringsAsFactors = FALSE
      )
      d_out <- d_out[order(d_out$cluster_id == "Singleton", d_out$cluster_id, d_out$`Adjusted P-value`), ]

      sheet_name <- substr(paste0(label, "_", ont), 1, 31)
      addWorksheet(wb_rr, sheet_name)
      writeData(wb_rr, sheet_name, d_out, headerStyle = header_style)
      addStyle(wb_rr, sheet_name, header_style, rows = 1, cols = 1:ncol(d_out), gridExpand = TRUE)
      setColWidths(wb_rr, sheet_name, cols = 1:ncol(d_out), widths = "auto")
      rr_sheet_summary[[sheet_name]] <- nrow(d_out)
    }
  }

  if (length(rr_sheet_summary) == 0) {
    addWorksheet(wb_rr, "No Results")
    writeData(wb_rr, "No Results", data.frame(Message = "No significant GO terms available to cluster for any group."))
  }
}

out_rr_xlsx <- file.path(output_dir, "final_go_rrvgo_clustered_results.xlsx")
saveWorkbook(wb_rr, out_rr_xlsx, overwrite = TRUE)
cat(paste("[11_run_group_enrichment] final_go_rrvgo_clustered_results.xlsx saved:", out_rr_xlsx,
          "(", length(rr_sheet_summary), "sheets,", rr_global_counter - 1L, "clusters total )\n"))

cat("[11_run_group_enrichment] Done.\n")
