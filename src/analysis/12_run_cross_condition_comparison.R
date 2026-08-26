# 파일 경로: src/analysis/12_run_cross_condition_comparison.R
#
# Cross-Condition GO 비교 (docs/plan_cross_condition_go.md 참고).
# 한 프로젝트 안에 pairwise 비교("조건")가 여러 개 있을 때(dose-response,
# time-course 등) 03_enrichment_analysis.R이 각 비교마다 독립적으로 만든
# go_enrichment_{up|down}_{ont}.csv를 서로 대조해서:
#   A. 여러 조건에서 공통으로 UP/DOWN인 term (find_common_direction)
#   B. 그룹 간 방향이 역전되는 term + gene Jaccard 검증 (find_direction_flip)
#   C. 특정 조건 조합에만(또는 없이) 나타나는 term (find_exclusive)
#   D. 같은 그룹 내에서 방향이 혼재하는 term (find_mixed)
# 을 찾는다. 03/11의 term_cluster·rrvgo(한 비교 안에서 term을 묶는 것)와는 다른
# 축으로, "이미 각각 실행된 독립적인 비교들" 사이를 대조하는 분석이다.
#
# 사용법: Rscript 12_run_cross_condition_comparison.R <config_path> <output_dir>
#   예) Rscript 12_run_cross_condition_comparison.R configs/config_X.yml \
#         output/X/cross_condition

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(AnnotationDbi)
  library(ggplot2)
  library(ggupset)
  library(pheatmap)
  library(dplyr)
  library(openxlsx)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 & config 로드 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 12_run_cross_condition_comparison.R <config_path> <output_dir>")
}
config_path <- args[1]
output_dir  <- args[2]

config <- yaml.load_file(config_path)
cc_cfg <- config$enrichment$cross_condition

if (is.null(cc_cfg) || !isTRUE(cc_cfg$enabled)) {
  cat("[12_run_cross_condition_comparison] cross_condition.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

project_dir <- config$output_dir
go_ontologies <- cc_cfg$go_ontologies %||% c("BP")
fdr_cutoff <- cc_cfg$fdr_cutoff %||% 0.05
fe_cutoff  <- cc_cfg$fold_enrichment_cutoff %||% 2.0
sparsity_warn_threshold <- cc_cfg$sparsity_warn_threshold %||% 20
jaccard_warn <- cc_cfg$gene_overlap_jaccard_warn %||% 0.3

conditions <- cc_cfg$conditions
if (is.null(conditions)) {
  pairs_cfg <- config$de_analysis$pairwise_comparisons %||% list()
  conditions <- vapply(pairs_cfg, function(p) paste0(p[[1]], "_vs_", p[[2]]), character(1))
}
if (length(conditions) < 2) {
  stop("[FATAL] cross_condition requires at least 2 conditions (pairwise comparisons).")
}
cat(paste("[12_run_cross_condition_comparison] Conditions (", length(conditions), "):",
          paste(conditions, collapse = ", "), "\n"))

min_conditions_common <- cc_cfg$min_conditions_common %||% round(0.75 * length(conditions))
cat(paste("[12_run_cross_condition_comparison] min_conditions_common =", min_conditions_common, "\n"))

viz_cfg <- cc_cfg$visualization %||% list()
viz_enabled <- if (is.null(viz_cfg$enabled)) TRUE else isTRUE(viz_cfg$enabled)
dotplot_max_terms <- viz_cfg$dotplot_max_terms %||% 60
up_color   <- config$plot_aesthetics$volcano$up_color   %||% "#FF5733"
down_color <- config$plot_aesthetics$volcano$down_color %||% "#3375FF"
# 조건별 signed -log10(FDR) 패턴이 비슷한 GO term끼리 축 위에서 인접하도록 묶어준다.
# heatmap은 pheatmap이 직접 클러스터링+덴드로그램을 그리고, dot plot은 ggplot2가
# 덴드로그램을 못 그리므로 여기서 미리 계산한 순서로 행(y축)만 재배열한다.
cluster_terms_enabled <- if (is.null(viz_cfg$cluster_terms)) TRUE else isTRUE(viz_cfg$cluster_terms)
cluster_method <- viz_cfg$cluster_method %||% "average"

# common_strict 결과(비교 "결과" term_id 목록)에 rrvgo를 재적용 — 03/11의 rrvgo와
# 같은 파라미터를 기본으로 공유하되(enrichment.rrvgo), cross_condition에서만 별도로
# 끄고 싶을 때를 위해 semantic_filtering.enabled로 개별 토글 가능.
semantic_cfg <- cc_cfg$semantic_filtering %||% list()
semantic_enabled <- if (is.null(semantic_cfg$enabled)) TRUE else isTRUE(semantic_cfg$enabled)
rrvgo_global_cfg <- config$enrichment$rrvgo %||% list()
rr_method    <- semantic_cfg$method    %||% rrvgo_global_cfg$method    %||% "Rel"
rr_threshold <- semantic_cfg$threshold %||% rrvgo_global_cfg$threshold %||% 0.7

species_info <- config$databases[[config$species]]
organism_db_name <- species_info$organism_db
if (!require(organism_db_name, character.only = TRUE, quietly = TRUE)) {
  stop(paste("[FATAL] Required organism DB package", organism_db_name, "is not installed."))
}
organism_db <- get(organism_db_name)

convert_entrez_to_symbols <- function(entrez_ids) {
  if (length(entrez_ids) == 0) return(character(0))
  lookup <- tryCatch({
    mapIds(organism_db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  }, error = function(e) setNames(entrez_ids, entrez_ids))
  out <- unname(lookup[entrez_ids])
  out[is.na(out)] <- entrez_ids[is.na(out)]
  out
}

parse_ratio <- function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1] / parts[2]
}
compute_fold_enrichment <- function(go_df) {
  mapply(function(gr, br) parse_ratio(gr) / parse_ratio(br), go_df$GeneRatio, go_df$BgRatio)
}

# --- 2. 조건별 go_enrichment_{up|down}_{ont}.csv 로드 & 유의 term만 남기기 ---
load_condition_direction <- function(condition, direction, ont) {
  # KEGG는 03_enrichment_analysis.R에서 ontology 구분 없이 kegg_enrichment_{direction}.csv
  # 하나로 저장되므로(BP/CC/MF 같은 하위분류가 없음) GO와 파일명 패턴이 다르다.
  path <- if (toupper(ont) == "KEGG") {
    file.path(project_dir, "pairwise", condition, "enrichment", paste0("kegg_enrichment_", direction, ".csv"))
  } else {
    file.path(project_dir, "pairwise", condition, "enrichment", paste0("go_enrichment_", direction, "_", ont, ".csv"))
  }
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d$FoldEnrichment <- compute_fold_enrichment(d)
  d <- d[!is.na(d$p.adjust) & d$p.adjust < fdr_cutoff &
         !is.na(d$FoldEnrichment) & d$FoldEnrichment > fe_cutoff, ]
  if (nrow(d) == 0) return(NULL)
  d$condition <- condition
  d$direction <- toupper(direction)
  d$GeneRatioNum <- sapply(d$GeneRatio, parse_ratio)
  d[, c("ID", "Description", "p.adjust", "FoldEnrichment", "GeneRatioNum", "Count", "geneID", "condition", "direction")]
}

condition_count_log <- character()

run_for_ontology <- function(ont) {
  cat(sprintf("\n[12_run_cross_condition_comparison] === Ontology: %s ===\n", ont))

  parts <- list()
  for (cond in conditions) {
    for (dir in c("up", "down")) {
      d <- load_condition_direction(cond, dir, ont)
      n <- if (is.null(d)) 0L else nrow(d)
      msg <- sprintf("  %s / %s(%s): %d significant terms", cond, dir, ont, n)
      cat(msg, "\n")
      condition_count_log[[length(condition_count_log) + 1]] <<- msg
      if (n < sparsity_warn_threshold) {
        warn_msg <- sprintf("  [WARN] %s/%s(%s) has only %d significant terms (< %d) — common-term results may be skewed by this sparsity.",
                             cond, dir, ont, n, sparsity_warn_threshold)
        cat(warn_msg, "\n")
        condition_count_log[[length(condition_count_log) + 1]] <<- warn_msg
      }
      if (!is.null(d)) parts[[paste(cond, dir)]] <- d
    }
  }
  if (length(parts) == 0) {
    cat(sprintf("  [12_run_cross_condition_comparison] No significant terms found for any condition/direction in %s — skipping.\n", ont))
    return(invisible(NULL))
  }
  all_df <- do.call(rbind, parts)
  rownames(all_df) <- NULL

  term_desc <- function(id) {
    rows <- all_df[all_df$ID == id, ]
    rows$Description[which.min(rows$p.adjust)]
  }
  gene_set_for <- function(sub_df) {
    if (nrow(sub_df) == 0) return(character(0))
    unique(unlist(strsplit(sub_df$geneID, "/", fixed = TRUE)))
  }

  # --- A. find_common_direction() ---
  common_strict_ids <- list(); common_loose_ids <- list()
  for (dir in c("UP", "DOWN")) {
    opp <- setdiff(c("UP", "DOWN"), dir)
    dir_df <- all_df[all_df$direction == dir, ]
    opp_ids <- unique(all_df$ID[all_df$direction == opp])
    if (nrow(dir_df) == 0) next

    agg <- do.call(rbind, lapply(split(dir_df, dir_df$ID), function(sub) {
      data.frame(
        `GO ID` = sub$ID[1],
        `GO Term` = term_desc(sub$ID[1]),
        `N Conditions` = length(unique(sub$condition)),
        `Conditions` = paste(sort(unique(sub$condition)), collapse = "; "),
        `Min Adjusted P-value` = min(sub$p.adjust),
        `Gene Symbols (union)` = paste(convert_entrez_to_symbols(gene_set_for(sub)), collapse = "/"),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }))
    agg <- agg[order(-agg$`N Conditions`, agg$`Min Adjusted P-value`), ]

    loose <- agg[agg$`N Conditions` >= min_conditions_common, ]
    strict <- loose[!(loose$`GO ID` %in% opp_ids), ]

    write.csv(loose,  file.path(output_dir, sprintf("common_%s_loose_%s.csv", tolower(dir), ont)), row.names = FALSE)
    write.csv(strict, file.path(output_dir, sprintf("common_%s_strict_%s.csv", tolower(dir), ont)), row.names = FALSE)
    cat(sprintf("  [common_%s] loose=%d terms, strict=%d terms (>= %d conditions)\n",
                dir, nrow(loose), nrow(strict), min_conditions_common))
    common_strict_ids[[dir]] <- strict$`GO ID`
    common_loose_ids[[dir]]  <- loose$`GO ID`
  }

  groups <- cc_cfg$groups %||% list()
  all_flip_ids <- character(0)
  all_excl_ids <- character(0)
  all_mixed_ids <- character(0)
  term_ids_for <- function(conds, dir) unique(all_df$ID[all_df$condition %in% conds & all_df$direction == dir])

  if (length(groups) == 0) {
    cat("  [12_run_cross_condition_comparison] No groups configured — skipping flip/exclusive/mixed analyses.\n")
  } else {
  gene_overlap_rows <- list()

  # --- B. find_direction_flip() ---
  flips <- cc_cfg$flips %||% list()
  for (flip in flips) {
    group_a_name <- flip[[1]]; group_b_name <- flip[[2]]; dir_from <- toupper(flip[[3]])
    dir_to <- setdiff(c("UP", "DOWN"), dir_from)
    group_a <- groups[[group_a_name]]; group_b <- groups[[group_b_name]]
    if (is.null(group_a) || is.null(group_b)) {
      cat(sprintf("  [flip] group '%s' or '%s' not found in groups — skipped.\n", group_a_name, group_b_name))
      next
    }

    flip_ids <- setdiff(
      intersect(term_ids_for(group_a, dir_from), term_ids_for(group_b, dir_to)),
      union(term_ids_for(group_a, dir_to), term_ids_for(group_b, dir_from))
    )

    label <- sprintf("flip_%s_%s_to_%s_%s", dir_from, group_a_name, dir_to, group_b_name)
    if (length(flip_ids) == 0) {
      cat(sprintf("  [%s] 0 terms.\n", label))
      next
    }
    all_flip_ids <- c(all_flip_ids, flip_ids)

    rows <- lapply(flip_ids, function(id) {
      sub_a <- all_df[all_df$ID == id & all_df$condition %in% group_a & all_df$direction == dir_from, ]
      sub_b <- all_df[all_df$ID == id & all_df$condition %in% group_b & all_df$direction == dir_to, ]
      genes_a <- gene_set_for(sub_a); genes_b <- gene_set_for(sub_b)
      jacc <- if (length(union(genes_a, genes_b)) == 0) NA_real_ else length(intersect(genes_a, genes_b)) / length(union(genes_a, genes_b))
      data.frame(
        `GO ID` = id,
        `GO Term` = term_desc(id),
        `Group A (from)` = group_a_name, `Direction A` = dir_from,
        `Conditions A` = paste(sort(unique(sub_a$condition)), collapse = "; "),
        `Min Adj P A` = min(sub_a$p.adjust),
        `Group B (to)` = group_b_name, `Direction B` = dir_to,
        `Conditions B` = paste(sort(unique(sub_b$condition)), collapse = "; "),
        `Min Adj P B` = min(sub_b$p.adjust),
        `Gene Jaccard (A vs B)` = round(jacc, 3),
        `Low Overlap` = !is.na(jacc) & jacc < jaccard_warn,
        `Gene Symbols A` = paste(convert_entrez_to_symbols(genes_a), collapse = "/"),
        `Gene Symbols B` = paste(convert_entrez_to_symbols(genes_b), collapse = "/"),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    flip_df <- do.call(rbind, rows)
    flip_df <- flip_df[order(flip_df$`Min Adj P A`), ]
    write.csv(flip_df, file.path(output_dir, paste0(label, "_", ont, ".csv")), row.names = FALSE)
    cat(sprintf("  [%s] %d terms (%d low-overlap flagged, Jaccard < %.1f)\n",
                label, nrow(flip_df), sum(flip_df$`Low Overlap`), jaccard_warn))
    gene_overlap_rows[[label]] <- data.frame(
      Comparison = label, `GO ID` = flip_df$`GO ID`, `GO Term` = flip_df$`GO Term`,
      Jaccard = flip_df$`Gene Jaccard (A vs B)`, `Low Overlap` = flip_df$`Low Overlap`,
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }

  if (length(gene_overlap_rows) > 0) {
    write.csv(do.call(rbind, gene_overlap_rows),
              file.path(output_dir, paste0("gene_overlap_summary_", ont, ".csv")), row.names = FALSE)
  }

  # --- C. find_exclusive() ---
  exclusives <- cc_cfg$exclusives %||% list()
  for (spec in exclusives) {
    label <- spec$label %||% paste(spec$target, collapse = "-")
    target <- spec$target; absent <- spec$absent; direction <- spec$direction

    if (is.null(direction)) {
      target_sets <- lapply(target, function(c) unique(all_df$ID[all_df$condition == c]))
    } else {
      target_sets <- lapply(target, function(c) unique(all_df$ID[all_df$condition == c & all_df$direction == toupper(direction)]))
    }
    target_ids <- Reduce(intersect, target_sets)
    absent_ids <- unique(all_df$ID[all_df$condition %in% absent])
    excl_ids <- setdiff(target_ids, absent_ids)

    if (length(excl_ids) == 0) {
      cat(sprintf("  [exclusive:%s] 0 terms.\n", label))
      next
    }
    all_excl_ids <- c(all_excl_ids, excl_ids)
    rows <- lapply(excl_ids, function(id) {
      sub <- all_df[all_df$ID == id & all_df$condition %in% target, ]
      data.frame(
        `GO ID` = id, `GO Term` = term_desc(id),
        `Target Conditions` = paste(sort(unique(sub$condition)), collapse = "; "),
        `Direction` = paste(sort(unique(sub$direction)), collapse = "/"),
        `Min Adjusted P-value` = min(sub$p.adjust),
        `Gene Symbols` = paste(convert_entrez_to_symbols(gene_set_for(sub)), collapse = "/"),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    excl_df <- do.call(rbind, rows)
    excl_df <- excl_df[order(excl_df$`Min Adjusted P-value`), ]
    write.csv(excl_df, file.path(output_dir, paste0("exclusive_", label, "_", ont, ".csv")), row.names = FALSE)
    cat(sprintf("  [exclusive:%s] %d terms.\n", label, nrow(excl_df)))
  }

  # --- D. find_mixed() ---
  for (group_name in names(groups)) {
    group_conds <- groups[[group_name]]
    up_ids <- term_ids_for(group_conds, "UP")
    down_ids <- term_ids_for(group_conds, "DOWN")
    mixed_ids <- intersect(up_ids, down_ids)
    if (length(mixed_ids) == 0) {
      cat(sprintf("  [mixed:%s] 0 terms.\n", group_name))
      next
    }
    all_mixed_ids <- c(all_mixed_ids, mixed_ids)
    rows <- lapply(mixed_ids, function(id) {
      sub <- all_df[all_df$ID == id & all_df$condition %in% group_conds, ]
      data.frame(
        `GO ID` = id, `GO Term` = term_desc(id),
        `UP in` = paste(sort(unique(sub$condition[sub$direction == "UP"])), collapse = "; "),
        `DOWN in` = paste(sort(unique(sub$condition[sub$direction == "DOWN"])), collapse = "; "),
        `Min Adjusted P-value` = min(sub$p.adjust),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    mixed_df <- do.call(rbind, rows)
    mixed_df <- mixed_df[order(mixed_df$`Min Adjusted P-value`), ]
    write.csv(mixed_df, file.path(output_dir, paste0("mixed_", group_name, "_", ont, ".csv")), row.names = FALSE)
    cat(sprintf("  [mixed:%s] %d terms.\n", group_name, nrow(mixed_df)))
  }
  }

  # --- Phase 4: rrvgo 의미론적 축약 재적용 (GO만 해당, KEGG는 term이 적어 생략) ---
  # 03/11이 "한 비교의 raw enrichment 결과"에 rrvgo를 적용하는 것과 달리, 여기서는
  # cross-condition 비교 "결과"(common_strict) term_id 목록에 재적용한다 — 조건마다
  # 서로 다른 representative term이 뽑히는 문제(§2)를 피하기 위해 원본 term_id로 먼저
  # 비교를 끝낸 뒤 마지막에 한 번만 축약한다.
  semantic_representative_ids <- character(0)
  if (isTRUE(semantic_enabled) && toupper(ont) != "KEGG") {
    for (dir in c("UP", "DOWN")) {
      strict_ids <- common_strict_ids[[dir]]
      if (length(strict_ids) < 2) next
      scores <- vapply(strict_ids, function(id) -log10(max(min(all_df$p.adjust[all_df$ID == id]), 1e-300)), numeric(1))
      names(scores) <- strict_ids
      reduced <- tryCatch({
        suppressPackageStartupMessages(library(rrvgo))
        simMatrix <- calculateSimMatrix(strict_ids, orgdb = organism_db_name, ont = ont, method = rr_method)
        reduceSimMatrix(simMatrix, scores = scores, threshold = rr_threshold, orgdb = organism_db_name)
      }, error = function(e) {
        cat(sprintf("  [semantic_filter] common_%s_strict(%s) failed: %s\n", tolower(dir), ont, e$message))
        NULL
      })
      if (is.null(reduced)) next
      out_csv <- file.path(output_dir, sprintf("common_%s_strict_rrvgo_%s.csv", tolower(dir), ont))
      write.csv(reduced[order(reduced$cluster, -reduced$score), ], out_csv, row.names = FALSE)
      # 클러스터(=parent)당 대표 term(가장 score가 높은 term)만 최종 시각화 후보로 사용
      rep_ids <- unlist(lapply(split(reduced, reduced$parent), function(g) g$go[which.max(g$score)]))
      semantic_representative_ids <- c(semantic_representative_ids, unname(rep_ids))
      cat(sprintf("  [semantic_filter] common_%s_strict(%s): %d terms -> %d representative parent groups\n",
                  tolower(dir), ont, length(strict_ids), length(unique(reduced$parent))))
    }
    semantic_representative_ids <- unique(semantic_representative_ids)
  }

  # --- Phase 3: 시각화 ---
  # groups 설정 여부와 무관하게(공통 term + 전체 유의 term 구조는 groups 없이도 계산됨)
  # 항상 시도한다. 개별 플롯이 실패해도 나머지 분석/출력은 영향받지 않도록 tryCatch로 감싼다.
  if (viz_enabled) {
    # 후보 term 집합 기준으로 dot plot(term 수 <= dotplot_max_terms) 또는 heatmap(초과 시)을
    # 그린다. out_name_dot/out_name_heatmap/title은 "raw"(원본 후보)와 "semantic"(rrvgo
    # 축약 후 대표 term만) 두 버전에서 재사용한다.
    # id -> 조건별 signed -log10(FDR) 벡터. dot plot 클러스터링 순서 계산과 heatmap 행렬
    # 구성 둘 다에 쓰는 공용 헬퍼(동일한 값 기준으로 묶어야 두 플롯의 "패턴"이 일치한다).
    signed_matrix_for <- function(ids) {
      mat <- matrix(0, nrow = length(ids), ncol = length(conditions), dimnames = list(ids, conditions))
      for (i in seq_along(ids)) {
        id <- ids[i]
        for (cond in conditions) {
          sub <- all_df[all_df$ID == id & all_df$condition == cond, ]
          if (nrow(sub) == 0) next
          best <- sub[which.min(sub$p.adjust), ]
          mat[i, cond] <- ifelse(best$direction == "UP", 1, -1) * -log10(best$p.adjust)
        }
      }
      mat
    }

    trunc_term <- function(x) ifelse(nchar(x) > 55, paste0(substr(x, 1, 52), "..."), x)

    # id 목록(주어진 순서 그대로) -> long-format 표(GO ID/Term/Condition/Direction/GeneRatio/
    # Adjusted P-value/signed -log10(FDR)/방향역전 여부). dot plot·heatmap이 그리는 값과
    # 정확히 같은 표를 그림 재현/검증용 source data로도 그대로 저장한다.
    build_long_df <- function(ids) {
      rows <- list()
      for (id in ids) {
        desc <- term_desc(id)
        for (cond in conditions) {
          sub <- all_df[all_df$ID == id & all_df$condition == cond, ]
          if (nrow(sub) == 0) next
          best <- sub[which.min(sub$p.adjust), ]
          rows[[length(rows) + 1]] <- data.frame(
            `GO ID` = id, `GO Term` = desc, Condition = cond, Direction = best$direction,
            `Gene Ratio` = best$GeneRatioNum, `Adjusted P-value` = best$p.adjust,
            `Signed -log10(FDR)` = ifelse(best$direction == "UP", 1, -1) * -log10(best$p.adjust),
            `Direction Flip` = id %in% all_flip_ids,
            check.names = FALSE, stringsAsFactors = FALSE
          )
        }
      }
      do.call(rbind, rows)
    }

    draw_heatmap <- function(ids, out_name_heatmap, title) {
      if (length(ids) < 2) return(invisible(NULL))
      cat(sprintf("  [viz] Heatmap: %d terms -> %s\n", length(ids), out_name_heatmap))
      tryCatch({
        raw_labels <- vapply(ids, term_desc, character(1))
        trunc_labels <- ifelse(nchar(raw_labels) > 55, paste0(substr(raw_labels, 1, 52), "..."), raw_labels)
        term_labels <- make.unique(trunc_labels)
        mat <- signed_matrix_for(ids)
        rownames(mat) <- term_labels
        # 연속값(signed -log10 FDR) 기준 거리 — 색은 아래에서 이산화하지만 클러스터링/
        # 덴드로그램은 원래의 연속적 유사도를 그대로 반영해야 하므로 따로 보관해둔다.
        dist_rows <- dist(mat)

        # FDR은 heavy-tailed라(예: 1e-10와 1e-62가 같은 매트릭스에 공존) 연속 컬러 스케일을
        # 쓰면 극단값 하나가 나머지 셀 전부를 흐리게 눌러버린다 — 그래서 유의성 구간
        # (0.05/0.01/0.001)으로 이산화해서 칠한다. 0은 "해당 조건에서 유의하지 않음(absent)".
        tier_code <- function(v) {
          if (v == 0) return(4L)
          if (v > 0) { if (v >= 3) return(7L); if (v >= 2) return(6L); return(5L) }
          if (v <= -3) return(1L); if (v <= -2) return(2L); return(3L)
        }
        mat_tier <- matrix(vapply(as.vector(mat), tier_code, integer(1)),
                            nrow = nrow(mat), dimnames = dimnames(mat))

        down_shades <- colorRampPalette(c(down_color, "white"))(4)[1:3]
        up_shades   <- colorRampPalette(c("white", up_color))(4)[2:4]
        tier_colors <- c(down_shades, "white", up_shades)
        tier_labels <- c("DOWN FDR<0.001", "DOWN FDR<0.01", "DOWN FDR<0.05", "absent",
                          "UP FDR<0.05", "UP FDR<0.01", "UP FDR<0.001")

        n_rows <- length(ids)
        row_fontsize <- if (n_rows > 100) 5 else if (n_rows > 60) 6 else 8
        # cluster_rows=TRUE면 pheatmap이 dist_rows(연속값 거리)로 재클러스터링해서
        # 덴드로그램을 그린다 — 표시되는 색은 이산화됐지만 묶는 기준은 연속값 그대로.
        img_width <- max(1400, 480 + max(nchar(term_labels)) * 8 + 130 * length(conditions) + 60)
        img_height <- max(900, 20 * n_rows + 250)
        png(file.path(output_dir, out_name_heatmap), width = img_width, height = img_height, res = 150)
        pheatmap(mat_tier, color = tier_colors, breaks = seq(0.5, 7.5, by = 1),
                 cluster_cols = FALSE, cluster_rows = cluster_terms_enabled,
                 clustering_distance_rows = dist_rows,
                 legend_breaks = 1:7, legend_labels = tier_labels,
                 fontsize_row = row_fontsize, fontsize_col = 9, main = title)
        dev.off()
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        cat(paste("  [viz] heatmap failed:", e$message, "\n"))
      })
    }

    draw_dotplot <- function(ids, out_name_dot, title) {
      if (length(ids) == 0) return(invisible(NULL))
      cat(sprintf("  [viz] Dot plot: %d terms -> %s\n", length(ids), out_name_dot))
      tryCatch({
        plot_df <- build_long_df(ids)
        plot_df$`GO Term` <- trunc_term(plot_df$`GO Term`)
        term_order <- trunc_term(vapply(ids, term_desc, character(1)))
        plot_df$`GO Term` <- factor(plot_df$`GO Term`, levels = rev(unique(term_order)))

        p <- ggplot(plot_df, aes(x = Condition, y = `GO Term`)) +
          geom_point(aes(size = `Gene Ratio`, color = `Signed -log10(FDR)`)) +
          scale_color_gradient2(low = down_color, mid = "grey85", high = up_color, midpoint = 0,
                                name = "sign(dir) x -log10(FDR)") +
          scale_size_continuous(name = "Gene Ratio", range = c(1, 8)) +
          theme_bw(base_size = 11) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = NULL, y = NULL, title = title,
               caption = if (any(plot_df$`Direction Flip`)) "* black outline = direction-flip term" else NULL)
        if (any(plot_df$`Direction Flip`)) {
          p <- p + geom_point(data = plot_df[plot_df$`Direction Flip`, ], shape = 21, size = 3.2,
                               color = "black", stroke = 1, show.legend = FALSE)
        }
        # 폭 = legend(~2.6in) + 조건당 패널 폭(~0.8in) + y축 라벨 폭(글자 수 비례, 최대 55자 기준)
        label_w <- 0.085 * max(nchar(as.character(levels(plot_df$`GO Term`))))
        w <- 2.6 + 0.8 * length(conditions) + label_w
        h <- max(5, 0.28 * length(ids) + 2)
        ggsave(file.path(output_dir, out_name_dot), plot = p, width = w, height = h, bg = "white", limitsize = FALSE)
      }, error = function(e) cat(paste("  [viz] dot plot failed:", e$message, "\n")))
    }

    # heatmap은 항상 전체 후보로, dot plot은 유의성 상위 dotplot_max_terms개만(별도로
    # 그 안에서 다시 클러스터링) — 둘 다 항상 같이 생성한다. source data CSV는 heatmap이
    # 쓰는 전체 후보 기준 long-format 표(= dot plot 표의 상위集合)로 하나만 저장한다.
    render_plot <- function(ids, out_name_dot, out_name_heatmap, out_name_data, title) {
      if (length(ids) == 0) return(invisible(NULL))
      term_rank <- vapply(ids, function(id) min(all_df$p.adjust[all_df$ID == id]), numeric(1))
      ids_by_sig <- ids[order(term_rank)]

      ids_full <- ids_by_sig
      if (cluster_terms_enabled && length(ids_full) >= 3) {
        hc <- tryCatch(hclust(dist(signed_matrix_for(ids_full)), method = cluster_method), error = function(e) NULL)
        if (!is.null(hc)) ids_full <- ids_full[hc$order]
      }

      write.csv(build_long_df(ids_full), file.path(output_dir, out_name_data), row.names = FALSE)

      draw_heatmap(ids_full, out_name_heatmap, sub("Dot Plot", "Heatmap", title, fixed = TRUE))

      # dot plot 상위 term은 순수 FDR로만 자르지 않고 방향별로 균형있게 뽑는다 — DOWN이
      # term 수/통계적 강도 모두 UP을 압도하는 경우(실측 사례 있음) 상위 N개를 그냥
      # FDR로만 자르면 UP이 하나도 안 남는 문제가 있었다. 각 term의 "대표 방향"은 그
      # term의 전체 조건 중 가장 작은 p.adjust를 낸 방향으로 정한다(term_rank와 동일 기준).
      best_dir <- vapply(ids_by_sig, function(id) {
        sub <- all_df[all_df$ID == id, ]
        sub$direction[which.min(sub$p.adjust)]
      }, character(1))
      half <- ceiling(dotplot_max_terms / 2)
      up_top   <- head(ids_by_sig[best_dir == "UP"],   half)
      down_top <- head(ids_by_sig[best_dir == "DOWN"], half)
      dot_ids <- ids_by_sig[ids_by_sig %in% c(up_top, down_top)]  # ids_by_sig 순서(유의성 순) 유지

      if (cluster_terms_enabled && length(dot_ids) >= 3) {
        hc2 <- tryCatch(hclust(dist(signed_matrix_for(dot_ids)), method = cluster_method), error = function(e) NULL)
        if (!is.null(hc2)) dot_ids <- dot_ids[hc2$order]
      }
      dot_title <- if (length(ids_by_sig) > length(dot_ids)) {
        sprintf("%s [top %d UP + %d DOWN of %d by FDR]", title, length(up_top), length(down_top), length(ids_by_sig))
      } else title
      draw_dotplot(dot_ids, out_name_dot, dot_title)
    }

    candidate_ids <- unique(c(common_strict_ids[["UP"]], common_strict_ids[["DOWN"]], all_flip_ids, all_excl_ids))
    if (length(candidate_ids) == 0) candidate_ids <- unique(c(common_loose_ids[["UP"]], common_loose_ids[["DOWN"]]))
    if (length(candidate_ids) == 0) {
      by_p <- all_df[order(all_df$p.adjust), ]
      candidate_ids <- unique(head(by_p$ID, 20))
    }
    render_plot(candidate_ids,
                sprintf("cross_condition_dotplot_%s.png", ont), sprintf("cross_condition_heatmap_%s.png", ont),
                sprintf("cross_condition_plot_data_%s.csv", ont),
                sprintf("Cross-Condition GO Dot Plot (%s)", ont))

    if (length(semantic_representative_ids) > 0) {
      render_plot(semantic_representative_ids,
                  sprintf("cross_condition_dotplot_semantic_%s.png", ont), sprintf("cross_condition_heatmap_semantic_%s.png", ont),
                  sprintf("cross_condition_plot_data_semantic_%s.csv", ont),
                  sprintf("Cross-Condition GO Dot Plot - Semantic Filtered (%s)", ont))
    }

    # 같은 그룹 안에서 방향이 혼재하는 term(D. find_mixed) 전용 패널. common_strict는
    # 반대 방향이 하나라도 있으면 제외하는 로직이라 이 term들은 원본/semantic 패널에는
    # 거의 안 나타난다 — 그래서 별도 후보 집합으로 항상 같은 render_plot()을 재사용해
    # heatmap(전체) + dot plot(상위 term) + source data를 만든다.
    mixed_ids_unique <- unique(all_mixed_ids)
    if (length(mixed_ids_unique) > 0) {
      render_plot(mixed_ids_unique,
                  sprintf("cross_condition_dotplot_mixed_%s.png", ont), sprintf("cross_condition_heatmap_mixed_%s.png", ont),
                  sprintf("cross_condition_plot_data_mixed_%s.csv", ont),
                  sprintf("Cross-Condition GO Dot Plot - Mixed Direction (%s)", ont))
    }

    # --- UpSet plot: 조건별 유의 term 중첩 구조 (UP/DOWN 각각) ---
    for (dir in c("UP", "DOWN")) {
      sub <- all_df[all_df$direction == dir, ]
      ids <- unique(sub$ID)
      if (length(ids) == 0) next
      membership <- lapply(ids, function(id) sort(unique(sub$condition[sub$ID == id])))

      membership_df <- data.frame(
        `GO ID` = ids, `GO Term` = vapply(ids, term_desc, character(1)),
        `N Conditions` = vapply(membership, length, integer(1)),
        Conditions = vapply(membership, paste, character(1), collapse = "; "),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      membership_df <- membership_df[order(-membership_df$`N Conditions`, membership_df$`GO ID`), ]
      write.csv(membership_df, file.path(output_dir, sprintf("upset_membership_%s_%s.csv", dir, ont)), row.names = FALSE)

      tryCatch({
        upset_df <- data.frame(ID = ids)
        upset_df$Conditions <- membership
        p <- ggplot(upset_df, aes(x = Conditions)) +
          geom_bar(fill = if (dir == "UP") up_color else down_color) +
          scale_x_upset(order_by = "freq") +
          labs(title = sprintf("UpSet: %s-regulated GO terms (%s)", dir, ont), x = NULL, y = "# GO terms") +
          theme_combmatrix(combmatrix.label.text = element_text(size = 9))
        ggsave(file.path(output_dir, sprintf("upset_%s_%s.png", dir, ont)),
               plot = p, width = max(8, 1.2 * length(conditions) + 3), height = 6, bg = "white", limitsize = FALSE)
        cat(sprintf("  [viz] UpSet plot saved: upset_%s_%s.png (%d terms)\n", dir, ont, length(ids)))
      }, error = function(e) cat(paste("  [viz] UpSet plot (", dir, ") failed:", e$message, "\n")))
    }
  }
}

for (ont in go_ontologies) run_for_ontology(ont)

# --- Phase 5: 카테고리별 CSV 취합 -> final_cross_condition_go_results.xlsx ---
# 11_run_group_enrichment.R과 같은 패턴(카테고리별 CSV를 다시 읽어 워크북 시트로 취합)이지만,
# common/flip/exclusive/mixed가 서로 다른 컬럼 스키마를 가진다는 점이 다르다(단일 비교
# GO 결과처럼 균일한 GeneSet x Ontology 그리드가 아님). 그래서 모든 시트에 공통으로
# Category/Ontology 컬럼을 추가해 시트 간 식별 및 이후 parquet 결합의 공통 키로 쓴다.
# flip_*/exclusive_*/mixed_* 라벨은 groups/flips/exclusives 설정에 따라 개수가 달라지므로
# 여기서 다시 생성한다(run_for_ontology() 내부의 라벨 생성 로직과 동일해야 함).
cat("\n[12_run_cross_condition_comparison] Building final_cross_condition_go_results.xlsx...\n")

read_csv_if_nonempty <- function(path) {
  if (!file.exists(path)) return(NULL)
  raw <- trimws(paste(readLines(path, warn = FALSE), collapse = ""))
  if (raw == "" || raw == '""') return(NULL)
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(d) == 0) return(NULL)
  d
}

flip_labels <- vapply(cc_cfg$flips %||% list(), function(flip) {
  dir_from <- toupper(flip[[3]]); dir_to <- setdiff(c("UP", "DOWN"), dir_from)
  sprintf("flip_%s_%s_to_%s_%s", dir_from, flip[[1]], dir_to, flip[[2]])
}, character(1))
exclusive_labels <- vapply(cc_cfg$exclusives %||% list(), function(spec) {
  paste0("exclusive_", spec$label %||% paste(spec$target, collapse = "-"))
}, character(1))
mixed_labels <- paste0("mixed_", names(cc_cfg$groups %||% list()))
category_labels <- c("common_up_strict", "common_up_loose", "common_down_strict", "common_down_loose",
                      "common_up_strict_rrvgo", "common_down_strict_rrvgo",
                      flip_labels, exclusive_labels, mixed_labels)

wb <- createWorkbook()
header_style <- createStyle(fontSize = 11, fontName = "Arial", textDecoration = "bold",
                             halign = "center", valign = "center", fgFill = "#4472C4",
                             fontColour = "#FFFFFF", border = "TopBottomLeftRight", borderColour = "#000000")
text_style <- createStyle(fontSize = 10, fontName = "Arial", halign = "left", valign = "center",
                           border = "TopBottomLeftRight", borderColour = "#D3D3D3")
pvalue_style <- createStyle(fontSize = 10, fontName = "Arial", halign = "right", valign = "center",
                             border = "TopBottomLeftRight", borderColour = "#D3D3D3", numFmt = "0.000")

sheet_summary <- list()
# parquet 결합용 사본(원본 xlsx 시트는 사람이 읽기 좋은 "GO ID"/"GO Term" 헤더를 유지하고,
# parquet만 term_id/description으로 통일 — 06_export_seqviewer.R의 GO+KEGG 결합과 동일한 이유).
sheet_data_for_parquet <- list()

for (label in category_labels) {
  for (ont in go_ontologies) {
    path <- file.path(output_dir, paste0(label, "_", ont, ".csv"))
    d <- read_csv_if_nonempty(path)
    if (is.null(d)) next
    d <- data.frame(Category = label, Ontology = ont, d, check.names = FALSE, stringsAsFactors = FALSE)

    sheet_name <- substr(paste0(label, "_", ont), 1, 31)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, d, headerStyle = header_style)
    addStyle(wb, sheet_name, header_style, rows = 1, cols = seq_len(ncol(d)), gridExpand = TRUE)
    pval_cols <- grep("P-value|P value|Adj P|Jaccard", colnames(d))
    text_cols <- setdiff(seq_len(ncol(d)), pval_cols)
    if (length(text_cols) > 0) addStyle(wb, sheet_name, text_style, rows = 2:(nrow(d) + 1), cols = text_cols, gridExpand = TRUE)
    if (length(pval_cols) > 0) addStyle(wb, sheet_name, pvalue_style, rows = 2:(nrow(d) + 1), cols = pval_cols, gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(d)), widths = "auto")
    sheet_summary[[sheet_name]] <- nrow(d)

    d_parquet <- d %>%
      rename(any_of(c(term_id = "GO ID", term_id = "go", description = "GO Term", description = "term")))
    sheet_data_for_parquet[[sheet_name]] <- d_parquet
  }
}

if (length(sheet_summary) == 0) {
  addWorksheet(wb, "No Results")
  writeData(wb, "No Results", data.frame(Message = "No cross-condition GO categories produced results."))
}

info_df <- data.frame(
  Parameter = c("Conditions", "GO ontologies", "FDR cutoff", "Fold enrichment cutoff",
                "min_conditions_common", "Groups configured", "Flips configured",
                "Exclusives configured", "Sheets with results"),
  Value = c(paste(conditions, collapse = ", "), paste(go_ontologies, collapse = ", "),
            fdr_cutoff, fe_cutoff, min_conditions_common, length(cc_cfg$groups %||% list()),
            length(cc_cfg$flips %||% list()), length(cc_cfg$exclusives %||% list()),
            length(sheet_summary)),
  stringsAsFactors = FALSE
)
addWorksheet(wb, "Analysis_Info")
writeData(wb, "Analysis_Info", info_df)
addStyle(wb, "Analysis_Info", header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "Analysis_Info", cols = 1:2, widths = c(25, 40))

cross_condition_xlsx <- file.path(output_dir, "final_cross_condition_go_results.xlsx")
saveWorkbook(wb, cross_condition_xlsx, overwrite = TRUE)
cat(paste("[12_run_cross_condition_comparison] final_cross_condition_go_results.xlsx saved:",
          cross_condition_xlsx, "(", length(sheet_summary), "sheets with data)\n"))

# --- Phase 6: CMG-SeqViewer export (parquet + staging JSON) — export_seqviewer: true 일 때만 ---
# 06_export_seqviewer.R의 GO/KEGG parquet 결합과 같은 디렉토리 규격을 쓰되, 이 스크립트
# 안에서 이미 만든 sheet_data_for_parquet를 그대로 재사용한다(01c/10처럼 인라인 헬퍼 중복 —
# 이 저장소의 기존 관행). Snakemake output 계약을 지키기 위해 결과가 0건이어도 빈 entries
# JSON은 항상 생성한다.
if (isTRUE(cc_cfg$export_seqviewer)) {
  suppressPackageStartupMessages({
    library(arrow)
    library(jsonlite)
  })

  make_alias_slug <- function(alias, max_len = 80) {
    slug <- gsub("[^\\w가-힣]+", "_", alias, perl = TRUE)
    slug <- gsub("^_|_$", "", slug)
    substr(slug, 1, max_len)
  }
  new_uuid <- function() {
    hex <- paste0(sample(c(0:9, letters[1:6]), 32, replace = TRUE), collapse = "")
    paste(substr(hex, 1, 8), substr(hex, 9, 12),
          paste0("4", substr(hex, 14, 16)),
          paste0(sample(c("8", "9", "a", "b"), 1), substr(hex, 18, 20)),
          substr(hex, 21, 32), sep = "-")
  }
  write_parquet_dataset <- function(df, alias, datasets_dir) {
    uid      <- new_uuid()
    slug     <- make_alias_slug(alias)
    filename <- paste0(slug, ".parquet")
    old_files <- list.files(datasets_dir, pattern = paste0("^", slug, "\\.parquet$"), full.names = TRUE)
    if (length(old_files) > 0) file.remove(old_files)
    write_parquet(df, file.path(datasets_dir, filename))
    list(uid = uid, filename = filename)
  }

  seqviewer_dir <- file.path(config$output_dir, "seqviewer")
  datasets_dir  <- file.path(seqviewer_dir, "datasets")
  staging_dir   <- file.path(seqviewer_dir, "staging")
  dir.create(datasets_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(staging_dir,  recursive = TRUE, showWarnings = FALSE)
  staging_path <- file.path(staging_dir, "cross_condition_entries.json")

  if (length(sheet_data_for_parquet) > 0) {
    cc_combined <- suppressWarnings(bind_rows(sheet_data_for_parquet))
    cc_alias <- paste(basename(config$output_dir), "Cross-Condition GO")
    cc_info  <- write_parquet_dataset(cc_combined, cc_alias, datasets_dir)

    cc_entry <- list(
      dataset_id           = cc_info$uid,
      alias                = cc_alias,
      original_filename    = cc_info$filename,
      dataset_type         = "cross_condition_go",
      experiment_condition = paste(conditions, collapse = ", "),
      organism             = config$species %||% "",
      cell_type            = "",
      tissue               = "",
      timepoint            = "",
      row_count            = nrow(cc_combined),
      gene_count           = 0L,
      significant_genes    = 0L,
      import_date          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      file_path            = cc_info$filename,
      notes                = sprintf("Cross-condition GO comparison (FDR<%.3g, FE>%.1f, %d conditions): common/flip/exclusive/mixed categories combined",
                                      fdr_cutoff, fe_cutoff, length(conditions)),
      tags                 = as.list(c("cross_condition", "GO", conditions))
    )
    write_json(list(cc_entry), staging_path, pretty = TRUE, auto_unbox = TRUE)
    cat(paste("[12_run_cross_condition_comparison] Seqviewer parquet saved:", cc_info$filename, "\n"))
    cat(paste("[12_run_cross_condition_comparison] Seqviewer staging JSON saved:", staging_path, "\n"))
  } else {
    write_json(list(), staging_path, pretty = TRUE, auto_unbox = TRUE)
    cat("[12_run_cross_condition_comparison] No cross-condition sheets with data — wrote empty staging JSON.\n")
  }
}

writeLines(condition_count_log, file.path(output_dir, "condition_count_log.txt"))
cat("\n[12_run_cross_condition_comparison] Done.\n")
