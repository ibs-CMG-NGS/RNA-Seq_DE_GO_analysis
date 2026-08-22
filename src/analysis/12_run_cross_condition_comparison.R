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
dotplot_max_terms <- viz_cfg$dotplot_max_terms %||% 40
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
    file.path(project_dir, "pairwise", condition, paste0("kegg_enrichment_", direction, ".csv"))
  } else {
    file.path(project_dir, "pairwise", condition, paste0("go_enrichment_", direction, "_", ont, ".csv"))
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

    render_plot <- function(ids, out_name_dot, out_name_heatmap, title) {
      if (length(ids) == 0) return(invisible(NULL))
      term_rank <- vapply(ids, function(id) min(all_df$p.adjust[all_df$ID == id]), numeric(1))
      ids <- ids[order(term_rank)]

      if (cluster_terms_enabled && length(ids) >= 3) {
        hc <- tryCatch(hclust(dist(signed_matrix_for(ids)), method = cluster_method), error = function(e) NULL)
        if (!is.null(hc)) ids <- ids[hc$order]
      }
      trunc_term <- function(x) ifelse(nchar(x) > 55, paste0(substr(x, 1, 52), "..."), x)

      if (length(ids) <= dotplot_max_terms) {
        cat(sprintf("  [viz] Dot plot: %d terms (<= %d) -> %s\n", length(ids), dotplot_max_terms, out_name_dot))
        tryCatch({
          plot_rows <- list()
          for (id in ids) {
            desc <- trunc_term(term_desc(id))
            for (cond in conditions) {
              sub <- all_df[all_df$ID == id & all_df$condition == cond, ]
              if (nrow(sub) == 0) next
              best <- sub[which.min(sub$p.adjust), ]
              plot_rows[[length(plot_rows) + 1]] <- data.frame(
                GO_ID = id, GO_Term = desc, Condition = cond,
                GeneRatio = best$GeneRatioNum,
                SignedLog10P = ifelse(best$direction == "UP", 1, -1) * -log10(best$p.adjust),
                IsFlip = id %in% all_flip_ids,
                stringsAsFactors = FALSE
              )
            }
          }
          plot_df <- do.call(rbind, plot_rows)
          term_order <- trunc_term(vapply(ids, term_desc, character(1)))
          plot_df$GO_Term <- factor(plot_df$GO_Term, levels = rev(unique(term_order)))

          p <- ggplot(plot_df, aes(x = Condition, y = GO_Term)) +
            geom_point(aes(size = GeneRatio, color = SignedLog10P)) +
            scale_color_gradient2(low = down_color, mid = "grey85", high = up_color, midpoint = 0,
                                  name = "sign(dir) x -log10(FDR)") +
            scale_size_continuous(name = "Gene Ratio", range = c(1, 8)) +
            theme_bw(base_size = 11) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL, title = title,
                 caption = if (any(plot_df$IsFlip)) "* black outline = direction-flip term" else NULL)
          if (any(plot_df$IsFlip)) {
            p <- p + geom_point(data = plot_df[plot_df$IsFlip, ], shape = 21, size = 3.2,
                                 color = "black", stroke = 1, show.legend = FALSE)
          }
          # 폭 = legend(~2.6in) + 조건당 패널 폭(~0.8in) + y축 라벨 폭(글자 수 비례, 최대 55자 기준)
          label_w <- 0.085 * max(nchar(as.character(levels(plot_df$GO_Term))))
          w <- 2.6 + 0.8 * length(conditions) + label_w
          h <- max(5, 0.28 * length(ids) + 2)
          ggsave(file.path(output_dir, out_name_dot), plot = p, width = w, height = h, bg = "white", limitsize = FALSE)
        }, error = function(e) cat(paste("  [viz] dot plot failed:", e$message, "\n")))
      } else {
        cat(sprintf("  [viz] %d terms (> %d) -> heatmap fallback %s\n", length(ids), dotplot_max_terms, out_name_heatmap))
        tryCatch({
          raw_labels <- vapply(ids, term_desc, character(1))
          trunc_labels <- ifelse(nchar(raw_labels) > 55, paste0(substr(raw_labels, 1, 52), "..."), raw_labels)
          term_labels <- make.unique(trunc_labels)
          mat <- signed_matrix_for(ids)
          rownames(mat) <- term_labels
          max_abs <- max(abs(mat), 1)
          breaks <- seq(-max_abs, max_abs, length.out = 102)
          colors <- colorRampPalette(c(down_color, "white", up_color))(101)
          n_rows <- length(ids)
          row_fontsize <- if (n_rows > 100) 5 else if (n_rows > 60) 6 else 8
          # cluster_rows=TRUE면 pheatmap이 자체적으로 재클러스터링해서 덴드로그램을 그린다
          # (dot plot과 같은 signed -log10(FDR) 행렬 기준이라 사실상 같은 패턴으로 묶임).
          img_width <- max(1400, 480 + max(nchar(term_labels)) * 8 + 130 * length(conditions) + 60)
          img_height <- max(900, 20 * n_rows + 250)
          png(file.path(output_dir, out_name_heatmap), width = img_width, height = img_height, res = 150)
          pheatmap(mat, color = colors, breaks = breaks, cluster_cols = FALSE,
                   cluster_rows = cluster_terms_enabled, clustering_method = cluster_method,
                   fontsize_row = row_fontsize, fontsize_col = 9, main = title)
          dev.off()
        }, error = function(e) {
          if (dev.cur() > 1) dev.off()
          cat(paste("  [viz] heatmap failed:", e$message, "\n"))
        })
      }
    }

    candidate_ids <- unique(c(common_strict_ids[["UP"]], common_strict_ids[["DOWN"]], all_flip_ids, all_excl_ids))
    if (length(candidate_ids) == 0) candidate_ids <- unique(c(common_loose_ids[["UP"]], common_loose_ids[["DOWN"]]))
    if (length(candidate_ids) == 0) {
      by_p <- all_df[order(all_df$p.adjust), ]
      candidate_ids <- unique(head(by_p$ID, 20))
    }
    render_plot(candidate_ids,
                sprintf("cross_condition_dotplot_%s.png", ont), sprintf("cross_condition_heatmap_%s.png", ont),
                sprintf("Cross-Condition GO Dot Plot (%s)", ont))

    if (length(semantic_representative_ids) > 0) {
      render_plot(semantic_representative_ids,
                  sprintf("cross_condition_dotplot_semantic_%s.png", ont), sprintf("cross_condition_heatmap_semantic_%s.png", ont),
                  sprintf("Cross-Condition GO Dot Plot - Semantic Filtered (%s)", ont))
    }

    # --- UpSet plot: 조건별 유의 term 중첩 구조 (UP/DOWN 각각) ---
    for (dir in c("UP", "DOWN")) {
      sub <- all_df[all_df$direction == dir, ]
      if (length(unique(sub$ID)) == 0) next
      tryCatch({
        ids <- unique(sub$ID)
        upset_df <- data.frame(ID = ids)
        upset_df$Conditions <- lapply(ids, function(id) sort(unique(sub$condition[sub$ID == id])))
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

writeLines(condition_count_log, file.path(output_dir, "condition_count_log.txt"))
cat("\n[12_run_cross_condition_comparison] Done.\n")
