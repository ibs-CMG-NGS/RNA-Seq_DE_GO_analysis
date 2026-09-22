#!/usr/bin/env Rscript
# 19_run_cluster_cross_dataset_comparison.R
#
# 18_run_cross_dataset_go_comparison.R은 "pairwise DE 비교"(조건 vs Control) 단위로
# 데이터셋 간 GO term을 비교했다. 이 스크립트는 그 축을 한 단계 더 확장해서
# maSigPro time-series 클러스터(01c) / DEGreport coexpression module(10) 각각을
# "같은 반응 패턴을 보이는 유전자 묶음"으로 보고, 데이터셋 쌍마다 모든 클러스터/모듈
# 사이의 GO term-set Jaccard 유사도를 전부 계산한다 — 어떤 데이터셋 A의 클러스터가
# 어떤 데이터셋 B의 클러스터와 기능적으로 가장 비슷한지 찾는 다대다 매칭 문제.
# 재계산 없이 11_run_group_enrichment.R이 이미 저장한
# go_termcluster_cluster{N}_BP.csv / go_termcluster_module{N}_BP.csv(raw enrichGO
# 결과, term_cluster 이전의 전체 유의 term 목록)만 재사용한다.
#
# GO ID는 species-agnostic 온톨로지라 ortholog 매핑 없이 term-set 그대로 Jaccard
# 비교 가능(18과 동일 전제) — 이 스크립트는 gene-level이 아니라 GO-ID-level Jaccard만
# 쓰므로 Entrez/Symbol 변환이 필요 없다.
#
# 참고 원형: atac-seq-da-analysis/src/analysis/28_run_cluster_cross_species_comparison.R
# (그쪽은 mouse 1개 + human 1개로 하드코딩되어 있었음 — 여기서는 데이터셋 2개 이상
# 임의 조합으로 일반화: combn(dataset_labels, 2)로 모든 쌍에 대해 반복 수행.)
#
# 사용법: Rscript 19_run_cluster_cross_dataset_comparison.R <cross_dataset_config.yaml>

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 19_run_cluster_cross_dataset_comparison.R <cross_dataset_config.yaml>")
}
cross_cfg <- yaml::read_yaml(args[1])
output_dir <- cross_cfg$output_dir %||% stop("[FATAL] output_dir not set in config.")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

datasets_cfg <- cross_cfg$datasets
if (is.null(datasets_cfg) || length(datasets_cfg) < 2) {
  stop("[FATAL] at least 2 datasets are required in config$datasets.")
}
dataset_labels <- vapply(datasets_cfg, function(d) d$label, character(1))
project_cfgs <- lapply(datasets_cfg, function(d) yaml::read_yaml(d$config))
names(project_cfgs) <- dataset_labels
project_dir_of <- setNames(vapply(dataset_labels, function(lbl) project_cfgs[[lbl]]$output_dir, character(1)), dataset_labels)

fdr_cutoff <- cross_cfg$fdr_cutoff %||% 0.05

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d
}

# --- {project_dir}/time_series[_{variant}]/ 또는 coexpression_modules[_{variant}]/ 를
# 스캔해 "{variant}_{cluster_id}" -> 유의 GO term ID 집합" 수집 ---
collect_cluster_term_sets <- function(project_dir, mode) {
  if (mode == "ts") {
    prefix <- "time_series"
    file_pattern <- "^go_termcluster_cluster(\\w+)_BP\\.csv$"
  } else {
    prefix <- "coexpression_modules"
    file_pattern <- "^go_termcluster_module(\\w+)_BP\\.csv$"
  }
  all_dirs <- list.dirs(project_dir, recursive = FALSE)
  is_match <- grepl(paste0("^", prefix, "($|_)"), basename(all_dirs))
  variant_dirs <- all_dirs[is_match]
  variant_labels <- sub(paste0("^", prefix, "_?"), "", basename(variant_dirs))
  variant_labels[variant_labels == ""] <- "base"

  result <- list()
  for (i in seq_along(variant_dirs)) {
    files <- list.files(variant_dirs[i], pattern = file_pattern, full.names = TRUE)
    for (f in files) {
      cid <- sub(file_pattern, "\\1", basename(f))
      d <- read_csv_safe(f)
      if (is.null(d) || !("ID" %in% colnames(d)) || !("p.adjust" %in% colnames(d))) next
      ids <- unique(d$ID[!is.na(d$p.adjust) & d$p.adjust < fdr_cutoff])
      if (length(ids) == 0) next
      key <- sprintf("%s_%s", variant_labels[i], cid)
      result[[key]] <- ids
    }
  }
  result
}

jaccard_matrix <- function(set_a, set_b) {
  mat <- matrix(0, nrow = length(set_a), ncol = length(set_b),
                 dimnames = list(names(set_a), names(set_b)))
  for (i in names(set_a)) for (j in names(set_b)) {
    u <- union(set_a[[i]], set_b[[j]])
    mat[i, j] <- if (length(u) == 0) 0 else length(intersect(set_a[[i]], set_b[[j]])) / length(u)
  }
  mat
}

run_comparison_for_pair <- function(label_a, label_b, mode, title_prefix, file_prefix) {
  message(sprintf("\n=== %s: %s vs %s ===", title_prefix, label_a, label_b))
  sets_a <- collect_cluster_term_sets(project_dir_of[[label_a]], mode)
  sets_b <- collect_cluster_term_sets(project_dir_of[[label_b]], mode)
  message(sprintf("  %s: %d개 유의 클러스터/모듈 (%s)", label_a, length(sets_a), paste(names(sets_a), collapse = ", ")))
  message(sprintf("  %s: %d개 유의 클러스터/모듈 (%s)", label_b, length(sets_b), paste(names(sets_b), collapse = ", ")))
  if (length(sets_a) == 0 || length(sets_b) == 0) {
    message("  비교 불가(한쪽에 유의 클러스터/모듈 없음) — 건너뜀")
    return(invisible(NULL))
  }

  mat <- jaccard_matrix(sets_a, sets_b)
  mat_df <- as.data.frame(mat)
  colnames(mat_df) <- names(sets_b)
  mat_df <- data.frame(setNames(list(rownames(mat)), sprintf("%s_cluster", label_a)),
                        mat_df, check.names = FALSE, stringsAsFactors = FALSE)
  out_prefix <- sprintf("%s_%s_vs_%s", file_prefix, label_a, label_b)
  write.csv(mat_df, file.path(output_dir, sprintf("%s_jaccard_matrix.csv", out_prefix)), row.names = FALSE)

  best_rows <- lapply(rownames(mat), function(i) {
    j_best <- colnames(mat)[which.max(mat[i, ])]
    score <- max(mat[i, ])
    n_shared <- length(intersect(sets_a[[i]], sets_b[[j_best]]))
    data.frame(a_cluster = i, b_best_match = j_best, jaccard = round(score, 3),
               n_shared_terms = n_shared, n_a_terms = length(sets_a[[i]]),
               n_b_terms = length(sets_b[[j_best]]), stringsAsFactors = FALSE)
  })
  best_df <- do.call(rbind, best_rows)
  colnames(best_df) <- c(sprintf("%s_cluster", label_a), sprintf("%s_best_match", label_b),
                          "jaccard", "n_shared_terms", sprintf("n_%s_terms", label_a), sprintf("n_%s_terms", label_b))
  write.csv(best_df, file.path(output_dir, sprintf("%s_best_match.csv", out_prefix)), row.names = FALSE)
  message("  best match:")
  for (i in seq_len(nrow(best_df))) {
    message(sprintf("    %s <-> %s : Jaccard=%.3f (공유 %d / %s %d / %s %d term)",
                     best_df[i, 1], best_df[i, 2], best_df$jaccard[i],
                     best_df$n_shared_terms[i], label_a, best_df[i, 5], label_b, best_df[i, 6]))
  }

  long_df <- as.data.frame(as.table(mat))
  colnames(long_df) <- c("a_cluster", "b_cluster", "Jaccard")
  p <- ggplot(long_df, aes(x = b_cluster, y = a_cluster, fill = Jaccard)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 3.5,
              color = ifelse(long_df$Jaccard > 0.3, "white", "black")) +
    scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, max(1e-6, max(long_df$Jaccard)))) +
    labs(x = sprintf("%s cluster/module", label_b), y = sprintf("%s cluster/module", label_a),
         title = sprintf("%s: GO BP Term-Set Jaccard Similarity (%s vs %s)", title_prefix, label_a, label_b),
         subtitle = "값 = |공통 유의 GO term| / |합집합 유의 GO term|") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(output_dir, sprintf("%s_jaccard_heatmap.png", out_prefix)),
         plot = p, width = max(8, 1.1 * length(sets_b) + 2), height = max(5, 0.9 * length(sets_a) + 2),
         bg = "white")
  message(sprintf("  heatmap 저장: %s_jaccard_heatmap.png", out_prefix))
}

dataset_pairs <- if (length(dataset_labels) >= 2) combn(dataset_labels, 2, simplify = FALSE) else list()
for (lp in dataset_pairs) {
  run_comparison_for_pair(lp[1], lp[2], "ts", "Time-Series Cluster Comparison", "ts_cluster")
  run_comparison_for_pair(lp[1], lp[2], "coexpr", "Co-expression Module Comparison", "coexpr_module")
}

prov_df <- data.frame(
  Comparison = c("ts_cluster", "coexpr_module"),
  `Dataset Pairs` = paste(vapply(dataset_pairs, paste, character(1), collapse = " vs "), collapse = "; "),
  `Source Pattern` = c(
    "{project}/time_series[_{variant}]/go_termcluster_cluster{N}_BP.csv",
    "{project}/coexpression_modules[_{variant}]/go_termcluster_module{N}_BP.csv"
  ),
  Processing = sprintf("각 파일의 유의 GO ID(p.adjust<%.3g) 집합 추출 -> 데이터셋 쌍마다 전체 조합 Jaccard 유사도 행렬 -> 클러스터별 최고 매칭 상대 클러스터 도출(재계산 없음, 01c/10/11 raw enrichGO 결과 재사용)", fdr_cutoff),
  check.names = FALSE, stringsAsFactors = FALSE
)
write.csv(prov_df, file.path(output_dir, "cluster_cross_dataset_provenance.csv"), row.names = FALSE)

message("\n[19_run_cluster_cross_dataset_comparison] Done.")
