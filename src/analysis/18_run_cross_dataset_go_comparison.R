#!/usr/bin/env Rscript
# 18_run_cross_dataset_go_comparison.R
#
# Cross-dataset(프로젝트 간) GO term 비교. 12_run_cross_condition_comparison.R의
# A(common)/B(flip)/C(exclusive)/D(mixed) 로직을 "조건"이 아니라 "데이터셋(프로젝트)"을
# 축으로 재사용한다 — 이미 완료된 프로젝트 2개 이상(다른 종이어도 무방)을 사후에
# 비교하는 메타분석이라 각 프로젝트의 Snakemake rule all에는 엮지 않고 독립 실행한다
# (참고: atac-seq-da-analysis/docs/plan_cross_dataset_atac_comparison.md와 동일 설계
# 결정, atac-seq-da-analysis/src/analysis/20_run_cross_species_go_comparison.R이 이
# 스크립트의 직접적인 원형).
#
# 12와의 핵심 차이:
#   - 12는 한 프로젝트의 project_dir가 고정이고 "조건"만 여러 개
#     (pairwise/{condition}/enrichment/go_enrichment_*.csv). 여기서는 프로젝트
#     (project_dir)가 여러 개이므로 "조건" 대신 (dataset_label, pair) 조합을 합성
#     condition id "{dataset_label}::{pair}"로 취급한다 — 나머지 common/flip/
#     exclusive/mixed/rrvgo/dot plot/heatmap/UpSet 로직은 12와 거의 동일하다.
#   - groups는 데이터셋 라벨로 고정(데이터셋 자체가 그룹) — "common"은 "모든
#     데이터셋에서 공통", "flip"은 "한 데이터셋에서는 UP인데 다른 데이터셋에서는 DOWN".
#   - KEGG pathway ID는 organism prefix(hsaNNNNN/mmuNNNNN)가 붙어 있어 그대로
#     비교하면 항상 불일치하므로, prefix를 제거한 pathway 번호를 비교 키로 정규화한다
#     (같은 reference pathway를 가리키는 서로 다른 organism instance라고 가정).
#   - 이 저장소의 enrichGO/enrichKEGG geneID 컬럼은 (03_enrichment_analysis.R이
#     readable=FALSE로 돌려서) Entrez ID다 — 종이 다르면 Entrez 번호 자체가 무관하므로
#     gene-set Jaccard(B. flip) 계산 전에 데이터셋별 organism_db로 SYMBOL 매핑을 먼저
#     한다. 그 다음 종간 대소문자 표기 차이(mouse "Bdnf" vs human "BDNF")를 흡수하려고
#     toupper()로 정규화한다 — 이는 문자 그대로의 심볼 일치일 뿐 정식 ortholog 매핑이
#     아니므로 참고용 지표로만 사용할 것.
#   - rrvgo(semantic_filtering)의 IC(information content)는 종별 GO 주석 빈도에서
#     계산되므로 원래 종마다 다르다 — 여기서는 datasets[1](첫 번째 데이터셋)의
#     organism_db 하나만 사용해 근사한다(한계, 주석 필요).
#
# 신규 시각화(12에는 없음): pair별 두 데이터셋 직접 산점도 — 각 pair(예:
# Acute_1D_vs_Control)마다 데이터셋 A의 signed -log10(FDR) vs 데이터셋 B의
# signed -log10(FDR)을 GO term 단위로 그려 concordant(같은 방향)/discordant(반대
# 방향)/데이터셋-특이 term을 구분한다.
#
# 사용법: Rscript 18_run_cross_dataset_go_comparison.R <cross_dataset_config.yaml>
#   예) Rscript 18_run_cross_dataset_go_comparison.R \
#         configs/cross_dataset_hiy-mouse-vs-human.yaml

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
  library(ggupset)
  library(pheatmap)
  library(dplyr)
  library(openxlsx)
  library(AnnotationDbi)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- assay 프리셋: RNA-Seq_DE_GO_analysis와 atac-seq-da-analysis는 GO enrichment
# 표(ID/Description/GeneRatio/BgRatio/p.adjust/geneID/Count) 자체는 동일하지만
# pairwise 결과를 어디서 찾을지·geneID가 어떤 포맷인지가 다르다(직접 대조 확인,
# docs/atac_pipeline_alignment_request.md 참고). datasets[].assay(기본 "rna")로
# 데이터셋별로 선택한다.
ASSAY_PRESETS <- list(
  rna = list(
    geneid_format = "entrez",
    pairwise_go_path = function(pd, pair, dir, ont)
      file.path(pd, "pairwise", pair, "enrichment", sprintf("go_enrichment_%s_%s.csv", dir, ont)),
    pairwise_kegg_path = function(pd, pair, dir)
      file.path(pd, "pairwise", pair, "enrichment", sprintf("kegg_enrichment_%s.csv", dir))
  ),
  atac = list(
    geneid_format = "symbol",
    pairwise_go_path = function(pd, pair, dir, ont)
      file.path(pd, "pairwise", pair, sprintf("go_enrichment_%s_%s.csv", dir, ont)),
    pairwise_kegg_path = function(pd, pair, dir)
      file.path(pd, "pairwise", pair, sprintf("kegg_enrichment_%s.csv", dir))
  )
)

# --- 1. 인자 파싱 & config 로드 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 18_run_cross_dataset_go_comparison.R <cross_dataset_config.yaml>")
}
cross_cfg <- yaml::read_yaml(args[1])

output_dir <- cross_cfg$output_dir %||% stop("[FATAL] output_dir not set in config.")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

datasets_cfg <- cross_cfg$datasets
if (is.null(datasets_cfg) || length(datasets_cfg) < 2) {
  stop("[FATAL] at least 2 datasets are required in config$datasets.")
}
pairs <- cross_cfg$pairs
if (is.null(pairs) || length(pairs) == 0) stop("[FATAL] config$pairs is empty.")

go_ontologies <- cross_cfg$go_ontologies %||% c("BP")
fdr_cutoff <- cross_cfg$fdr_cutoff %||% 0.05
fe_cutoff  <- cross_cfg$fold_enrichment_cutoff %||% 2.0
sparsity_warn_threshold <- cross_cfg$sparsity_warn_threshold %||% 20
jaccard_warn <- cross_cfg$gene_overlap_jaccard_warn %||% 0.3

semantic_cfg <- cross_cfg$semantic_filtering %||% list()
semantic_enabled <- if (is.null(semantic_cfg$enabled)) TRUE else isTRUE(semantic_cfg$enabled)

viz_cfg <- cross_cfg$visualization %||% list()
viz_enabled <- if (is.null(viz_cfg$enabled)) TRUE else isTRUE(viz_cfg$enabled)
dotplot_max_terms <- viz_cfg$dotplot_max_terms %||% 60
cluster_terms_enabled <- if (is.null(viz_cfg$cluster_terms)) TRUE else isTRUE(viz_cfg$cluster_terms)
cluster_method <- viz_cfg$cluster_method %||% "average"

# --- 2. 데이터셋 메타(각 프로젝트 config에서 output_dir/species/organism_db 로드) ---
dataset_labels <- vapply(datasets_cfg, function(d) d$label, character(1))
if (length(dataset_labels) != length(unique(dataset_labels))) {
  stop("[FATAL] duplicate dataset labels in config$datasets.")
}
project_cfgs <- lapply(datasets_cfg, function(d) yaml::read_yaml(d$config))
names(project_cfgs) <- dataset_labels

# 프로젝트 config의 output_dir은 그 프로젝트의 "본가" 레포(RNA-Seq_DE_GO_analysis
# 또는 atac-seq-da-analysis) 작업 디렉토리를 기준으로 한 상대경로로 적혀 있다(각
# 레포의 Snakefile/스크립트가 항상 그 레포 루트에서 실행된다는 전제). 이 cross-dataset
# 스크립트는 항상 RNA-Seq_DE_GO_analysis에서 실행되므로, 다른 레포(ATAC)의 상대
# output_dir을 그대로 file.path()하면 엉뚱한(RNA 레포 기준) 경로가 되어 조용히 0건
# 매칭 실패로 이어진다 — datasets[].config 파일이 위치한 레포 루트(config 경로의
# 조부모 디렉토리, "<repo>/configs/config_X.yml" 관례)를 기준으로 resolve한다.
is_abs_path <- function(p) grepl("^/", p)
config_repo_root <- function(config_path) dirname(dirname(normalizePath(config_path)))

dataset_meta <- setNames(lapply(seq_along(dataset_labels), function(i) {
  lbl <- dataset_labels[i]
  pcfg <- project_cfgs[[lbl]]
  species <- pcfg$species
  assay <- datasets_cfg[[i]]$assay %||% "rna"
  if (!assay %in% names(ASSAY_PRESETS)) {
    stop(sprintf("[FATAL] dataset '%s': unknown assay '%s' (must be one of: %s)",
                  lbl, assay, paste(names(ASSAY_PRESETS), collapse = ", ")))
  }
  raw_output_dir <- pcfg$output_dir
  project_dir <- if (is_abs_path(raw_output_dir)) {
    raw_output_dir
  } else {
    file.path(config_repo_root(datasets_cfg[[i]]$config), raw_output_dir)
  }
  list(label = lbl, project_dir = project_dir, species = species,
       organism_db = pcfg$databases[[species]]$organism_db,
       assay = assay, preset = ASSAY_PRESETS[[assay]])
}), dataset_labels)

for (lbl in dataset_labels) {
  if (!suppressPackageStartupMessages(
        requireNamespace(dataset_meta[[lbl]]$organism_db, quietly = TRUE))) {
    stop(sprintf("[FATAL] Required organism DB package %s (dataset '%s') is not installed.",
                  dataset_meta[[lbl]]$organism_db, lbl))
  }
  # get()으로 바로 쓰려면(Entrez->SYMBOL 매핑) requireNamespace만으로는 부족 —
  # 네임스페이스를 실제로 attach한다.
  suppressPackageStartupMessages(
    library(dataset_meta[[lbl]]$organism_db, character.only = TRUE))
}

primary_organism_db <- dataset_meta[[dataset_labels[1]]]$organism_db
rrvgo_global_cfg <- project_cfgs[[dataset_labels[1]]]$enrichment$rrvgo %||% list()
rr_method    <- semantic_cfg$method    %||% rrvgo_global_cfg$method    %||% "Rel"
rr_threshold <- semantic_cfg$threshold %||% rrvgo_global_cfg$threshold %||% 0.7

up_color   <- project_cfgs[[dataset_labels[1]]]$plot_aesthetics$volcano$up_color   %||% "#FF5733"
down_color <- project_cfgs[[dataset_labels[1]]]$plot_aesthetics$volcano$down_color %||% "#3375FF"

message(sprintf("[18_run_cross_dataset_go_comparison] Datasets (%d): %s",
                 length(dataset_labels), paste(dataset_labels, collapse = ", ")))
message(sprintf("[18_run_cross_dataset_go_comparison] Pairs (%d): %s",
                 length(pairs), paste(pairs, collapse = ", ")))
message(sprintf("[18_run_cross_dataset_go_comparison] GO ontologies: %s",
                 paste(go_ontologies, collapse = ", ")))

# --- 2b. 데이터셋별 Entrez -> SYMBOL 매핑 테이블 (gene-level Jaccard용) ---
# 프로젝트별로 한 번만 organism_db 전체 ENTREZID<->SYMBOL 매핑을 읽어 캐시한다
# (per-row mapIds 호출은 느리고, 어차피 필요한 ID는 그 프로젝트가 다루는 전체 유전자
# 범위 안에 있으므로 keys=keytype 없이 keys 지정 없는 전체 dump가 오히려 더 빠르다).
# geneid_format이 이미 "symbol"인 데이터셋(예: atac 프리셋)은 매핑이 필요 없으므로
# 건너뛴다 — 애초에 org db 전체를 훑는 이 작업 자체가 그 데이터셋엔 낭비.
entrez2symbol <- setNames(lapply(dataset_labels, function(lbl) {
  if (dataset_meta[[lbl]]$preset$geneid_format != "entrez") return(NULL)
  db <- get(dataset_meta[[lbl]]$organism_db)
  suppressMessages(AnnotationDbi::select(db, keys = keys(db, keytype = "ENTREZID"),
                                          keytype = "ENTREZID", columns = "SYMBOL"))
}), dataset_labels)
for (lbl in dataset_labels) {
  tab <- entrez2symbol[[lbl]]
  if (is.null(tab)) next
  entrez2symbol[[lbl]] <- setNames(tab$SYMBOL, tab$ENTREZID)
}

# --- 3. condition 축 = (dataset_label, pair) 합성 id, groups = 데이터셋 ---
# pairs는 "정식(canonical)" 비교 이름이지만, 프로젝트마다 실제 pairwise 폴더명이
# 다를 수 있다(예: 같은 Acute 1D인데 한 프로젝트는 "Acute_1D_vs_Control", 다른
# 프로젝트는 "D1_vs_Control"). datasets_cfg[[i]]$pair_map(선택)으로 "정식 이름 ->
# 그 데이터셋의 실제 폴더명"을 매핑할 수 있다 — condition id(dataset::pair)는 항상
# 정식 이름을 쓰고, 실제 파일 경로 생성에만 매핑된 이름을 쓴다.
pair_maps <- setNames(lapply(datasets_cfg, function(d) d$pair_map %||% list()), dataset_labels)
resolve_actual_pair <- function(lbl, p) {
  v <- pair_maps[[lbl]][[p]]
  if (is.null(v)) p else v
}

condition_meta <- do.call(rbind, lapply(dataset_labels, function(lbl) {
  actual <- vapply(pairs, resolve_actual_pair, character(1), lbl = lbl)
  data.frame(id = paste0(lbl, "::", pairs), dataset = lbl, pair = pairs, actual_pair = actual,
             stringsAsFactors = FALSE)
}))
conditions <- condition_meta$id
if (length(conditions) < 2) stop("[FATAL] cross-dataset comparison requires at least 2 (dataset, pair) conditions.")

min_conditions_common <- cross_cfg$min_datasets_common %||% round(0.75 * length(conditions))
message(sprintf("[18_run_cross_dataset_go_comparison] min_conditions_common = %d (of %d)",
                 min_conditions_common, length(conditions)))

groups <- split(condition_meta$id, condition_meta$dataset)[dataset_labels]

# flips: 데이터셋 라벨의 모든 조합(unordered pair)마다 UP->DOWN/DOWN->UP 양방향
label_pairs <- if (length(dataset_labels) >= 2) combn(dataset_labels, 2, simplify = FALSE) else list()
flips <- do.call(c, lapply(label_pairs, function(lp) {
  list(list(lp[1], lp[2], "UP"), list(lp[1], lp[2], "DOWN"))
}))

# exclusives: 각 데이터셋에서만(다른 모든 데이터셋에는 전혀 없이) 일관되게 UP/DOWN인 term
exclusives <- do.call(c, lapply(dataset_labels, function(lbl) {
  lapply(c("UP", "DOWN"), function(dir) {
    list(label = sprintf("%s_only_%s", lbl, tolower(dir)),
         target = groups[[lbl]],
         absent = unlist(groups[dataset_labels != lbl], use.names = FALSE),
         direction = dir)
  })
}))

# --- 4. 조건별 pairwise/{pair}/enrichment/go_enrichment_{dir}_{ont}.csv 로드 & 유의 term만 남기기 ---
parse_ratio <- function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1] / parts[2]
}
compute_fold_enrichment <- function(go_df) {
  mapply(function(gr, br) parse_ratio(gr) / parse_ratio(br), go_df$GeneRatio, go_df$BgRatio)
}
# KEGG pathway ID의 organism prefix(hsa/mmu 등) 제거 -> 같은 reference pathway 번호로 정규화
normalize_kegg_id <- function(id) sub("^[a-zA-Z]+", "", id)

entrez_list_to_symbols <- function(gene_id_str, lookup) {
  ids <- strsplit(gene_id_str, "/", fixed = TRUE)[[1]]
  syms <- unname(lookup[ids])
  syms <- syms[!is.na(syms)]
  paste(syms, collapse = "/")
}

load_condition_direction <- function(condition_id, direction, ont) {
  meta_row <- condition_meta[condition_meta$id == condition_id, ]
  ds <- dataset_meta[[meta_row$dataset]]
  pair <- meta_row$actual_pair
  path <- if (toupper(ont) == "KEGG") {
    ds$preset$pairwise_kegg_path(ds$project_dir, pair, direction)
  } else {
    ds$preset$pairwise_go_path(ds$project_dir, pair, direction, ont)
  }
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d$FoldEnrichment <- compute_fold_enrichment(d)
  d <- d[!is.na(d$p.adjust) & d$p.adjust < fdr_cutoff &
         !is.na(d$FoldEnrichment) & d$FoldEnrichment > fe_cutoff, ]
  if (nrow(d) == 0) return(NULL)
  if (toupper(ont) == "KEGG") d$ID <- normalize_kegg_id(d$ID)
  if (ds$preset$geneid_format == "entrez") {
    lookup <- entrez2symbol[[meta_row$dataset]]
    d$geneID <- vapply(d$geneID, entrez_list_to_symbols, character(1), lookup = lookup)
  }
  d$condition <- condition_id
  d$direction <- toupper(direction)
  d$GeneRatioNum <- sapply(d$GeneRatio, parse_ratio)
  d[, c("ID", "Description", "p.adjust", "FoldEnrichment", "GeneRatioNum", "Count", "geneID", "condition", "direction")]
}

condition_count_log <- character()

# --- 5. pair별 두 데이터셋 직접 산점도(신규, 12에는 없음) ---
draw_pair_scatter <- function(all_df, ont) {
  if (length(label_pairs) == 0) return(invisible(NULL))
  thr <- -log10(fdr_cutoff)
  signed_val <- function(id, condition_id) {
    sub <- all_df[all_df$ID == id & all_df$condition == condition_id, ]
    if (nrow(sub) == 0) return(0)
    best <- sub[which.min(sub$p.adjust), ]
    ifelse(best$direction == "UP", 1, -1) * -log10(best$p.adjust)
  }
  for (lp in label_pairs) {
    ds_a <- lp[1]; ds_b <- lp[2]
    for (pr in pairs) {
      cond_a <- paste0(ds_a, "::", pr); cond_b <- paste0(ds_b, "::", pr)
      ids <- unique(all_df$ID[all_df$condition %in% c(cond_a, cond_b)])
      if (length(ids) == 0) next
      val_a <- vapply(ids, signed_val, numeric(1), condition_id = cond_a)
      val_b <- vapply(ids, signed_val, numeric(1), condition_id = cond_b)
      class_lbl <- mapply(function(a, b) {
        sig_a <- abs(a) >= thr; sig_b <- abs(b) >= thr
        if (sig_a && sig_b) {
          if (sign(a) == sign(b)) "concordant" else "discordant"
        } else if (sig_a) {
          sprintf("%s-specific", ds_a)
        } else {
          sprintf("%s-specific", ds_b)
        }
      }, val_a, val_b)
      desc <- vapply(ids, function(id) {
        rows <- all_df[all_df$ID == id, ]
        rows$Description[which.min(rows$p.adjust)]
      }, character(1))

      plot_df <- data.frame(`GO ID` = ids, `GO Term` = desc, a = val_a, b = val_b,
                             Class = class_lbl, check.names = FALSE, stringsAsFactors = FALSE)
      colnames(plot_df)[colnames(plot_df) == "a"] <- ds_a
      colnames(plot_df)[colnames(plot_df) == "b"] <- ds_b

      csv_path <- file.path(output_dir, sprintf("pair_scatter_data_%s_%s_vs_%s_%s.csv", pr, ds_a, ds_b, ont))
      write.csv(plot_df, csv_path, row.names = FALSE)
      message(sprintf("  [pair_scatter] %s (%s vs %s, %s): %d terms (%d concordant, %d discordant)",
                       pr, ds_a, ds_b, ont, nrow(plot_df),
                       sum(class_lbl == "concordant"), sum(class_lbl == "discordant")))

      if (!viz_enabled) next
      tryCatch({
        class_levels <- c("concordant", "discordant", sprintf("%s-specific", ds_a), sprintf("%s-specific", ds_b))
        class_colors <- setNames(c("#2166AC", "#B2182B", "#999999", "#4DAF4A"), class_levels)
        p <- ggplot(plot_df, aes(x = .data[[ds_a]], y = .data[[ds_b]], color = Class)) +
          geom_hline(yintercept = c(-thr, thr), linetype = "dashed", color = "grey75") +
          geom_vline(xintercept = c(-thr, thr), linetype = "dashed", color = "grey75") +
          geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
          geom_point(alpha = 0.75, size = 2) +
          scale_color_manual(values = class_colors, breaks = class_levels, name = NULL) +
          labs(x = sprintf("%s: sign(direction) x -log10(FDR)", ds_a),
               y = sprintf("%s: sign(direction) x -log10(FDR)", ds_b),
               title = sprintf("Cross-Dataset GO Term Comparison: %s (%s)", pr, ont),
               caption = sprintf("dashed lines = FDR %.2g threshold; dotted = y=x", fdr_cutoff)) +
          theme_bw(base_size = 12)
        ggsave(file.path(output_dir, sprintf("pair_scatter_%s_%s_vs_%s_%s.png", pr, ds_a, ds_b, ont)),
               plot = p, width = 7.5, height = 6.5, bg = "white")
      }, error = function(e) message(paste("  [pair_scatter] plot failed:", conditionMessage(e))))
    }
  }
}

# --- 6. 온톨로지별 실행 (A. common / B. flip / C. exclusive / D. mixed + rrvgo + 시각화) ---
run_for_ontology <- function(ont) {
  message(sprintf("\n[18_run_cross_dataset_go_comparison] === Ontology: %s ===", ont))

  parts <- list()
  for (cond in conditions) {
    for (dir in c("up", "down")) {
      d <- load_condition_direction(cond, dir, ont)
      n <- if (is.null(d)) 0L else nrow(d)
      msg <- sprintf("  %s / %s(%s): %d significant terms", cond, dir, ont, n)
      message(msg)
      condition_count_log[[length(condition_count_log) + 1]] <<- msg
      if (n < sparsity_warn_threshold) {
        warn_msg <- sprintf("  [WARN] %s/%s(%s) has only %d significant terms (< %d) — common-term results may be skewed by this sparsity.",
                             cond, dir, ont, n, sparsity_warn_threshold)
        message(warn_msg)
        condition_count_log[[length(condition_count_log) + 1]] <<- warn_msg
      }
      if (!is.null(d)) parts[[paste(cond, dir)]] <- d
    }
  }
  if (length(parts) == 0) {
    message(sprintf("  [18_run_cross_dataset_go_comparison] No significant terms found for any condition/direction in %s — skipping.", ont))
    return(invisible(NULL))
  }
  all_df <- do.call(rbind, parts)
  rownames(all_df) <- NULL

  term_desc <- function(id) {
    rows <- all_df[all_df$ID == id, ]
    rows$Description[which.min(rows$p.adjust)]
  }
  # geneID는 이미 load_condition_direction()에서 SYMBOL로 변환됨 — 종간 대소문자
  # 표기 차이(mouse "Bdnf" vs human "BDNF")를 흡수하려고 toupper()만 추가 적용.
  # 정식 ortholog 매핑은 아니라 참고용 지표.
  gene_set_for <- function(sub_df) {
    if (nrow(sub_df) == 0) return(character(0))
    unique(toupper(unlist(strsplit(sub_df$geneID, "/", fixed = TRUE))))
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
        `Gene Symbols (union)` = paste(gene_set_for(sub), collapse = "/"),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }))
    agg <- agg[order(-agg$`N Conditions`, agg$`Min Adjusted P-value`), ]

    loose <- agg[agg$`N Conditions` >= min_conditions_common, ]
    strict <- loose[!(loose$`GO ID` %in% opp_ids), ]

    write.csv(loose,  file.path(output_dir, sprintf("common_%s_loose_%s.csv", tolower(dir), ont)), row.names = FALSE)
    write.csv(strict, file.path(output_dir, sprintf("common_%s_strict_%s.csv", tolower(dir), ont)), row.names = FALSE)
    message(sprintf("  [common_%s] loose=%d terms, strict=%d terms (>= %d conditions)",
                dir, nrow(loose), nrow(strict), min_conditions_common))
    common_strict_ids[[dir]] <- strict$`GO ID`
    common_loose_ids[[dir]]  <- loose$`GO ID`
  }

  all_flip_ids <- character(0)
  all_excl_ids <- character(0)
  all_mixed_ids <- character(0)
  term_ids_for <- function(conds, dir) unique(all_df$ID[all_df$condition %in% conds & all_df$direction == dir])

  gene_overlap_rows <- list()

  # --- B. find_direction_flip() ---
  for (flip in flips) {
    group_a_name <- flip[[1]]; group_b_name <- flip[[2]]; dir_from <- toupper(flip[[3]])
    dir_to <- setdiff(c("UP", "DOWN"), dir_from)
    group_a <- groups[[group_a_name]]; group_b <- groups[[group_b_name]]

    flip_ids <- setdiff(
      intersect(term_ids_for(group_a, dir_from), term_ids_for(group_b, dir_to)),
      union(term_ids_for(group_a, dir_to), term_ids_for(group_b, dir_from))
    )

    label <- sprintf("flip_%s_%s_to_%s_%s", dir_from, group_a_name, dir_to, group_b_name)
    if (length(flip_ids) == 0) {
      message(sprintf("  [%s] 0 terms.", label))
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
        `Gene Symbols A` = paste(genes_a, collapse = "/"),
        `Gene Symbols B` = paste(genes_b, collapse = "/"),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    flip_df <- do.call(rbind, rows)
    flip_df <- flip_df[order(flip_df$`Min Adj P A`), ]
    write.csv(flip_df, file.path(output_dir, paste0(label, "_", ont, ".csv")), row.names = FALSE)
    message(sprintf("  [%s] %d terms (%d low-overlap flagged, Jaccard < %.1f)",
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
      message(sprintf("  [exclusive:%s] 0 terms.", label))
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
        `Gene Symbols` = paste(gene_set_for(sub), collapse = "/"),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    excl_df <- do.call(rbind, rows)
    excl_df <- excl_df[order(excl_df$`Min Adjusted P-value`), ]
    write.csv(excl_df, file.path(output_dir, paste0("exclusive_", label, "_", ont, ".csv")), row.names = FALSE)
    message(sprintf("  [exclusive:%s] %d terms.", label, nrow(excl_df)))
  }

  # --- D. find_mixed() ---
  for (group_name in names(groups)) {
    group_conds <- groups[[group_name]]
    up_ids <- term_ids_for(group_conds, "UP")
    down_ids <- term_ids_for(group_conds, "DOWN")
    mixed_ids <- intersect(up_ids, down_ids)
    if (length(mixed_ids) == 0) {
      message(sprintf("  [mixed:%s] 0 terms.", group_name))
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
    message(sprintf("  [mixed:%s] %d terms.", group_name, nrow(mixed_df)))
  }

  # --- 신규: pair별 두 데이터셋 직접 산점도 ---
  draw_pair_scatter(all_df, ont)

  # --- rrvgo 의미론적 축약 재적용 (GO만 해당, KEGG는 term이 적어 생략) ---
  semantic_representative_ids <- character(0)
  if (isTRUE(semantic_enabled) && toupper(ont) != "KEGG") {
    for (dir in c("UP", "DOWN")) {
      strict_ids <- common_strict_ids[[dir]]
      if (length(strict_ids) < 2) next
      scores <- vapply(strict_ids, function(id) -log10(max(min(all_df$p.adjust[all_df$ID == id]), 1e-300)), numeric(1))
      names(scores) <- strict_ids
      reduced <- tryCatch({
        suppressPackageStartupMessages(library(rrvgo))
        simMatrix <- calculateSimMatrix(strict_ids, orgdb = primary_organism_db, ont = ont, method = rr_method)
        reduceSimMatrix(simMatrix, scores = scores, threshold = rr_threshold, orgdb = primary_organism_db)
      }, error = function(e) {
        message(sprintf("  [semantic_filter] common_%s_strict(%s) failed: %s", tolower(dir), ont, conditionMessage(e)))
        NULL
      })
      if (is.null(reduced)) next
      out_csv <- file.path(output_dir, sprintf("common_%s_strict_rrvgo_%s.csv", tolower(dir), ont))
      write.csv(reduced[order(reduced$cluster, -reduced$score), ], out_csv, row.names = FALSE)
      rep_ids <- unlist(lapply(split(reduced, reduced$parent), function(g) g$go[which.max(g$score)]))
      semantic_representative_ids <- c(semantic_representative_ids, unname(rep_ids))
      message(sprintf("  [semantic_filter] common_%s_strict(%s): %d terms -> %d representative parent groups",
                  tolower(dir), ont, length(strict_ids), length(unique(reduced$parent))))
    }
    semantic_representative_ids <- unique(semantic_representative_ids)
  }

  # --- 시각화: dot plot / heatmap / UpSet (12와 동일 로직, 축만 conditions=dataset::pair) ---
  if (viz_enabled) {
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
      message(sprintf("  [viz] Heatmap: %d terms -> %s", length(ids), out_name_heatmap))
      tryCatch({
        raw_labels <- vapply(ids, term_desc, character(1))
        trunc_labels <- ifelse(nchar(raw_labels) > 55, paste0(substr(raw_labels, 1, 52), "..."), raw_labels)
        term_labels <- make.unique(trunc_labels)
        mat <- signed_matrix_for(ids)
        rownames(mat) <- term_labels
        dist_rows <- dist(mat)

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
        message(paste("  [viz] heatmap failed:", conditionMessage(e)))
      })
    }

    draw_dotplot <- function(ids, out_name_dot, title) {
      if (length(ids) == 0) return(invisible(NULL))
      message(sprintf("  [viz] Dot plot: %d terms -> %s", length(ids), out_name_dot))
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
        label_w <- 0.085 * max(nchar(as.character(levels(plot_df$`GO Term`))))
        w <- 2.6 + 0.8 * length(conditions) + label_w
        h <- max(5, 0.28 * length(ids) + 2)
        ggsave(file.path(output_dir, out_name_dot), plot = p, width = w, height = h, bg = "white", limitsize = FALSE)
      }, error = function(e) message(paste("  [viz] dot plot failed:", conditionMessage(e))))
    }

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

      best_dir <- vapply(ids_by_sig, function(id) {
        sub <- all_df[all_df$ID == id, ]
        sub$direction[which.min(sub$p.adjust)]
      }, character(1))
      half <- ceiling(dotplot_max_terms / 2)
      up_top   <- head(ids_by_sig[best_dir == "UP"],   half)
      down_top <- head(ids_by_sig[best_dir == "DOWN"], half)
      dot_ids <- ids_by_sig[ids_by_sig %in% c(up_top, down_top)]

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
                sprintf("cross_dataset_dotplot_%s.png", ont), sprintf("cross_dataset_heatmap_%s.png", ont),
                sprintf("cross_dataset_plot_data_%s.csv", ont),
                sprintf("Cross-Dataset GO Dot Plot (%s)", ont))

    if (length(semantic_representative_ids) > 0) {
      render_plot(semantic_representative_ids,
                  sprintf("cross_dataset_dotplot_semantic_%s.png", ont), sprintf("cross_dataset_heatmap_semantic_%s.png", ont),
                  sprintf("cross_dataset_plot_data_semantic_%s.csv", ont),
                  sprintf("Cross-Dataset GO Dot Plot - Semantic Filtered (%s)", ont))
    }

    mixed_ids_unique <- unique(all_mixed_ids)
    if (length(mixed_ids_unique) > 0) {
      render_plot(mixed_ids_unique,
                  sprintf("cross_dataset_dotplot_mixed_%s.png", ont), sprintf("cross_dataset_heatmap_mixed_%s.png", ont),
                  sprintf("cross_dataset_plot_data_mixed_%s.csv", ont),
                  sprintf("Cross-Dataset GO Dot Plot - Mixed Direction (%s)", ont))
    }

    # --- UpSet plot: (dataset::pair)별 유의 term 중첩 구조 (UP/DOWN 각각) ---
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
        message(sprintf("  [viz] UpSet plot saved: upset_%s_%s.png (%d terms)", dir, ont, length(ids)))
      }, error = function(e) message(paste("  [viz] UpSet plot (", dir, ") failed:", conditionMessage(e))))
    }
  }
}

for (ont in go_ontologies) run_for_ontology(ont)

# --- 카테고리별 CSV 취합 -> final_cross_dataset_go_results.xlsx ---
message("\n[18_run_cross_dataset_go_comparison] Building final_cross_dataset_go_results.xlsx...")

read_csv_if_nonempty <- function(path) {
  if (!file.exists(path)) return(NULL)
  raw <- trimws(paste(readLines(path, warn = FALSE), collapse = ""))
  if (raw == "" || raw == '""') return(NULL)
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(d) == 0) return(NULL)
  d
}

flip_labels <- vapply(flips, function(flip) {
  dir_from <- toupper(flip[[3]]); dir_to <- setdiff(c("UP", "DOWN"), dir_from)
  sprintf("flip_%s_%s_to_%s_%s", dir_from, flip[[1]], dir_to, flip[[2]])
}, character(1))
exclusive_labels <- vapply(exclusives, function(spec) paste0("exclusive_", spec$label), character(1))
mixed_labels <- paste0("mixed_", names(groups))
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
used_sheet_names <- character(0)

for (label in category_labels) {
  for (ont in go_ontologies) {
    path <- file.path(output_dir, paste0(label, "_", ont, ".csv"))
    d <- read_csv_if_nonempty(path)
    if (is.null(d)) next
    d <- data.frame(Category = label, Ontology = ont, d, check.names = FALSE, stringsAsFactors = FALSE)

    # Excel 시트명 31자 제한. 데이터셋 라벨이 길면(예: "mouse_rna"/"mouse_atac")
    # label만으로 이미 31자를 넘어 ont 접미사가 통째로 잘려나가 BP/KEGG 시트가
    # 서로 같은 이름이 되는 충돌이 생길 수 있다(실측 확인) — ont 접미사가 항상
    # 살아남도록 label 쪽을 먼저 줄이고, 그래도 겹치면 숫자를 붙여 최종 방어.
    max_label_len <- 31 - nchar(ont) - 1
    sheet_name <- paste0(substr(label, 1, max_label_len), "_", ont)
    if (sheet_name %in% used_sheet_names) {
      suffix <- 2
      while (paste0(substr(sheet_name, 1, 31 - nchar(suffix) - 1), suffix) %in% used_sheet_names) suffix <- suffix + 1
      sheet_name <- paste0(substr(sheet_name, 1, 31 - nchar(suffix) - 1), suffix)
    }
    used_sheet_names <- c(used_sheet_names, sheet_name)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, d, headerStyle = header_style)
    addStyle(wb, sheet_name, header_style, rows = 1, cols = seq_len(ncol(d)), gridExpand = TRUE)
    pval_cols <- grep("P-value|P value|Adj P|Jaccard", colnames(d))
    text_cols <- setdiff(seq_len(ncol(d)), pval_cols)
    if (length(text_cols) > 0) addStyle(wb, sheet_name, text_style, rows = 2:(nrow(d) + 1), cols = text_cols, gridExpand = TRUE)
    if (length(pval_cols) > 0) addStyle(wb, sheet_name, pvalue_style, rows = 2:(nrow(d) + 1), cols = pval_cols, gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(d)), widths = "auto")
    sheet_summary[[sheet_name]] <- nrow(d)
  }
}

if (length(sheet_summary) == 0) {
  addWorksheet(wb, "No Results")
  writeData(wb, "No Results", data.frame(Message = "No cross-dataset GO categories produced results."))
}

info_df <- data.frame(
  Parameter = c("Datasets", "Pairs", "GO ontologies", "FDR cutoff", "Fold enrichment cutoff",
                "min_conditions_common", "Flips", "Exclusives", "Sheets with results"),
  Value = c(paste(dataset_labels, collapse = ", "), paste(pairs, collapse = ", "),
            paste(go_ontologies, collapse = ", "), fdr_cutoff, fe_cutoff, min_conditions_common,
            length(flips), length(exclusives), length(sheet_summary)),
  stringsAsFactors = FALSE
)
addWorksheet(wb, "Analysis_Info")
writeData(wb, "Analysis_Info", info_df)
addStyle(wb, "Analysis_Info", header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "Analysis_Info", cols = 1:2, widths = c(25, 40))

cross_dataset_xlsx <- file.path(output_dir, "final_cross_dataset_go_results.xlsx")
saveWorkbook(wb, cross_dataset_xlsx, overwrite = TRUE)
message(sprintf("[18_run_cross_dataset_go_comparison] final_cross_dataset_go_results.xlsx saved: %s (%d sheets with data)",
                 cross_dataset_xlsx, length(sheet_summary)))

writeLines(condition_count_log, file.path(output_dir, "condition_count_log.txt"))
message("\n[18_run_cross_dataset_go_comparison] Done.")
