# 파일 경로: src/analysis/10_run_coexpression_modules.R
# Coexpression module 분석 — omnibus test로 걸러진 유의 유전자 서브셋에 한해
# DEGreport::degPatterns로 발현 패턴 기반 클러스터링 수행.
# (전체 transcriptome 규모의 WGCNA/CEMiTool은 다루지 않음 — 향후 별도 도입 검토)
# + (선택) CMG-SeqViewer용 parquet + staging JSON export
#
# 사용법: Rscript 10_run_coexpression_modules.R [config_path] [omnibus_csv_path] [output_dir]

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(DESeq2)
  library(DEGreport)
  library(openxlsx)
  library(RColorBrewer)
  library(pheatmap)
  library(dplyr)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. 인자 파싱 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript 10_run_coexpression_modules.R [config_path] [omnibus_csv_path] [output_dir]")
}
config_path      <- args[1]
omnibus_csv_path <- args[2]
output_dir       <- args[3]

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. Config 로드 & coexpression_modules 설정 확인 ---
config <- yaml.load_file(config_path)
cm_cfg <- config$de_analysis$coexpression_modules

if (is.null(cm_cfg) || !isTRUE(cm_cfg$enabled)) {
  cat("[10_run_coexpression_modules] coexpression_modules.enabled is not true — skipping.\n")
  quit(save = "no", status = 0)
}

padj_cutoff     <- cm_cfg$padj_cutoff %||% config$de_analysis$padj_cutoff %||% 0.05
min_genes       <- cm_cfg$min_genes %||% 10
min_cluster_size <- cm_cfg$min_cluster_size %||% 5
export_sv       <- isTRUE(cm_cfg$export_seqviewer)
group_var       <- config$de_analysis$group_variable

cat("[10_run_coexpression_modules] Starting coexpression module analysis\n")
cat(paste("  padj_cutoff      :", padj_cutoff, "\n"))
cat(paste("  min_genes        :", min_genes, "\n"))
cat(paste("  min_cluster_size :", min_cluster_size, "\n"))

# --- 3. Omnibus 결과에서 유의 유전자 필터 ---
omnibus_df <- read.csv(omnibus_csv_path, row.names = 1, check.names = FALSE)
sig_gene_ids <- rownames(omnibus_df)[!is.na(omnibus_df$padj) & omnibus_df$padj <= padj_cutoff]
cat(paste("[10_run_coexpression_modules] Significant genes (padj <=", padj_cutoff, "):", length(sig_gene_ids), "\n"))

if (length(sig_gene_ids) < min_genes) {
  cat(paste("[10_run_coexpression_modules] Only", length(sig_gene_ids), "significant genes (< min_genes =",
            min_genes, ") — skipping coexpression clustering.\n"))
  empty_df <- data.frame(gene_id = character(0), module_id = character(0))
  write.csv(empty_df, file.path(output_dir, "coexpression_module_assignments.csv"), row.names = FALSE)
  file.copy(config_path, file.path(output_dir, "config_used.yml"), overwrite = TRUE)
  cat("[10_run_coexpression_modules] Done (skipped — too few significant genes).\n")
  quit(save = "no", status = 0)
}

# --- 4. 데이터 로드 & VST 변환 (유의 유전자만) ---
counts <- read.csv(here(config$count_data_path), row.names = 1, check.names = FALSE)
meta   <- read.csv(here(config$metadata_path),   row.names = 1)

meta_samples      <- rownames(meta)
missing_in_counts <- meta_samples[!meta_samples %in% colnames(counts)]
if (length(missing_in_counts) > 0) {
  stop(paste("Metadata samples not found in count data:", paste(missing_in_counts, collapse = ", ")))
}
counts <- counts[, meta_samples, drop = FALSE]
meta[[group_var]] <- as.factor(meta[[group_var]])

# --- 4b. (선택) include_groups로 클러스터링에 쓸 샘플만 제한 ---
# 유의 유전자 목록(sig_gene_ids)은 그대로 전체 그룹 기준 omnibus 결과를 사용하고,
# VST 계산 및 degPatterns 클러스터링에 들어가는 샘플만 지정한 그룹으로 제한한다.
# null/미설정이면 기존과 동일하게 전체 그룹을 사용한다(하위 호환).
include_groups <- cm_cfg$include_groups
if (!is.null(include_groups)) {
  include_groups <- as.character(include_groups)
  all_groups <- unique(as.character(meta[[group_var]]))
  unknown_groups <- setdiff(include_groups, all_groups)
  if (length(unknown_groups) > 0) {
    stop(paste0("[10_run_coexpression_modules] include_groups에 '", group_var, "' 컬럼에 없는 값이 있습니다: ",
                paste(unknown_groups, collapse = ", ")))
  }
  keep_samples <- rownames(meta)[as.character(meta[[group_var]]) %in% include_groups]
  if (length(keep_samples) < 2) {
    stop("[10_run_coexpression_modules] include_groups 필터 적용 후 샘플이 2개 미만입니다.")
  }
  meta <- meta[keep_samples, , drop = FALSE]
  meta[[group_var]] <- droplevels(meta[[group_var]])
  counts <- counts[, keep_samples, drop = FALSE]
  cat(paste("[10_run_coexpression_modules] include_groups filter:", paste(include_groups, collapse = ", "),
            "->", length(keep_samples), "samples\n"))
}

sig_gene_ids <- intersect(sig_gene_ids, rownames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = as.formula(config$de_analysis$design_formula))
dds <- estimateSizeFactors(dds)
vst_blind <- config$de_analysis$advanced_options$vst_blind %||% FALSE
vsd <- vst(dds, blind = vst_blind)
vst_mat <- assay(vsd)[sig_gene_ids, , drop = FALSE]

# --- 4c. 그룹평균 기준 분산이 0인 유전자 제거 ---
# degPatterns 내부에서 Kendall correlation으로 거리행렬을 만드는데, 그룹 수가 적을 때
# (특히 include_groups로 좁힌 경우) 어떤 유전자의 그룹평균 벡터가 상수(sd=0)이면
# correlation이 NA가 되어 diana()가 "NA values in the dissimilarity matrix not allowed"로
# 죽는다(실제 확인된 실패 사례). 이런 유전자는 그룹 간 패턴 자체가 없으므로 클러스터링
# 입력에서 제외한다.
group_labels_for_var <- as.character(meta[[group_var]])
group_means <- sapply(split(seq_len(ncol(vst_mat)), group_labels_for_var), function(idx) {
  if (length(idx) == 1) vst_mat[, idx] else rowMeans(vst_mat[, idx, drop = FALSE])
})
zero_var_genes <- rownames(vst_mat)[apply(group_means, 1, function(x) stats::sd(x) == 0)]
if (length(zero_var_genes) > 0) {
  cat(paste("[10_run_coexpression_modules] Dropping", length(zero_var_genes),
            "genes with zero variance across group means (would break degPatterns clustering)\n"))
  vst_mat <- vst_mat[!(rownames(vst_mat) %in% zero_var_genes), , drop = FALSE]
}

# --- 5. degPatterns로 발현 패턴 클러스터링 ---
# minc를 명시적으로 낮추지 않으면(기본 15) 유의 유전자가 적을 때 클러스터가 조용히 병합/폐기됨
patterns <- degPatterns(vst_mat, metadata = meta, time = group_var,
                         minc = min_cluster_size, plot = FALSE)

module_df <- patterns$df
colnames(module_df)[colnames(module_df) == "genes"] <- "gene_id"
colnames(module_df)[colnames(module_df) == "cluster"] <- "module_id"
module_df <- module_df[, c("gene_id", "module_id")]
n_modules <- length(unique(module_df$module_id))
cat(paste("[10_run_coexpression_modules]", nrow(module_df), "genes assigned to", n_modules, "modules\n"))

# --- 5b. Module 패턴 플롯 (그룹별 평균 z-score +/- SE, 모듈별 facet) ---
# degPatterns의 내부 plot 객체 구조에 의존하지 않고, 이미 계산 가능한 z-score로 직접 구성
z_mat_all <- t(scale(t(vst_mat[module_df$gene_id, , drop = FALSE])))
sample_group <- as.character(meta[colnames(z_mat_all), group_var])

long_rows <- do.call(rbind, lapply(seq_len(nrow(z_mat_all)), function(i) {
  data.frame(gene_id   = rownames(z_mat_all)[i],
             module_id = module_df$module_id[match(rownames(z_mat_all)[i], module_df$gene_id)],
             group     = sample_group,
             z         = z_mat_all[i, ],
             stringsAsFactors = FALSE)
}))
# 그룹 순서: 대조군(pairwise_comparisons의 base 그룹)을 맨 앞에, 이후는
# pairwise_comparisons에 나열된 순서를 따른다(기존의 단순 알파벳순 대신).
# include_groups로 그룹이 제한된 경우 존재하는 그룹만 남기고, pairwise_comparisons에
# 없는 그룹은 알파벳순으로 맨 뒤에 붙인다(빠짐 방지).
comparisons     <- config$de_analysis$pairwise_comparisons
base_groups     <- unique(vapply(comparisons, function(x) x[[2]], character(1)))
compare_groups  <- unique(vapply(comparisons, function(x) x[[1]], character(1)))
preferred_order <- unique(c(base_groups, compare_groups))
present_groups  <- unique(sample_group)
group_levels    <- c(intersect(preferred_order, present_groups),
                      sort(setdiff(present_groups, preferred_order)))
long_rows$group <- factor(long_rows$group, levels = group_levels)

pattern_summary <- long_rows %>%
  group_by(module_id, group) %>%
  summarise(mean_z = mean(z), se_z = sd(z) / sqrt(n()), .groups = "drop")

pattern_plot <- ggplot(pattern_summary, aes(x = group, y = mean_z, group = module_id)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue") +
  geom_errorbar(aes(ymin = mean_z - se_z, ymax = mean_z + se_z), width = 0.15) +
  facet_wrap(~ module_id, scales = "free_y", labeller = label_both) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = group_var, y = "Mean z-score (+/- SE)", title = "Coexpression module patterns")

pattern_plot_path <- file.path(output_dir, "coexpression_pattern_plot.png")
ggsave(pattern_plot_path, plot = pattern_plot, width = 10, height = 8, dpi = 300, bg = "white")
cat(paste("[10_run_coexpression_modules] Pattern plot saved:", pattern_plot_path, "\n"))

# --- 6. Module heatmap (기존 02b_generate_pairwise_qc_plots.R z-score+RdBu 패턴 재사용) ---
annotation_row <- data.frame(module = factor(module_df$module_id), row.names = module_df$gene_id)
heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)

heatmap_path <- file.path(output_dir, "coexpression_module_heatmap.png")
png(heatmap_path, width = 12, height = 10, units = "in", res = 300, bg = "white")
pheatmap(z_mat_all, color = heatmap_colors, annotation_row = annotation_row,
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
dev.off()
cat(paste("[10_run_coexpression_modules] Module heatmap saved:", heatmap_path, "\n"))

# --- 7. gene_symbol 주석 (선택) ---
species     <- config$species %||% "human"
org_db_name <- config$databases[[species]]$organism_db
final_df <- module_df
final_df$padj <- omnibus_df[final_df$gene_id, "padj"]
if (!is.null(org_db_name)) {
  tryCatch({
    suppressPackageStartupMessages(requireNamespace(org_db_name, quietly = TRUE))
    org_db <- get(org_db_name, envir = loadNamespace(org_db_name))
    gene_id_type <- config$gene_id_type %||% "ENSEMBL"
    symbols <- AnnotationDbi::mapIds(org_db, keys = final_df$gene_id, column = "SYMBOL",
                                      keytype = gene_id_type, multiVals = "first")
    symbols[is.na(symbols)] <- names(symbols)[is.na(symbols)]
    final_df <- cbind(gene_symbol = symbols[final_df$gene_id], final_df)
  }, error = function(e) {
    cat(paste("[10_run_coexpression_modules] gene_symbol annotation skipped:", conditionMessage(e), "\n"))
  })
}

# 그룹별로 정렬된 샘플 컬럼 (VST) 추가
sample_order <- unlist(lapply(sort(unique(as.character(meta[[group_var]]))), function(grp)
  sort(rownames(meta)[meta[[group_var]] == grp])))
vst_ordered <- as.data.frame(vst_mat[final_df$gene_id, sample_order, drop = FALSE])
final_df <- cbind(final_df, vst_ordered)

# --- 8. 결과 저장 (CSV + xlsx + config_used.yml) ---
output_csv_path <- file.path(output_dir, "coexpression_module_assignments.csv")
write.csv(final_df, output_csv_path, row.names = FALSE)
cat(paste("[10_run_coexpression_modules] CSV saved:", output_csv_path, "\n"))

if (isTRUE(config$export$export_to_excel)) {
  wb <- createWorkbook()
  addWorksheet(wb, "Module_Assignments")
  writeData(wb, "Module_Assignments", final_df, rowNames = FALSE)
  saveWorkbook(wb, file.path(output_dir, "coexpression_module_assignments.xlsx"), overwrite = TRUE)
  cat("[10_run_coexpression_modules] Excel file saved.\n")
}

file.copy(config_path, file.path(output_dir, "config_used.yml"), overwrite = TRUE)

# ─────────────────────────────────────────────────────────────
# 9. CMG-SeqViewer export (parquet + staging JSON) — export_seqviewer: true 일 때만
# ─────────────────────────────────────────────────────────────
if (export_sv) {
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
    paste(substr(hex,1,8), substr(hex,9,12),
          paste0("4", substr(hex,14,16)),
          paste0(sample(c("8","9","a","b"),1), substr(hex,18,20)),
          substr(hex,21,32), sep = "-")
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

  # 같은 프로젝트 안에서 include_groups로 여러 서브셋(variant)을 따로 export할 수 있도록
  # variant_label이 있으면 alias/parquet slug/staging 파일명에 반영해 서로 겹치지 않게 한다.
  # 미설정이면 기존과 동일한 이름(하위 호환).
  variant_label <- cm_cfg$variant_label
  cm_alias <- if (!is.null(variant_label)) {
    paste(basename(config$output_dir), "Coexpression-Modules", variant_label)
  } else {
    paste(basename(config$output_dir), "Coexpression-Modules")
  }
  cm_info  <- write_parquet_dataset(final_df, cm_alias, datasets_dir)

  cm_entry <- list(
    dataset_id           = cm_info$uid,
    alias                = cm_alias,
    original_filename    = cm_info$filename,
    dataset_type         = "coexpression_module",
    experiment_condition = paste(sort(unique(as.character(meta[[group_var]]))), collapse = " / "),
    organism             = config$species %||% "",
    cell_type            = "",
    tissue               = "",
    timepoint            = "",
    row_count            = nrow(final_df),
    gene_count           = nrow(final_df),
    significant_genes    = nrow(final_df),
    import_date          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    file_path            = cm_info$filename,
    notes                = paste0("DEGreport::degPatterns clustering (minc=", min_cluster_size,
                                    ") on omnibus-significant genes (padj<=", padj_cutoff, "), ",
                                    n_modules, " modules",
                                    if (!is.null(variant_label)) paste0(" [subset: ", variant_label, "]") else ""),
    tags                 = as.list(c("coexpression", "module",
                                      if (!is.null(variant_label)) variant_label else NULL,
                                      sort(unique(as.character(meta[[group_var]])))))
  )

  # 06b_aggregate_seqviewer.R가 "_entries.json"으로 끝나는 파일만 수집하므로
  # variant_label을 접미사가 아니라 접두어 쪽에 넣어야 한다.
  staging_filename <- if (!is.null(variant_label)) {
    paste0("coexpression_modules_", make_alias_slug(variant_label), "_entries.json")
  } else {
    "coexpression_modules_entries.json"
  }
  staging_path <- file.path(staging_dir, staging_filename)
  write_json(list(cm_entry), staging_path, pretty = TRUE, auto_unbox = TRUE)
  cat(paste("[10_run_coexpression_modules] Seqviewer parquet saved:", cm_info$filename, "\n"))
  cat(paste("[10_run_coexpression_modules] Seqviewer staging JSON saved:", staging_path, "\n"))
}

cat("[10_run_coexpression_modules] Done.\n")
