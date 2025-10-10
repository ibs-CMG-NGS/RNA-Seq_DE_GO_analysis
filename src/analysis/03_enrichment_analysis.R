# 파일 경로: src/analysis/03_enrichment_analysis.R

# 필요한 라이브러리 로드
library(clusterProfiler)
library(enrichplot)
library(here)
library(ggplot2)

# --- 1. 셋업 ---
species_info <- config$databases[[config$species]]
organism_db_name <- species_info$organism_db
kegg_organism <- species_info$kegg_code

if (!require(organism_db_name, character.only = TRUE)) {
  stop(paste("Genome package", organism_db_name, "is not installed."))
}
organism_db <- get(organism_db_name)
dp_aes <- config$plot_aesthetics$dotplot


# DE 분석 결과를 불러옵니다. 
res_path <- file.path(output_path, "final_de_results.csv")
res <- read.csv(res_path, row.names = 1) # row.names=1을 유지하여 유전자 ID를 행 이름으로 사용


# --- 3. 유전자 목록 및 Ontology 조합에 따른 반복 분석 ---
cat("Starting enrichment analysis based on config settings...\n")

# config에서 분석할 유전자 목록 그룹을 가져옴
for (gene_set in config$enrichment$gene_lists) {
  
  cat(paste("\n--- Preparing gene list for:", gene_set, "regulated genes ---\n"))
  
  # 유전자 목록 선택
  significant_genes <- subset(res, !is.na(padj) & padj < config$de_analysis$padj_cutoff) # 수정
  
  gene_list <- if (gene_set == "up") {
    subset(significant_genes, log2FoldChange > config$de_analysis$log2fc_cutoff) # 수정
} else if (gene_set == "down") {
    subset(significant_genes, log2FoldChange < -config$de_analysis$log2fc_cutoff) # 수정
  } else { # "total"
    significant_genes
  }
  
  # 분석할 유전자가 없으면 다음 목록으로 넘어감
  if (nrow(gene_list) == 0) {
    cat(paste("No significant genes found for the '", gene_set, "' set. Skipping.\n"))
    next
  }

  # ENSEMBL ID를 ENTREZ ID로 변환
  entrez_ids <- mapIds(organism_db,
                       keys = rownames(gene_list),
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)

  if (length(entrez_ids) == 0) {
    cat(paste("No Entrez IDs could be mapped for the '", gene_set, "' set. Skipping.\n"))
    next
  }

  # --- GO Enrichment Analysis (내부 반복문) ---
  for (ont in config$enrichment$go_ontologies) {
    cat(paste("Running GO analysis for", gene_set, "genes, Ontology:", ont, "\n"))
    
    go_results <- enrichGO(gene = entrez_ids, OrgDb = organism_db, keyType = 'ENTREZID', ont = ont,
                           pAdjustMethod = "BH", pvalueCutoff = config$enrichment$pvalue_cutoff,
                           qvalueCutoff = config$enrichment$qvalue_cutoff)
    
    # 동적 파일 이름 생성
    out_csv <- paste0("go_enrichment_", gene_set, "_", ont, ".csv")
    out_plot <- paste0("go_dotplot_", gene_set, "_", ont, ".png")
    
    write.csv(as.data.frame(go_results), file.path(output_path, out_csv))
    
    if (nrow(go_results) > 0) {
      # config 설정에 따라 x축 변수를 동적으로 선택
      x_var <- dp_aes$x_axis_variable
      
      # y축 정렬을 위해 Description을 factor로 변환
      plot_df <- as.data.frame(go_results)
      plot_df$Description <- factor(plot_df$Description, levels = rev(unique(plot_df$Description[order(plot_df[[x_var]])])))
      
      # color aesthetic을 -log10(p.adjust)로 변경
      # ggplot 코드를 최종 형태로 업그레이드
      go_dotplot <- ggplot(plot_df, aes_string(x = x_var, y = "Description", 
                                               color = "-log10(p.adjust)", size = "Count")) +
        geom_point() +
        scale_color_gradient(low = "blue", high = "red") +
        labs(
          title = paste("GO Enrichment -", ont, "(", gene_set, "regulated )"),
          x = x_var, # 동적으로 x축 라벨 설정
          y = "GO Term",
          color = "-log10(p.adjust)",
          size = "Gene Count"
        ) +
        theme_minimal(base_size = 14) +
        theme(axis.text.y = element_text(size = 10))
      ggsave(file.path(output_path, out_plot), plot = go_dotplot, width = 10, height = 8)
    }
  }

  # --- KEGG Pathway Analysis ---
  cat(paste("Running KEGG analysis for", gene_set, "genes\n"))
  kegg_results <- enrichKEGG(gene = entrez_ids, organism = kegg_organism, pvalueCutoff = config$enrichment$pvalue_cutoff)
  
  out_csv_kegg <- paste0("kegg_enrichment_", gene_set, ".csv")
  out_plot_kegg <- paste0("kegg_dotplot_", gene_set, ".png")
  
  write.csv(as.data.frame(kegg_results), file.path(output_path, out_csv_kegg))
  
  if (nrow(kegg_results) > 0) {
    # config 설정에 따라 x축 변수를 동적으로 선택
      x_var <- dp_aes$x_axis_variable
      
      # y축 정렬을 위해 Description을 factor로 변환
      plot_df <- as.data.frame(kegg_results)
      plot_df$Description <- factor(plot_df$Description, levels = rev(unique(plot_df$Description[order(plot_df[[x_var]])])))
    
    # color aesthetic을 -log10(p.adjust)로 변경
    kegg_dotplot <- ggplot(plot_df, aes_string(x = x_var, y = "Description", 
                                               color = "-log10(p.adjust)", size = "Count")) +
        geom_point() +
        scale_color_gradient(low = "blue", high = "red") +
        labs(
          title = paste("KEGG Enrichment -", ont, "(", gene_set, "regulated )"),
          x = x_var, # 동적으로 x축 라벨 설정
          y = "KEGG Pathway",
          color = "-log10(p.adjust)",
          size = "Gene Count"
        ) +
        theme_minimal(base_size = 14) +
        theme(axis.text.y = element_text(size = 10))
    ggsave(file.path(output_path, out_plot_kegg), plot = kegg_dotplot, width = 10, height = 8)
  }
}


cat("\nEnrichment analysis pipeline finished successfully! 🚀\n")
