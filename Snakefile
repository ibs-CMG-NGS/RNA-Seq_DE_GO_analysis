# Snakemake 워크플로우 예제
# RNA-Seq 분석 파이프라인을 Snakemake로 실행하는 예제입니다.
# 
# 사용법:
#   snakemake --cores 4                    # 전체 파이프라인 실행
#   snakemake --cores 4 de_analysis        # DE 분석만 실행
#   snakemake --cores 4 enrichment         # Enrichment 분석까지만 실행

import os
from datetime import datetime

# 설정
CONFIGFILE = "config.yml"
TIMESTAMP = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
OUTPUT_DIR = f"output/{TIMESTAMP}"

# 전체 워크플로우 - 모든 단계의 출력물
rule all:
    input:
        f"{OUTPUT_DIR}/final_de_results.csv",
        f"{OUTPUT_DIR}/pca_plot.png",
        f"{OUTPUT_DIR}/volcano_plot.png",
        # enrichment 파일들은 config에 따라 달라질 수 있으므로 동적으로 처리
        expand(f"{OUTPUT_DIR}/go_enrichment_{{geneset}}_BP.csv", 
               geneset=["total", "up", "down"])

# Step 1: Differential Expression 분석
rule de_analysis:
    input:
        config=CONFIGFILE,
        data="data/raw/airway_scaledcounts.csv",
        metadata="data/raw/airway_metadata.csv"
    output:
        results=f"{OUTPUT_DIR}/final_de_results.csv",
        config_copy=f"{OUTPUT_DIR}/config_used.yml"
    log:
        f"logs/{TIMESTAMP}_de_analysis.log"
    shell:
        """
        Rscript run_pipeline.R \
            --config {input.config} \
            --output {OUTPUT_DIR} \
            --step de \
            2>&1 | tee {log}
        """

# Step 2: 시각화 플롯 생성
rule generate_plots:
    input:
        results=f"{OUTPUT_DIR}/final_de_results.csv",
        config=CONFIGFILE
    output:
        pca=f"{OUTPUT_DIR}/pca_plot.png",
        volcano=f"{OUTPUT_DIR}/volcano_plot.png"
    log:
        f"logs/{TIMESTAMP}_plots.log"
    shell:
        """
        Rscript run_pipeline.R \
            --config {input.config} \
            --output {OUTPUT_DIR} \
            --step plots \
            2>&1 | tee {log}
        """

# Step 3: Enrichment 분석
rule enrichment:
    input:
        results=f"{OUTPUT_DIR}/final_de_results.csv",
        config=CONFIGFILE
    output:
        go_total=f"{OUTPUT_DIR}/go_enrichment_total_BP.csv",
        go_up=f"{OUTPUT_DIR}/go_enrichment_up_BP.csv",
        go_down=f"{OUTPUT_DIR}/go_enrichment_down_BP.csv"
    log:
        f"logs/{TIMESTAMP}_enrichment.log"
    shell:
        """
        Rscript run_pipeline.R \
            --config {input.config} \
            --output {OUTPUT_DIR} \
            --step enrichment \
            2>&1 | tee {log}
        """

# Step 4: GO 플롯 생성
rule go_plots:
    input:
        go_results=expand(f"{OUTPUT_DIR}/go_enrichment_{{geneset}}_BP.csv",
                         geneset=["total", "up", "down"]),
        config=CONFIGFILE
    output:
        f"{OUTPUT_DIR}/go_barplot_total.png",
        f"{OUTPUT_DIR}/go_barplot_up.png",
        f"{OUTPUT_DIR}/go_barplot_down.png"
    log:
        f"logs/{TIMESTAMP}_go_plots.log"
    shell:
        """
        Rscript run_pipeline.R \
            --config {input.config} \
            --output {OUTPUT_DIR} \
            --step go_plots \
            2>&1 | tee {log}
        """

# 전체 파이프라인을 한 번에 실행하는 간단한 규칙
rule run_full_pipeline:
    input:
        config=CONFIGFILE,
        data="data/raw/airway_scaledcounts.csv",
        metadata="data/raw/airway_metadata.csv"
    output:
        f"{OUTPUT_DIR}/pipeline_complete.txt"
    log:
        f"logs/{TIMESTAMP}_full_pipeline.log"
    shell:
        """
        Rscript run_pipeline.R \
            --config {input.config} \
            --output {OUTPUT_DIR} \
            --step all \
            2>&1 | tee {log} && \
        echo "Pipeline completed at $(date)" > {output}
        """

# 클린업 규칙
rule clean:
    shell:
        "rm -rf output/*/ logs/*.log"
