# Snakemake Workflow for Multi-Group DE & GO Analysis

import yaml
from pathlib import Path

# --- 1. Load Configuration ---
configfile: "config.yml"

# 고정된 출력 디렉토리 (예: "results" 또는 "output")
OUTPUT_DIR = Path(config["output_dir"])
# R 분석 환경 이름
R_ENV_NAME = "rna-seq-de-go-analysis" 

# --- 2. Helper Function: Get all comparison pair strings ---
# config.yml에서 ["TreatA", "Control"] 쌍을 읽어 "TreatA_vs_Control" 문자열 리스트로 변환
def get_pairs(config):
    pairs = []
    if "pairwise_comparisons" in config["de_analysis"]:
        for pair_list in config["de_analysis"]["pairwise_comparisons"]:
            compare, base = pair_list
            pairs.append(f"{compare}_vs_{base}")
    return pairs

PAIRS = get_pairs(config) # 예: ["treated_vs_control"]

# --- 3. Target Rule: Define all final outputs ---
rule all:
    input:
        # 1a. Omnibus test result (if requested)
        expand(OUTPUT_DIR / "omnibus_test_results.csv", allow_missing=True) if config["de_analysis"]["run_omnibus_test"] else [],
        
        # 1b. Global PCA Plot (runs once)
        OUTPUT_DIR / "global_pca_plot.png",

        # 2. All Pairwise results
        # 각 1:1 비교 쌍에 대해 DE 분석, Volcano, Enrichment, Bar plot이 모두 생성되어야 함
        expand(OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv", pair=PAIRS),
        expand(OUTPUT_DIR / "pairwise/{pair}/volcano_plot.png", pair=PAIRS),
        expand(OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag", pair=PAIRS), # 03번 스크립트 완료 플래그
        expand(OUTPUT_DIR / "pairwise/{pair}/.go_barplots_done.flag", pair=PAIRS), # 04번 스크립트 완료 플래그

# --- 4. Analysis Rules ---

# [신규] Rule 1a: Run Omnibus Test (전체 그룹 비교)
rule run_omnibus_test:
    input:
        script = "src/analysis/01a_run_omnibus_test.R", # ★ 신규 스크립트
        config_file = "config.yml",
        counts = config["count_data_path"],
        meta = config["metadata_path"]
    output:
        csv = OUTPUT_DIR / "omnibus_test_results.csv"
    log:
        OUTPUT_DIR / "logs/01a_run_omnibus_test.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {output.csv} > {log} 2>&1"

# [신규] Rule 1b: Run Pairwise DE (각 1:1 비교)
rule run_pairwise_de:
    input:
        script = "src/analysis/01b_run_pairwise_de.R", # ★ 신규 스크립트
        config_file = "config.yml",
        counts = config["count_data_path"],
        meta = config["metadata_path"]
    output:
        # 각 비교 쌍을 위한 하위 폴더 및 결과 파일
        out_dir = directory(OUTPUT_DIR / "pairwise/{pair}"),
        csv = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv",
        config_copy = OUTPUT_DIR / "pairwise/{pair}/config_used.yml"
    params:
        compare = "{pair.split('_vs_')[0]}", # 와일드카드에서 'treated' 추출
        base = "{pair.split('_vs_')[1]}"     # 와일드카드에서 'control' 추출
    log:
        OUTPUT_DIR / "logs/01b_run_pairwise_de_{pair}.log"
    conda:
        R_ENV_NAME
    shell:
        # R 스크립트에 config 경로, 비교군, 기준군, 출력 폴더를 인자로 전달
        "Rscript {input.script} {input.config_file} {params.compare} {params.base} {output.out_dir} > {log} 2>&1"

# [수정] Rule 2a: Generate Global PCA Plot (전체 샘플 1회 실행)
rule generate_global_pca:
    input:
        script = "src/analysis/02_generate_plots.R", # 기존 스크립트 재활용
        config_file = "config.yml",
        counts = config["count_data_path"],
        meta = config["metadata_path"]
    output:
        pca = OUTPUT_DIR / "global_pca_plot.png"
    log:
        OUTPUT_DIR / "logs/02_generate_global_pca.log"
    conda:
        R_ENV_NAME
    shell:
        # R 스크립트에 "pca" 작업만 수행하도록 인자 전달
        "Rscript {input.script} {input.config_file} --task pca --output_file {output.pca} > {log} 2>&1"

# [수정] Rule 2b: Generate Pairwise Volcano Plot (각 1:1 비교)
rule generate_pairwise_volcano:
    input:
        script = "src/analysis/02_generate_plots.R", # 기존 스크립트 재활용
        config_file = "config.yml",
        de_results = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv"
    output:
        volcano = OUTPUT_DIR / "pairwise/{pair}/volcano_plot.png"
    params:
        out_file = OUTPUT_DIR / "pairwise/{pair}/volcano_plot.png"
    log:
        OUTPUT_DIR / "logs/02_generate_volcano_{pair}.log"
    conda:
        R_ENV_NAME
    shell:
        # R 스크립트에 "volcano" 작업만 수행하도록 인자 전달
        "Rscript {input.script} {input.config_file} --task volcano --input_file {input.de_results} --output_file {output.volcano} > {log} 2>&1"

# [수정] Rule 3: Run Pairwise Enrichment (각 1:1 비교)
rule enrichment_analysis:
    input:
        script = "src/analysis/03_enrichment_analysis.R",
        config_file = "config.yml",
        de_results = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv"
    output:
        flag = touch(OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag")
    params:
        input_csv = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv",
        output_dir = OUTPUT_DIR / "pairwise/{pair}"
    log:
        OUTPUT_DIR / "logs/03_enrichment_{pair}.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {params.input_csv} {params.output_dir} > {log} 2>&1 && touch {output.flag}"

# [수정] Rule 4: Generate Pairwise GO Barplots (각 1:1 비교)
rule go_barplots:
    input:
        script = "src/analysis/04_generate_go_plots.R",
        config_file = "config.yml",
        enrichment_flag = OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag",
        # 03번 규칙에서 생성된 CSV들에 의존
        go_csvs = expand(
            OUTPUT_DIR / "pairwise/{pair}/go_enrichment_{geneset}_{ontology}.csv",
            pair=r"{{pair}}", # {pair} 와일드카드는 상위 규칙에서 받음
            geneset=config.get("enrichment", {}).get("gene_lists", []),
            ontology=config.get("enrichment", {}).get("go_ontologies", [])
        )
    output:
        flag = touch(OUTPUT_DIR / "pairwise/{pair}/.go_barplots_done.flag")
    params:
        output_dir = OUTPUT_DIR / "pairwise/{pair}"
    log:
        OUTPUT_DIR / "logs/04_go_barplots_{pair}.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {params.output_dir} > {log} 2>&1 && touch {output.flag}"