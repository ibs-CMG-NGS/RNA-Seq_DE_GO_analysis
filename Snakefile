import yaml
from pathlib import Path

# --- 1. Load Configuration ---
configfile: "config_H2O2.yml"

OUTPUT_DIR = Path(config["output_dir"])
R_ENV_NAME = "rna-seq-de-go-analysis" 

# --- 2. Helper Function: Get all comparison pair strings ---
def get_pairs(config):
    pairs = []
    if "de_analysis" in config and "pairwise_comparisons" in config["de_analysis"]:
        for pair_list in config["de_analysis"]["pairwise_comparisons"]:
            compare, base = pair_list
            pairs.append(f"{compare}_vs_{base}")
    return pairs

PAIRS = get_pairs(config)

# --- 3. Target Rule: Define all final outputs ---
rule all:
    input:
        # 1a. Omnibus test result (if requested)
        expand(OUTPUT_DIR / "omnibus_test_results.csv", allow_missing=True) if config["de_analysis"]["run_omnibus_test"] else [],
        
        # 1b. Global PCA Plot (runs once)
        OUTPUT_DIR / "global_pca_plot.png",

        # 2. All Pairwise results
        expand(OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv", pair=PAIRS),
        expand(OUTPUT_DIR / "pairwise/{pair}/volcano_plot.png", pair=PAIRS),
        # Enrichment 완료 플래그 (CSV, Dotplot 포함)
        expand(OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag", pair=PAIRS),
        # Barplot 완료 플래그
        expand(OUTPUT_DIR / "pairwise/{pair}/.go_barplots_done.flag", pair=PAIRS)

# --- 4. Analysis Rules ---

# Rule 1a: Omnibus Test
rule run_omnibus_test:
    input:
        script = "src/analysis/01a_run_omnibus_test.R",
        config_file = "config_H2O2.yml",
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

# Rule 1b: Run Pairwise DE
rule run_pairwise_de:
    input:
        script = "src/analysis/01b_run_pairwise_de.R",
        config_file = "config_H2O2.yml",
        counts = config["count_data_path"],
        meta = config["metadata_path"]
    output:
        csv = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv",
        config_copy = OUTPUT_DIR / "pairwise/{pair}/config_used.yml"
    params:
        compare = lambda wildcards: wildcards.pair.split('_vs_')[0],
        base = lambda wildcards: wildcards.pair.split('_vs_')[1],
        out_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair)
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/01b_run_pairwise_de.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {params.compare} {params.base} {params.out_dir} > {log} 2>&1"

# Rule 2a: Generate Global PCA Plot
rule generate_global_pca:
    input:
        script = "src/analysis/02_generate_plots.R",
        config_file = "config_H2O2.yml",
        counts = config["count_data_path"],
        meta = config["metadata_path"]
    output:
        pca = OUTPUT_DIR / "global_pca_plot.png"
    log:
        OUTPUT_DIR / "logs/02_generate_global_pca.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --task pca --output_file {output.pca} > {log} 2>&1"

# Rule 2b: Generate Pairwise Volcano Plot
rule generate_pairwise_volcano:
    input:
        script = "src/analysis/02_generate_plots.R",
        config_file = "config_H2O2.yml",
        de_results = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv"
    output:
        volcano = OUTPUT_DIR / "pairwise/{pair}/volcano_plot.png"
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/02_generate_volcano.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --task volcano --input_file {input.de_results} --output_file {output.volcano} > {log} 2>&1"

# Rule 3a: Pairwise GO Enrichment ({pair}, {geneset}, {ontology})
rule go_enrichment:
    input:
        script = "src/analysis/03_enrichment_analysis.R",
        config_file = "config_H2O2.yml",
        # 각 pair별 DE 결과에 의존
        de_results = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv"
    output:
        go_csv = OUTPUT_DIR / "pairwise/{pair}/go_enrichment_{geneset}_{ontology}.csv",
        go_plot = OUTPUT_DIR / "pairwise/{pair}/go_dotplot_{geneset}_{ontology}.png"
    params:
        output_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair)
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/03a_go_{geneset}_{ontology}.log"
    conda:
        R_ENV_NAME
    shell:
        # 필요한 모든 인자(--task, --geneset, --ontology)를 전달합니다.
        "Rscript {input.script} --config {input.config_file} --input_csv {input.de_results} \
         --output_dir {params.output_dir} --task go \
         --geneset {wildcards.geneset} --ontology {wildcards.ontology} > {log} 2>&1"

# Rule 3b: Pairwise KEGG Enrichment ({pair}, {geneset})
rule kegg_enrichment:
    input:
        script = "src/analysis/03_enrichment_analysis.R",
        config_file = "config_H2O2.yml",
        de_results = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv"
    output:
        kegg_csv = OUTPUT_DIR / "pairwise/{pair}/kegg_enrichment_{geneset}.csv",
        kegg_plot = OUTPUT_DIR / "pairwise/{pair}/kegg_dotplot_{geneset}.png"
    params:
        output_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair)
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/03b_kegg_{geneset}.log"
    conda:
        R_ENV_NAME
    shell:
        # 필요한 모든 인자(--task, --geneset)를 전달합니다. (--ontology는 필요 없음)
        "Rscript {input.script} --config {input.config_file} --input_csv {input.de_results} \
         --output_dir {params.output_dir} --task kegg \
         --geneset {wildcards.geneset} > {log} 2>&1"

# Rule 3c: Mark enrichment as done for this pair
# 개별 GO/KEGG 작업들이 모두 완료되었음을 확인하는 중간 규칙
rule enrichment_done:
    input:
        # 모든 GO 결과물
        expand(
            OUTPUT_DIR / "pairwise/{{pair}}/go_enrichment_{geneset}_{ontology}.csv",
            geneset = config.get("enrichment", {}).get("gene_lists", []),
            ontology = config.get("enrichment", {}).get("go_ontologies", [])
        ),
        # 모든 KEGG 결과물
        expand(
            OUTPUT_DIR / "pairwise/{{pair}}/kegg_enrichment_{geneset}.csv",
            geneset = config.get("enrichment", {}).get("gene_lists", [])
        )
    output:
        flag = touch(OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag")
    shell:
        "echo 'Enrichment analysis complete for {wildcards.pair}'"

# Rule 4: Pairwise GO Barplots
rule go_barplots:
    input:
        script = "src/analysis/04_generate_go_plots.R",
        config_file = "config_H2O2.yml",
        # [수정] 플래그 파일 대신 실제 CSV 파일들을 입력으로 받음
        # enrichment_flag = OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag",
        go_csvs = lambda wildcards: expand(
            OUTPUT_DIR / "pairwise/{pair}/go_enrichment_{geneset}_{ontology}.csv",
            pair=wildcards.pair,
            geneset=config.get("enrichment", {}).get("gene_lists", []),
            ontology=config.get("enrichment", {}).get("go_ontologies", [])
        )
    output:
        flag = touch(OUTPUT_DIR / "pairwise/{pair}/.go_barplots_done.flag")
    params:
        output_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair)
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/04_go_barplots.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --output_dir {params.output_dir} > {log} 2>&1 && touch {output.flag}"