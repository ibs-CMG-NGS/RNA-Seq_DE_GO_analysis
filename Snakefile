import yaml
from pathlib import Path

# --- 1. Load Configuration ---
# ★ Config 파일 경로 (여기서만 수정하면 전체 파이프라인에 적용됨)
# Template: configs/template/config.yml
# User configs: configs/config_ACAS.yml, configs/config_H2O2.yml, etc.
CONFIG_FILE = workflow.configfiles[0] if workflow.configfiles else "configs/config_ACAS.yml"

configfile: CONFIG_FILE

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
        # 1a. Omnibus test result (if requested and count data available)
        expand(OUTPUT_DIR / "omnibus_test_results.csv", allow_missing=True) if (config.get("de_analysis", {}).get("run_omnibus_test", False) and config.get("count_data_path")) else [],
        
        # 1b. Global PCA Plot (if count data available)
        ([OUTPUT_DIR / "global_pca_plot.png"] if config.get("count_data_path") else []),

        # 1c. Global QC Plots (if enabled and count data available)
        expand(OUTPUT_DIR / "qc_plots/.global_qc_done.flag", allow_missing=True) if (config.get("qc_plots", {}).get("generate_global_qc", False) and config.get("count_data_path")) else [],
        
        # 1d. Global QC Report (if enabled and count data available)
        expand(OUTPUT_DIR / "qc_plots/global_qc_report.html", allow_missing=True) if (config.get("qc_plots", {}).get("generate_global_qc", False) and config.get("count_data_path")) else [],

        # 2. All Pairwise results
        expand(OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv", pair=PAIRS),
        expand(OUTPUT_DIR / "pairwise/{pair}/volcano_plot.png", pair=PAIRS),
        # Enrichment 완료 플래그 (CSV, Dotplot 포함)
        expand(OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag", pair=PAIRS),
        # Barplot 완료 플래그
        expand(OUTPUT_DIR / "pairwise/{pair}/.go_barplots_done.flag", pair=PAIRS),
        # GO Summary Table (논문용 통합 Excel 파일)
        expand(OUTPUT_DIR / "pairwise/{pair}/final_go_results.xlsx", pair=PAIRS),
        
        # 2b. Pairwise QC Plots (if enabled)
        expand(OUTPUT_DIR / "pairwise/{pair}/qc_plots/.pairwise_qc_done.flag", pair=PAIRS) if config.get("qc_plots", {}).get("generate_pairwise_qc", False) else [],
        
        # 2c. Pairwise QC Reports (if enabled)
        expand(OUTPUT_DIR / "pairwise/{pair}/qc_plots/pairwise_qc_report.html", pair=PAIRS) if config.get("qc_plots", {}).get("generate_pairwise_qc", False) else [],

        # 3. cmg-seqviewer export (강력 권장, 선택 사항 — export.seqviewer: true 로 활성화)
        [OUTPUT_DIR / "seqviewer/.seqviewer_done.flag"] if config.get("export", {}).get("seqviewer", False) else [],

        # 4. Summary report (모든 pairwise 완료 후 자동 생성)
        OUTPUT_DIR / "summary_report.html",

        # 5. Methods section Markdown
        OUTPUT_DIR / "methods_section.md",

        # 6. Multi-group result (omnibus 확장 — run_omnibus_test + multi_group_export.enabled 시)
        #    CSV는 root에, parquet + staging JSON은 seqviewer/에 저장됨
        ([OUTPUT_DIR / "seqviewer/staging/multi_group_entries.json"]
            if (config.get("de_analysis", {}).get("run_omnibus_test", False)
                and config.get("de_analysis", {}).get("multi_group_export", {}).get("enabled", False))
            else [])

# --- 4. Analysis Rules ---

# Rule 1a: Omnibus Test
rule run_omnibus_test:
    input:
        script = "src/analysis/01a_run_omnibus_test.R",
        config_file = CONFIG_FILE,
        counts = lambda wildcards: config["count_data_path"] if "count_data_path" in config else [],
        meta = lambda wildcards: config["metadata_path"] if "metadata_path" in config else []
    output:
        csv = OUTPUT_DIR / "omnibus_test_results.csv"
    log:
        OUTPUT_DIR / "logs/01a_run_omnibus_test.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {output.csv} > {log} 2>&1"

# Rule 1b: Multi-group result CSV — omnibus 통계 + normalized counts 통합
rule export_multi_group:
    input:
        omnibus_csv = OUTPUT_DIR / "omnibus_test_results.csv",
        script      = "src/analysis/09_export_multi_group.R",
        config_file = CONFIG_FILE,
        counts = lambda wildcards: config["count_data_path"] if "count_data_path" in config else [],
        meta   = lambda wildcards: config["metadata_path"]   if "metadata_path"   in config else []
    output:
        csv     = OUTPUT_DIR / "multi_group_result.csv",
        staging = OUTPUT_DIR / "seqviewer/staging/multi_group_entries.json"
    log:
        OUTPUT_DIR / "logs/09_export_multi_group.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {input.omnibus_csv} {output.csv} > {log} 2>&1"

# Rule 1c: Run Pairwise DE
# Note: If final_de_results.csv already exists, Snakemake will skip this rule
rule run_pairwise_de:
    input:
        script = "src/analysis/01b_run_pairwise_de.R",
        config_file = CONFIG_FILE,
        counts = lambda wildcards: config["count_data_path"] if "count_data_path" in config else [],
        meta = lambda wildcards: config["metadata_path"] if "metadata_path" in config else []
    output:
        csv = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv",
        xlsx = OUTPUT_DIR / "pairwise/{pair}/final_de_results.xlsx",
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
        config_file = CONFIG_FILE,
        counts = lambda wildcards: config["count_data_path"] if "count_data_path" in config else [],
        meta = lambda wildcards: config["metadata_path"] if "metadata_path" in config else []
    output:
        pca = OUTPUT_DIR / "global_pca_plot.png"
    log:
        OUTPUT_DIR / "logs/02_generate_global_pca.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --task pca --output_file {output.pca} > {log} 2>&1"

# Rule 2a-1: Generate Global QC Plots
rule generate_global_qc_plots:
    input:
        script = "src/analysis/02a_generate_qc_plots.R",
        config_file = CONFIG_FILE,
        counts = lambda wildcards: config["count_data_path"] if "count_data_path" in config else [],
        meta = lambda wildcards: config["metadata_path"] if "metadata_path" in config else []
    output:
        flag = touch(OUTPUT_DIR / "qc_plots/.global_qc_done.flag"),
        sample_dist = OUTPUT_DIR / "qc_plots/sample_distance_heatmap.png",
        dispersion = OUTPUT_DIR / "qc_plots/dispersion_plot.png",
        pca = OUTPUT_DIR / "qc_plots/pca_plot.png",
        scree = OUTPUT_DIR / "qc_plots/pca_scree_plot.png",
        boxplot = OUTPUT_DIR / "qc_plots/count_distribution_boxplot.png"
    params:
        output_dir = str(OUTPUT_DIR / "qc_plots")
    log:
        OUTPUT_DIR / "logs/02a_generate_global_qc_plots.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --output_dir {params.output_dir} > {log} 2>&1"

# Rule 2a-2: Generate Global QC Report (HTML)
rule generate_global_qc_report:
    input:
        script = "src/analysis/02c_generate_global_qc_report.R",
        config_file = CONFIG_FILE,
        flag = OUTPUT_DIR / "qc_plots/.global_qc_done.flag",
        plots = [
            OUTPUT_DIR / "qc_plots/sample_distance_heatmap.png",
            OUTPUT_DIR / "qc_plots/dispersion_plot.png",
            OUTPUT_DIR / "qc_plots/pca_plot.png",
            OUTPUT_DIR / "qc_plots/pca_scree_plot.png",
            OUTPUT_DIR / "qc_plots/count_distribution_boxplot.png"
        ]
    output:
        html = OUTPUT_DIR / "qc_plots/global_qc_report.html"
    params:
        qc_plots_dir = str(OUTPUT_DIR / "qc_plots")
    log:
        OUTPUT_DIR / "logs/02c_generate_global_qc_report.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --qc_plots_dir {params.qc_plots_dir} --output_file {output.html} > {log} 2>&1"

# Rule 2b: Generate Pairwise Volcano Plot
rule generate_pairwise_volcano:
    input:
        script = "src/analysis/02_generate_plots.R",
        config_file = CONFIG_FILE,
        de_results = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv"
    output:
        volcano = OUTPUT_DIR / "pairwise/{pair}/volcano_plot.png"
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/02_generate_volcano.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --task volcano --input_file {input.de_results} --output_file {output.volcano} > {log} 2>&1"

# Rule 2b-1: Generate Pairwise QC Plots
rule generate_pairwise_qc_plots:
    input:
        script = "src/analysis/02b_generate_pairwise_qc_plots.R",
        config_file = CONFIG_FILE,
        de_results = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv"
    output:
        flag = touch(OUTPUT_DIR / "pairwise/{pair}/qc_plots/.pairwise_qc_done.flag"),
        ma_plot = OUTPUT_DIR / "pairwise/{pair}/qc_plots/ma_plot.png",
        pval_hist = OUTPUT_DIR / "pairwise/{pair}/qc_plots/pvalue_histogram.png",
        padj_hist = OUTPUT_DIR / "pairwise/{pair}/qc_plots/padj_histogram.png",
        heatmap = OUTPUT_DIR / "pairwise/{pair}/qc_plots/top_genes_heatmap.png",
        fc_dist = OUTPUT_DIR / "pairwise/{pair}/qc_plots/log2fc_distribution.png",
        effect_plot = OUTPUT_DIR / "pairwise/{pair}/qc_plots/effect_size_vs_significance.png"
    params:
        output_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair / "qc_plots"),
        comparison = lambda wildcards: wildcards.pair
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/02b_generate_pairwise_qc_plots.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --comparison {params.comparison} --de_results {input.de_results} --output_dir {params.output_dir} > {log} 2>&1"

# Rule 2b-2: Generate Pairwise QC Report (HTML)
rule generate_pairwise_qc_report:
    input:
        script = "src/analysis/02d_generate_pairwise_qc_report.R",
        config_file = CONFIG_FILE,
        de_results = OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv",
        flag = OUTPUT_DIR / "pairwise/{pair}/qc_plots/.pairwise_qc_done.flag",
        plots = [
            OUTPUT_DIR / "pairwise/{pair}/qc_plots/ma_plot.png",
            OUTPUT_DIR / "pairwise/{pair}/qc_plots/pvalue_histogram.png",
            OUTPUT_DIR / "pairwise/{pair}/qc_plots/padj_histogram.png",
            OUTPUT_DIR / "pairwise/{pair}/qc_plots/top_genes_heatmap.png",
            OUTPUT_DIR / "pairwise/{pair}/qc_plots/log2fc_distribution.png",
            OUTPUT_DIR / "pairwise/{pair}/qc_plots/effect_size_vs_significance.png"
        ]
    output:
        html = OUTPUT_DIR / "pairwise/{pair}/qc_plots/pairwise_qc_report.html"
    params:
        qc_plots_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair / "qc_plots"),
        comparison = lambda wildcards: wildcards.pair
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/02d_generate_pairwise_qc_report.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --comparison {params.comparison} --de_results {input.de_results} --qc_plots_dir {params.qc_plots_dir} --output_file {output.html} > {log} 2>&1"

# Rule 3a: Pairwise GO Enrichment ({pair}, {geneset}, {ontology})
rule go_enrichment:
    input:
        script = "src/analysis/03_enrichment_analysis.R",
        config_file = CONFIG_FILE,
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
        config_file = CONFIG_FILE,
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
        config_file = CONFIG_FILE,
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

# Rule 5: Generate GO Summary Table for Publication
rule generate_go_summary_table:
    input:
        script = "src/analysis/05_generate_go_table.R",
        config_file = CONFIG_FILE,
        enrichment_flag = OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag",
        go_csvs = lambda wildcards: expand(
            OUTPUT_DIR / "pairwise/{pair}/go_enrichment_{geneset}_{ontology}.csv",
            pair=wildcards.pair,
            geneset=config.get("enrichment", {}).get("gene_lists", []),
            ontology=config.get("enrichment", {}).get("go_ontologies", [])
        )
    output:
        excel = OUTPUT_DIR / "pairwise/{pair}/final_go_results.xlsx"
    params:
        compare = lambda wildcards: wildcards.pair.split('_vs_')[0],
        base = lambda wildcards: wildcards.pair.split('_vs_')[1],
        output_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair)
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/05_generate_go_table.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {params.compare} {params.base} {params.output_dir} > {log} 2>&1"


# Rule 6: cmg-seqviewer용 parquet + staging JSON 생성 (per pair)
rule export_seqviewer_pair:
    input:
        script = "src/analysis/06_export_seqviewer.R",
        config_file = CONFIG_FILE,
        de_xlsx = OUTPUT_DIR / "pairwise/{pair}/final_de_results.xlsx",
        enrichment_flag = OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag",
        go_xlsx = OUTPUT_DIR / "pairwise/{pair}/final_go_results.xlsx"
    output:
        flag = touch(OUTPUT_DIR / "pairwise/{pair}/.seqviewer_export_done.flag")
    params:
        compare = lambda wildcards: wildcards.pair.split('_vs_')[0],
        base = lambda wildcards: wildcards.pair.split('_vs_')[1],
        pair_output_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair)
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/06_export_seqviewer.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {params.compare} {params.base} {params.pair_output_dir} > {log} 2>&1"


# Rule 8: Summary Report — 모든 pairwise 완료 후 통합 HTML 리포트 생성
rule generate_summary_report:
    input:
        script      = "src/analysis/07_generate_summary_report.R",
        config_file = CONFIG_FILE,
        de_results  = expand(OUTPUT_DIR / "pairwise/{pair}/final_de_results.csv", pair=PAIRS),
        enrich_done = expand(OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag", pair=PAIRS),
        volcanos    = expand(OUTPUT_DIR / "pairwise/{pair}/volcano_plot.png", pair=PAIRS),
    output:
        html = OUTPUT_DIR / "summary_report.html"
    params:
        output_dir = str(OUTPUT_DIR)
    log:
        OUTPUT_DIR / "logs/07_generate_summary_report.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --output-dir {params.output_dir} --output {output.html} > {log} 2>&1"


# Rule 7: 모든 pair의 staging JSON을 합쳐 metadata.json 생성
rule aggregate_seqviewer:
    input:
        script = "src/analysis/06b_aggregate_seqviewer.R",
        flags = expand(OUTPUT_DIR / "pairwise/{pair}/.seqviewer_export_done.flag", pair=PAIRS),
        mg_staging = lambda wildcards: [OUTPUT_DIR / "seqviewer/staging/multi_group_entries.json"]
            if (config.get("de_analysis", {}).get("run_omnibus_test", False)
                and config.get("de_analysis", {}).get("multi_group_export", {}).get("enabled", False))
            else []
    output:
        flag = touch(OUTPUT_DIR / "seqviewer/.seqviewer_done.flag")
    params:
        seqviewer_dir = str(OUTPUT_DIR / "seqviewer")
    log:
        OUTPUT_DIR / "logs/06b_aggregate_seqviewer.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {params.seqviewer_dir} > {log} 2>&1"


# Rule 9: Methods Section Markdown — 파이프라인 완료 후 분석 방법 문서 자동 생성
rule generate_methods_section:
    input:
        script      = "src/analysis/08_generate_methods_section.R",
        config_file = CONFIG_FILE,
        summary     = OUTPUT_DIR / "summary_report.html",
    output:
        md = OUTPUT_DIR / "methods_section.md"
    params:
        output_dir = str(OUTPUT_DIR)
    log:
        OUTPUT_DIR / "logs/08_generate_methods_section.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} --config {input.config_file} --output-dir {params.output_dir} --output {output.md} > {log} 2>&1"
