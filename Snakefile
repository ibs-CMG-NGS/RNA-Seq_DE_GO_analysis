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

        # GO term clustering 결과 (CMG-SeqViewer Clustered-GO 포맷, term_cluster.enabled 시)
        (expand(OUTPUT_DIR / "pairwise/{pair}/final_go_clustered_results.xlsx", pair=PAIRS)
            if config.get("enrichment", {}).get("term_cluster", {}).get("enabled", True) else []),

        # rrvgo 의미론적 축약 결과 (CMG-SeqViewer Clustered-GO 포맷, rrvgo.enabled 시)
        (expand(OUTPUT_DIR / "pairwise/{pair}/final_go_rrvgo_clustered_results.xlsx", pair=PAIRS)
            if config.get("enrichment", {}).get("rrvgo", {}).get("enabled", True) else []),

        # GO Slim overview (up/down 대칭 bar chart, go_slim.enabled 시)
        (expand(OUTPUT_DIR / "pairwise/{pair}/go_slim_overview_BP.png", pair=PAIRS)
            if config.get("enrichment", {}).get("go_slim", {}).get("enabled", True) else []),

        # Cross-Condition GO 비교 (조건 3개 이상 dose-response/time-course에서 유용, 옵트인)
        ([OUTPUT_DIR / "cross_condition/condition_count_log.txt"]
            if config.get("enrichment", {}).get("cross_condition", {}).get("enabled", False) else []),

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
            else []),

        # 6a. Time-series result (maSigPro — de_analysis.time_series.enabled 시)
        (expand(OUTPUT_DIR / "time_series/time_series_significant_genes.csv", allow_missing=True)
            if config.get("de_analysis", {}).get("time_series", {}).get("enabled", False)
            else []),

        # 6a-1. Time-series 클러스터별 GO/KEGG enrichment (+ term_cluster/rrvgo 요약)
        (expand([OUTPUT_DIR / "time_series/final_go_results.xlsx",
                 OUTPUT_DIR / "time_series/final_go_clustered_results.xlsx",
                 OUTPUT_DIR / "time_series/final_go_rrvgo_clustered_results.xlsx"], allow_missing=True)
            if (config.get("de_analysis", {}).get("time_series", {}).get("enabled", False)
                and config.get("de_analysis", {}).get("time_series", {}).get("enrichment_enabled", True))
            else []),

        # 6b. Coexpression module result (omnibus 유의 유전자 서브셋 — coexpression_modules.enabled 시)
        (expand(OUTPUT_DIR / "coexpression_modules/coexpression_module_assignments.csv", allow_missing=True)
            if (config.get("de_analysis", {}).get("run_omnibus_test", False)
                and config.get("de_analysis", {}).get("coexpression_modules", {}).get("enabled", False))
            else []),

        # 6b-1. Coexpression 모듈별 GO/KEGG enrichment (+ term_cluster/rrvgo 요약)
        (expand([OUTPUT_DIR / "coexpression_modules/final_go_results.xlsx",
                 OUTPUT_DIR / "coexpression_modules/final_go_clustered_results.xlsx",
                 OUTPUT_DIR / "coexpression_modules/final_go_rrvgo_clustered_results.xlsx"], allow_missing=True)
            if (config.get("de_analysis", {}).get("run_omnibus_test", False)
                and config.get("de_analysis", {}).get("coexpression_modules", {}).get("enabled", False)
                and config.get("de_analysis", {}).get("coexpression_modules", {}).get("enrichment_enabled", True))
            else []),

        # 7. Google Drive 업로드 (upload.gdrive: true 시)
        ([OUTPUT_DIR / ".gdrive_upload_done.flag"]
            if config.get("upload", {}).get("gdrive", False)
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

# Rule 1b-1: Time-series 분석 (maSigPro) — de_analysis.time_series.enabled 시에만 rule all에 편입
rule run_masigpro_timeseries:
    input:
        script = "src/analysis/01c_run_masigpro_timeseries.R",
        config_file = CONFIG_FILE,
        counts = lambda wildcards: config["count_data_path"] if "count_data_path" in config else [],
        meta = lambda wildcards: config["metadata_path"] if "metadata_path" in config else []
    output:
        csv = OUTPUT_DIR / "time_series/time_series_significant_genes.csv",
        config_copy = OUTPUT_DIR / "time_series/config_used.yml",
        staging = ([OUTPUT_DIR / "seqviewer/staging/time_series_entries.json"]
            if config.get("de_analysis", {}).get("time_series", {}).get("export_seqviewer", True)
            else [])
    params:
        out_dir = str(OUTPUT_DIR / "time_series")
    log:
        OUTPUT_DIR / "logs/01c_run_masigpro_timeseries.log"
    conda:
        R_ENV_NAME
    shell:
        # 큰 유전자 집합에서 see.genes의 덴드로그램 처리가 깊은 재귀를 유발해
        # "C stack usage ... too close to the limit"로 죽는 사례가 확인되어
        # 실행 전 스택 크기를 늘림.
        "ulimit -s unlimited; Rscript {input.script} {input.config_file} {params.out_dir} > {log} 2>&1"

# Rule 1b-2: Coexpression module 분석 — omnibus 유의 유전자 서브셋에 한해
# DEGreport::degPatterns로 클러스터링 (coexpression_modules.enabled 시에만 rule all에 편입)
rule run_coexpression_modules:
    input:
        omnibus_csv = OUTPUT_DIR / "omnibus_test_results.csv",
        script      = "src/analysis/10_run_coexpression_modules.R",
        config_file = CONFIG_FILE,
        counts = lambda wildcards: config["count_data_path"] if "count_data_path" in config else [],
        meta   = lambda wildcards: config["metadata_path"]   if "metadata_path"   in config else []
    output:
        csv = OUTPUT_DIR / "coexpression_modules/coexpression_module_assignments.csv",
        config_copy = OUTPUT_DIR / "coexpression_modules/config_used.yml",
        staging = ([OUTPUT_DIR / "seqviewer/staging/coexpression_modules_entries.json"]
            if config.get("de_analysis", {}).get("coexpression_modules", {}).get("export_seqviewer", True)
            else [])
    params:
        out_dir = str(OUTPUT_DIR / "coexpression_modules")
    log:
        OUTPUT_DIR / "logs/10_run_coexpression_modules.log"
    conda:
        R_ENV_NAME
    shell:
        # degPatterns의 덴드로그램 처리가 유의 유전자 수가 많을 때 깊은 재귀를 유발해
        # "C stack usage ... too close to the limit"로 죽는 사례가 확인되어
        # 실행 전 스택 크기를 늘림.
        "ulimit -s unlimited; Rscript {input.script} {input.config_file} {input.omnibus_csv} {params.out_dir} > {log} 2>&1"

# Rule 1b-3: Time-series 클러스터별(+전체) GO/KEGG enrichment
# (time_series.enrichment_enabled 시에만 rule all에 편입)
rule run_timeseries_enrichment:
    input:
        script = "src/analysis/11_run_group_enrichment.R",
        config_file = CONFIG_FILE,
        csv = OUTPUT_DIR / "time_series/time_series_significant_genes.csv"
    output:
        xlsx = OUTPUT_DIR / "time_series/final_go_results.xlsx",
        # term_cluster(Jaccard)/rrvgo(의미론적 축약) 그룹별 요약 — 05b/05d와 동일한 컬럼
        # 계약, 클러스터/모듈 + TOTAL 전체를 대상으로 함(11_run_group_enrichment.R 내부에서
        # enrichment.term_cluster.enabled / enrichment.rrvgo.enabled=false여도 placeholder를
        # 항상 만들어서 이 output 선언을 항상 만족시킴)
        clustered_xlsx = OUTPUT_DIR / "time_series/final_go_clustered_results.xlsx",
        rrvgo_clustered_xlsx = OUTPUT_DIR / "time_series/final_go_rrvgo_clustered_results.xlsx"
    params:
        out_dir = str(OUTPUT_DIR / "time_series")
    log:
        OUTPUT_DIR / "logs/11_run_timeseries_enrichment.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {input.csv} cluster_id cluster {params.out_dir} > {log} 2>&1"

# Rule 1b-4: Coexpression 모듈별(+전체) GO/KEGG enrichment
# (coexpression_modules.enrichment_enabled 시에만 rule all에 편입)
rule run_coexpression_enrichment:
    input:
        script = "src/analysis/11_run_group_enrichment.R",
        config_file = CONFIG_FILE,
        csv = OUTPUT_DIR / "coexpression_modules/coexpression_module_assignments.csv"
    output:
        xlsx = OUTPUT_DIR / "coexpression_modules/final_go_results.xlsx",
        clustered_xlsx = OUTPUT_DIR / "coexpression_modules/final_go_clustered_results.xlsx",
        rrvgo_clustered_xlsx = OUTPUT_DIR / "coexpression_modules/final_go_rrvgo_clustered_results.xlsx"
    params:
        out_dir = str(OUTPUT_DIR / "coexpression_modules")
    log:
        OUTPUT_DIR / "logs/11_run_coexpression_enrichment.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {input.csv} module_id module {params.out_dir} > {log} 2>&1"

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


# Rule 5b: GO term clustering 결과를 CMG-SeqViewer Clustered-GO 포맷으로 집계
# (03_enrichment_analysis.R의 term-clustering 단계는 geneset=up/down에서만 동작 — "total" 제외)
rule generate_clustered_go_table:
    input:
        script = "src/analysis/05b_generate_clustered_go_table.R",
        config_file = CONFIG_FILE,
        # go_termcluster_{geneset}_{ontology}.csv는 03_enrichment_analysis.R(go_enrichment rule)의
        # 조건부 부산물이라 정식 output으로 선언돼 있지 않음(term_cluster.enabled=false거나 유의
        # term이 부족하면 생성되지 않을 수 있어 고정 output 선언이 불가능). 대신 go_enrichment가
        # 모든 geneset/ontology에 대해 끝난 뒤에만 세팅되는 이 flag로 실행 순서를 보장한다.
        enrichment_flag = OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag"
    output:
        excel = OUTPUT_DIR / "pairwise/{pair}/final_go_clustered_results.xlsx"
    params:
        compare = lambda wildcards: wildcards.pair.split('_vs_')[0],
        base = lambda wildcards: wildcards.pair.split('_vs_')[1],
        output_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair)
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/05b_generate_clustered_go_table.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {params.compare} {params.base} {params.output_dir} > {log} 2>&1"


# Rule 5d: rrvgo 의미론적 축약 결과를 CMG-SeqViewer Clustered-GO 포맷으로 집계
# (05b의 Jaccard/유전자중복도 클러스터링과는 별개 파일 — docs/RRVGO_CLUSTERED_GO_FORMAT.md 참고)
rule generate_rrvgo_clustered_go_table:
    input:
        script = "src/analysis/05d_generate_rrvgo_clustered_go_table.R",
        config_file = CONFIG_FILE,
        # go_rrvgo_{geneset}_{ontology}.csv도 05b와 동일한 이유로 조건부 부산물이라 고정
        # output 선언이 불가능 — go_enrichment가 끝난 뒤 세팅되는 이 flag로 순서만 보장한다.
        enrichment_flag = OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag"
    output:
        excel = OUTPUT_DIR / "pairwise/{pair}/final_go_rrvgo_clustered_results.xlsx"
    params:
        compare = lambda wildcards: wildcards.pair.split('_vs_')[0],
        base = lambda wildcards: wildcards.pair.split('_vs_')[1],
        output_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair)
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/05d_generate_rrvgo_clustered_go_table.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {params.compare} {params.base} {params.output_dir} > {log} 2>&1"


# Rule 5c: GO Slim rollup(level 3 기본) 결과를 ontology별로 up/down 합쳐 대칭 bar chart로 집계
rule generate_go_slim_overview:
    input:
        script = "src/analysis/05c_generate_go_slim_overview.R",
        config_file = CONFIG_FILE,
        # go_slim_{geneset}_{ontology}.csv도 go_termcluster_*.csv와 동일한 이유로 조건부
        # 부산물이라 정식 output 선언이 불가능 — enrichment_done flag로 순서만 보장한다.
        enrichment_flag = OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag"
    output:
        # BP는 스크립트가 항상 생성을 보장(결과가 없어도 placeholder). CC/MF는 보너스 산출물이라
        # 존재 여부와 무관하게 이 규칙의 output으로는 추적하지 않는다.
        overview_bp = OUTPUT_DIR / "pairwise/{pair}/go_slim_overview_BP.png"
    params:
        compare = lambda wildcards: wildcards.pair.split('_vs_')[0],
        base = lambda wildcards: wildcards.pair.split('_vs_')[1],
        output_dir = lambda wildcards: str(OUTPUT_DIR / "pairwise" / wildcards.pair)
    log:
        OUTPUT_DIR / "pairwise/{pair}/logs/05c_generate_go_slim_overview.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {params.compare} {params.base} {params.output_dir} > {log} 2>&1"


# Rule 5e: Cross-Condition GO 비교 (docs/plan_cross_condition_go.md 참고)
# pairwise 비교("조건")가 3개 이상일 때 각 비교의 GO/KEGG enrichment 결과를 서로 대조해
# 공통 term/방향 역전/조건 특이 term을 찾는다 — enrichment.cross_condition.enabled: true로
# 명시적으로 켜야 하는 옵트인 분석(조건이 2개뿐이면 의미가 없어 기본 false).
rule run_cross_condition_comparison:
    input:
        script = "src/analysis/12_run_cross_condition_comparison.R",
        config_file = CONFIG_FILE,
        # 모든 pairwise 비교의 GO/KEGG enrichment가 끝난 뒤에만 실행되도록 전체 pair의
        # enrichment 완료 flag를 의존성으로 건다(cross-condition은 pair 단위가 아니라
        # 프로젝트 전체를 대상으로 한 번만 실행됨).
        enrichment_flags = expand(OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag", pair=PAIRS)
    output:
        # condition_count_log.txt는 스크립트가 항상 마지막에 생성을 보장하는 유일한 고정
        # 파일이다(그 외 common_*/flip_*/exclusive_* 등은 groups/flips/exclusives 설정에
        # 따라 개수가 달라지는 조건부 산출물이라 정식 output으로 선언할 수 없음).
        log_summary = OUTPUT_DIR / "cross_condition/condition_count_log.txt"
    params:
        output_dir = str(OUTPUT_DIR / "cross_condition")
    log:
        OUTPUT_DIR / "logs/12_run_cross_condition_comparison.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {input.config_file} {params.output_dir} > {log} 2>&1"


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
        time_series = ([OUTPUT_DIR / "time_series/time_series_significant_genes.csv"]
            if config.get("de_analysis", {}).get("time_series", {}).get("enabled", False)
            else []),
        coexpression_modules = ([OUTPUT_DIR / "coexpression_modules/coexpression_module_assignments.csv"]
            if (config.get("de_analysis", {}).get("run_omnibus_test", False)
                and config.get("de_analysis", {}).get("coexpression_modules", {}).get("enabled", False))
            else []),
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
            else [],
        ts_staging = lambda wildcards: [OUTPUT_DIR / "seqviewer/staging/time_series_entries.json"]
            if (config.get("de_analysis", {}).get("time_series", {}).get("enabled", False)
                and config.get("de_analysis", {}).get("time_series", {}).get("export_seqviewer", True))
            else [],
        cm_staging = lambda wildcards: [OUTPUT_DIR / "seqviewer/staging/coexpression_modules_entries.json"]
            if (config.get("de_analysis", {}).get("run_omnibus_test", False)
                and config.get("de_analysis", {}).get("coexpression_modules", {}).get("enabled", False)
                and config.get("de_analysis", {}).get("coexpression_modules", {}).get("export_seqviewer", True))
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


# Rule 10: Upload results to Google Drive via rclone
rule upload_to_gdrive:
    input:
        summary   = OUTPUT_DIR / "summary_report.html",
        methods   = OUTPUT_DIR / "methods_section.md",
        go_tables = expand(OUTPUT_DIR / "pairwise/{pair}/final_go_results.xlsx", pair=PAIRS),
        go_clustered_tables = (expand(OUTPUT_DIR / "pairwise/{pair}/final_go_clustered_results.xlsx", pair=PAIRS)
            if config.get("enrichment", {}).get("term_cluster", {}).get("enabled", True) else []),
        go_rrvgo_clustered_tables = (expand(OUTPUT_DIR / "pairwise/{pair}/final_go_rrvgo_clustered_results.xlsx", pair=PAIRS)
            if config.get("enrichment", {}).get("rrvgo", {}).get("enabled", True) else []),
        go_slim_overviews = (expand(OUTPUT_DIR / "pairwise/{pair}/go_slim_overview_BP.png", pair=PAIRS)
            if config.get("enrichment", {}).get("go_slim", {}).get("enabled", True) else []),
        cross_condition = ([OUTPUT_DIR / "cross_condition/condition_count_log.txt"]
            if config.get("enrichment", {}).get("cross_condition", {}).get("enabled", False) else []),
        # Force ordering: without these, upload_to_gdrive is a DAG sibling of
        # the seqviewer export jobs (not a dependent), so Snakemake can finish
        # (and touch) the upload flag before seqviewer/multi-group artifacts
        # exist. Once touched the flag never re-triggers, silently shipping an
        # incomplete upload.
        seqviewer = ([OUTPUT_DIR / "seqviewer/.seqviewer_done.flag"]
            if config.get("export", {}).get("seqviewer", False) else []),
        multi_group = ([OUTPUT_DIR / "seqviewer/staging/multi_group_entries.json"]
            if (config.get("de_analysis", {}).get("run_omnibus_test", False)
                and config.get("de_analysis", {}).get("multi_group_export", {}).get("enabled", False))
            else []),
        time_series = ([OUTPUT_DIR / "time_series/time_series_significant_genes.csv"]
            if config.get("de_analysis", {}).get("time_series", {}).get("enabled", False)
            else []),
        time_series_enrichment = ([OUTPUT_DIR / "time_series/final_go_results.xlsx",
                                    OUTPUT_DIR / "time_series/final_go_clustered_results.xlsx",
                                    OUTPUT_DIR / "time_series/final_go_rrvgo_clustered_results.xlsx"]
            if (config.get("de_analysis", {}).get("time_series", {}).get("enabled", False)
                and config.get("de_analysis", {}).get("time_series", {}).get("enrichment_enabled", True))
            else []),
        coexpression_modules = ([OUTPUT_DIR / "coexpression_modules/coexpression_module_assignments.csv"]
            if (config.get("de_analysis", {}).get("run_omnibus_test", False)
                and config.get("de_analysis", {}).get("coexpression_modules", {}).get("enabled", False))
            else []),
        coexpression_enrichment = ([OUTPUT_DIR / "coexpression_modules/final_go_results.xlsx",
                                     OUTPUT_DIR / "coexpression_modules/final_go_clustered_results.xlsx",
                                     OUTPUT_DIR / "coexpression_modules/final_go_rrvgo_clustered_results.xlsx"]
            if (config.get("de_analysis", {}).get("run_omnibus_test", False)
                and config.get("de_analysis", {}).get("coexpression_modules", {}).get("enabled", False)
                and config.get("de_analysis", {}).get("coexpression_modules", {}).get("enrichment_enabled", True))
            else []),
    output:
        flag = touch(OUTPUT_DIR / ".gdrive_upload_done.flag")
    params:
        remote      = lambda wildcards: config.get("upload", {}).get("remote", "gd_pargilbong:"),
        dest_folder = lambda wildcards: config.get("upload", {}).get("dest_folder", "RNA-Seq_Results"),
        src_dir     = str(OUTPUT_DIR),
        basename    = OUTPUT_DIR.name
    log:
        OUTPUT_DIR / "logs/10_upload_to_gdrive.log"
    shell:
        "/home/ngs/program/anaconda3/bin/rclone copy {params.src_dir} {params.remote}{params.dest_folder}/{params.basename} "
        "--exclude 'logs/10_upload_to_gdrive.log' "
        "--log-level INFO > {log} 2>&1"
