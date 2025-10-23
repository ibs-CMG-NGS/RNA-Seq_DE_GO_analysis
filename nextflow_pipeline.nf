#!/usr/bin/env nextflow

/*
 * RNA-Seq 분석 파이프라인 - Nextflow 예제
 * 
 * 사용법:
 *   nextflow run nextflow_pipeline.nf
 *   nextflow run nextflow_pipeline.nf --config my_config.yml
 *   nextflow run nextflow_pipeline.nf --outdir results/exp1
 */

// 파라미터 정의
params.config = "config.yml"
params.outdir = "output/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}"
params.install_deps = false

// 작업 디렉토리
log.info """
=========================================
RNA-Seq Analysis Pipeline - Nextflow
=========================================
Config file: ${params.config}
Output dir : ${params.outdir}
=========================================
"""

// 의존성 설치 프로세스 (선택적)
process installDependencies {
    when:
    params.install_deps

    output:
    path '.dependencies_installed' into deps_ch

    script:
    """
    Rscript run_pipeline.R --install-deps --step de
    touch .dependencies_installed
    """
}

// DE 분석 프로세스
process differentialExpression {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path config from file(params.config)
    path data from file('data/raw/airway_scaledcounts.csv')
    path metadata from file('data/raw/airway_metadata.csv')

    output:
    path 'final_de_results.csv' into de_results_ch
    path 'config_used.yml'

    script:
    """
    Rscript run_pipeline.R \\
        --config ${config} \\
        --output . \\
        --step de
    """
}

// 플롯 생성 프로세스
process generatePlots {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path de_results from de_results_ch
    path config from file(params.config)

    output:
    path 'pca_plot.png'
    path 'volcano_plot.png'

    script:
    """
    Rscript run_pipeline.R \\
        --config ${config} \\
        --output . \\
        --step plots
    """
}

// Enrichment 분석 프로세스
process enrichmentAnalysis {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path de_results from de_results_ch
    path config from file(params.config)

    output:
    path 'go_enrichment_*.csv' into go_results_ch
    path 'kegg_enrichment_*.csv'
    path 'go_dotplot_*.png'
    path 'kegg_dotplot_*.png'

    script:
    """
    Rscript run_pipeline.R \\
        --config ${config} \\
        --output . \\
        --step enrichment
    """
}

// GO 플롯 생성 프로세스
process generateGOPlots {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path go_results from go_results_ch.collect()
    path config from file(params.config)

    output:
    path 'go_barplot_*.png'

    script:
    """
    Rscript run_pipeline.R \\
        --config ${config} \\
        --output . \\
        --step go_plots
    """
}

// 완료 메시지
workflow.onComplete {
    log.info """
    =========================================
    Pipeline completed!
    Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Results  : ${params.outdir}
    Duration : ${workflow.duration}
    =========================================
    """
}
