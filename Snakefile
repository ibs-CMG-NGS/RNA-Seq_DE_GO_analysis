# Snakemake Workflow for RNA-Seq DE & GO Analysis

import yaml
from pathlib import Path

# --- 1. Load Configuration ---
configfile: "config.yml"

# Output directory from config.yml (should be fixed, e.g., "results")
OUTPUT_DIR = Path(config["output_dir"])

# R analysis conda environment name (from environment.yml)
R_ENV_NAME = "rna-seq-de-go-analysis" # environment.yml의 name과 일치해야 함

# --- 2. Helper Function: Define Final Output Files ---
# Dynamically generates the list of all expected final output files
def get_final_outputs(config, output_dir):
    files = []
    output_dir = Path(output_dir) # Ensure Path object

    # 1. DE analysis results
    files.append(str(output_dir / "final_de_results.csv"))
    files.append(str(output_dir / "config_used.yml"))
    if config.get("export", {}).get("export_to_excel", False):
        files.append(str(output_dir / "final_de_results.xlsx"))

    # 2. Basic plots
    files.append(str(output_dir / "pca_plot.png"))
    files.append(str(output_dir / "volcano_plot.png"))

    # 3. Enrichment analysis results (DotPlots are representative)
    gene_sets = config.get("enrichment", {}).get("gene_lists", [])
    ontologies = config.get("enrichment", {}).get("go_ontologies", [])

    for gs in gene_sets:
        files.append(str(output_dir / f"kegg_dotplot_{gs}.png"))
        for ont in ontologies:
            # Check if ontology results are expected based on namespaces in go_barplot
            if ont in config.get("go_barplot", {}).get("namespaces", []):
                 files.append(str(output_dir / f"go_dotplot_{gs}_{ont}.png"))


    # 4. GO Barplot results
    for gs in gene_sets:
         # Check if barplots are expected based on namespaces
         if config.get("go_barplot", {}).get("namespaces", []):
              files.append(str(output_dir / f"go_barplot_{gs}.png"))

    return files

# --- 3. Target Rule: Define final desired output ---
# [CORRECTED] Depend on the flag files from the final rules
rule all:
    input:
        # Core results from rule run_de_analysis
        OUTPUT_DIR / "final_de_results.csv",
        OUTPUT_DIR / "config_used.yml",
        # Core plots from rule generate_plots
        OUTPUT_DIR / "pca_plot.png",
        OUTPUT_DIR / "volcano_plot.png",
        # Flag file indicating enrichment analysis completion (includes dotplots)
        OUTPUT_DIR / ".enrichment_done.flag",
        # Flag file indicating go barplot completion
        OUTPUT_DIR / ".go_barplots_done.flag"

# --- 4. Analysis Rules ---

# Rule 1: Run Differential Expression Analysis
rule run_de_analysis:
    input:
        script = "src/analysis/01_run_de_analysis.R",
        config_file = "config.yml",
        # Explicitly mention input data files for dependency tracking
        counts = config["count_data_path"],
        meta = config["metadata_path"]
    output:
        csv = OUTPUT_DIR / "final_de_results.csv",
        config_copy = OUTPUT_DIR / "config_used.yml"
    params:
        config_path = "config.yml" # Argument for R script
    log:
        OUTPUT_DIR / "logs/01_run_de_analysis.log"
    conda:
        R_ENV_NAME # Specify the R analysis environment
    shell:
        "Rscript {input.script} {params.config_path} > {log} 2>&1"

# Rule 2: Generate Basic Plots (PCA, Volcano)
rule generate_plots:
    input:
        script = "src/analysis/02_generate_plots.R",
        config_file = "config.yml",
        de_results = OUTPUT_DIR / "final_de_results.csv",
        # Also depends on raw data for PCA re-calculation if needed
        counts = config["count_data_path"],
        meta = config["metadata_path"]
    output:
        pca = OUTPUT_DIR / "pca_plot.png",
        volcano = OUTPUT_DIR / "volcano_plot.png"
    params:
        config_path = "config.yml"
    log:
        OUTPUT_DIR / "logs/02_generate_plots.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {params.config_path} > {log} 2>&1"

# Rule 3: Run Enrichment Analysis (GO & KEGG)
# Creates many output files dynamically, use a flag file.
rule enrichment_analysis:
    input:
        script = "src/analysis/03_enrichment_analysis.R",
        config_file = "config.yml",
        de_results = OUTPUT_DIR / "final_de_results.csv"
    output:
        # Flag file to indicate completion
        flag = touch(OUTPUT_DIR / ".enrichment_done.flag")
    params:
        config_path = "config.yml"
    log:
        OUTPUT_DIR / "logs/03_enrichment_analysis.log"
    conda:
        R_ENV_NAME
    shell:
        # Run script and then create the flag file upon success
        "Rscript {input.script} {params.config_path} > {log} 2>&1 && touch {output.flag}"

# Rule 4: Generate GO Bar Plots
# Depends on the completion flag from the enrichment rule.
rule go_barplots:
    input:
        script = "src/analysis/04_generate_go_plots.R",
        config_file = "config.yml",
        enrichment_flag = OUTPUT_DIR / ".enrichment_done.flag",
        # Needs the enrichment CSVs created by the previous rule
        # Use an expand to list potential input CSVs dynamically
        go_csvs = lambda wildcards: [
            OUTPUT_DIR / f"go_enrichment_{gs}_{ont}.csv"
            for gs in config.get("enrichment", {}).get("gene_lists", [])
            for ont in config.get("enrichment", {}).get("go_ontologies", [])
        ]
    output:
        # Flag file to indicate completion
        flag = touch(OUTPUT_DIR / ".go_barplots_done.flag")
    params:
        config_path = "config.yml"
    log:
        OUTPUT_DIR / "logs/04_generate_go_plots.log"
    conda:
        R_ENV_NAME
    shell:
        "Rscript {input.script} {params.config_path} > {log} 2>&1 && touch {output.flag}"

# --- 5. Optional: Rule to clean results ---
rule clean:
    shell:
        "rm -rf {OUTPUT_DIR}"
