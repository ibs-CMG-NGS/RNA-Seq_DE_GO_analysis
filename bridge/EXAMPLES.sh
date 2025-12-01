#!/bin/bash
################################################################################
# Example Usage: Pipeline Bridge Scripts
#
# This script demonstrates how to use the bridge layer to connect
# RNA-Seq_DE_GO_analysis with RNA-Seq_GO_GSEA_analysis
#
# Run this script to see example commands (does not execute them)
################################################################################

cat << 'EOF'
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║                   Pipeline Bridge - Example Usage                            ║
║                                                                              ║
║  Connecting RNA-Seq_DE_GO_analysis → RNA-Seq_GO_GSEA_analysis              ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝

SCENARIO 1: Convert Single Comparison
═══════════════════════════════════════════════════════════════════════════════
You have completed DE analysis and want to run GSEA on H2O2_vs_Control:

  cd bridge
  
  python3 convert_de_to_gsea.py \
    --input ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
    --output ../../RNA-Seq_GO_GSEA_analysis/data/H2O2_vs_Control.xlsx \
    --comparison-name H2O2_vs_Control

Expected Output:
  ✓ ../RNA-Seq_GO_GSEA_analysis/data/H2O2_vs_Control.xlsx
    - Sheet 1: DE_Results (all genes)
    - Sheet 2: Metadata (statistics)
    - Sheet 3: Significant_Only (padj < 0.05)


SCENARIO 2: Batch Convert All Comparisons
═══════════════════════════════════════════════════════════════════════════════
You want to convert all pairwise comparisons in H2O2_Neuron experiment:

  cd bridge
  
  python3 convert_de_to_gsea.py \
    --batch \
    --de-output-dir ../output/H2O2_Neuron \
    --gsea-input-dir ../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline

This will convert all comparisons found in:
  ../output/H2O2_Neuron/pairwise/*/final_de_results.csv

Expected Output:
  ../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/
    ├── H2O2_vs_Control_DE_results.xlsx
    ├── GABA_vs_Control_DE_results.xlsx
    └── ...


SCENARIO 3: Automated Workflow with Wrapper Script
═══════════════════════════════════════════════════════════════════════════════
Use the automation script for a streamlined workflow:

  cd bridge
  
  # Make script executable (first time only)
  chmod +x run_downstream_analysis.sh
  
  # Run conversion for single comparison
  ./run_downstream_analysis.sh H2O2_vs_Control --experiment-dir H2O2_Neuron

The script will:
  1. Check dependencies (Python, converter script, GSEA pipeline)
  2. Convert DE results to Excel format
  3. Display summary and next steps


SCENARIO 4: Full Pipeline Integration
═══════════════════════════════════════════════════════════════════════════════
Complete workflow from DE analysis to GSEA:

  # Step 1: Run DE analysis
  cd /path/to/RNA-Seq_DE_GO_analysis
  snakemake --configfile config_H2O2_Neuron.yml --cores 4

  # Step 2: Convert all results
  cd bridge
  ./run_downstream_analysis.sh --batch --experiment-dir H2O2_Neuron

  # Step 3: Use converted files in GSEA pipeline
  cd ../../RNA-Seq_GO_GSEA_analysis
  # Open notebooks/GO_Pipeline.ipynb or notebooks/GSEA_Pipeline.ipynb
  # Load data from: data/from_de_pipeline/H2O2_vs_Control_DE_results.xlsx


SCENARIO 5: Using edgeR or limma Results
═══════════════════════════════════════════════════════════════════════════════
If you used edgeR or limma instead of DESeq2:

  cd bridge
  
  # For edgeR results
  python3 convert_de_to_gsea.py \
    --source-pipeline edger \
    --input ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
    --output ../../RNA-Seq_GO_GSEA_analysis/data/H2O2_vs_Control.xlsx
  
  # For limma results
  python3 convert_de_to_gsea.py \
    --source-pipeline limma \
    --input ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
    --output ../../RNA-Seq_GO_GSEA_analysis/data/H2O2_vs_Control.xlsx


SCENARIO 6: Testing the Conversion
═══════════════════════════════════════════════════════════════════════════════
Test the converter before processing all your data:

  cd bridge
  
  # Create test directory
  mkdir -p test_output
  
  # Test single file conversion
  python3 convert_de_to_gsea.py \
    --input ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
    --output test_output/test_conversion.xlsx \
    --comparison-name "TEST_H2O2_vs_Control"
  
  # Verify output
  ls -lh test_output/test_conversion.xlsx
  
  # Open in Excel/LibreOffice to inspect
  # Check sheets: DE_Results, Metadata, Significant_Only


SCENARIO 7: Custom Output Directory
═══════════════════════════════════════════════════════════════════════════════
Save converted files to a custom location:

  cd bridge
  
  ./run_downstream_analysis.sh H2O2_vs_Control \
    --experiment-dir H2O2_Neuron \
    --gsea-input-dir /custom/path/to/gsea/input

Or with the Python script:

  python3 convert_de_to_gsea.py \
    --batch \
    --de-output-dir ../output/H2O2_Neuron \
    --gsea-input-dir /custom/path/to/gsea/input


SCENARIO 8: Checking Conversion Results
═══════════════════════════════════════════════════════════════════════════════
After conversion, verify the output:

  # Method 1: Python pandas
  python3 << PYTHON
import pandas as pd
excel_file = '../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/H2O2_vs_Control_DE_results.xlsx'

# Check sheets
xls = pd.ExcelFile(excel_file)
print("Sheets:", xls.sheet_names)

# Read DE results
df = pd.read_excel(excel_file, sheet_name='DE_Results')
print(f"\nTotal genes: {len(df)}")
print(f"Columns: {list(df.columns)}")
print(f"\nFirst few rows:\n{df.head()}")

# Read metadata
meta = pd.read_excel(excel_file, sheet_name='Metadata')
print(f"\nMetadata:\n{meta}")

# Check significant genes
sig = pd.read_excel(excel_file, sheet_name='Significant_Only')
print(f"\nSignificant genes (padj < 0.05): {len(sig)}")
PYTHON

  # Method 2: Command line (requires csvkit or similar)
  in2csv ../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/H2O2_vs_Control_DE_results.xlsx | head


SCENARIO 9: Debugging Conversion Issues
═══════════════════════════════════════════════════════════════════════════════
If conversion fails, enable detailed logging:

  cd bridge
  
  # Run with output redirection
  python3 convert_de_to_gsea.py \
    --input ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
    --output test.xlsx 2>&1 | tee conversion.log
  
  # Check log file
  cat conversion.log
  
  # Common issues:
  # 1. File not found → Check path
  # 2. Missing columns → Check source pipeline parameter
  # 3. Permission denied → Check write permissions


SCENARIO 10: Integration with Snakemake
═══════════════════════════════════════════════════════════════════════════════
Add conversion as a Snakemake rule (optional advanced usage):

Add this rule to your Snakefile:

  rule convert_to_gsea:
      input:
          "output/{experiment}/pairwise/{comparison}/final_de_results.csv"
      output:
          "../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/{comparison}_DE_results.xlsx"
      log:
          "output/{experiment}/logs/convert_{comparison}.log"
      shell:
          """
          python3 bridge/convert_de_to_gsea.py \
            --input {input} \
            --output {output} \
            --comparison-name {wildcards.comparison} \
            2> {log}
          """

Then run:
  snakemake --configfile config_H2O2_Neuron.yml --cores 4


╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║                            Quick Reference                                   ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝

Most Common Commands:
─────────────────────────────────────────────────────────────────────────────

1. Convert single comparison:
   python3 convert_de_to_gsea.py -i INPUT.csv -o OUTPUT.xlsx

2. Batch convert all:
   python3 convert_de_to_gsea.py --batch \
     --de-output-dir ../output/H2O2_Neuron \
     --gsea-input-dir ../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline

3. Automated workflow:
   ./run_downstream_analysis.sh H2O2_vs_Control --experiment-dir H2O2_Neuron


File Locations:
─────────────────────────────────────────────────────────────────────────────
Input:   output/{experiment}/pairwise/{comparison}/final_de_results.csv
Output:  ../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/{comparison}_DE_results.xlsx


Need Help?
─────────────────────────────────────────────────────────────────────────────
python3 convert_de_to_gsea.py --help
./run_downstream_analysis.sh --help
cat bridge/README.md


╔══════════════════════════════════════════════════════════════════════════════╗
║  For detailed documentation, see: bridge/README.md                           ║
╚══════════════════════════════════════════════════════════════════════════════╝

EOF
