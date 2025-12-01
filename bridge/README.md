# Pipeline Bridge: DE Analysis → GSEA Analysis

This directory contains bridge scripts that connect the **RNA-Seq_DE_GO_analysis** pipeline with the **RNA-Seq_GO_GSEA_analysis** pipeline.

## 📋 Overview

The bridge layer enables seamless data flow between two independent RNA-seq analysis pipelines:

```
RNA-Seq_DE_GO_analysis         Bridge Layer              RNA-Seq_GO_GSEA_analysis
(R-based DESeq2/edgeR)    →    Conversion Scripts   →    (Python-based GO/GSEA)
```

### Workflow

1. **Run DE Analysis**: Execute RNA-Seq_DE_GO_analysis pipeline
2. **Convert Results**: Use bridge scripts to convert CSV → Excel format
3. **Run Advanced Analysis**: Execute RNA-Seq_GO_GSEA_analysis with converted data

---

## 📁 Files in This Directory

| File | Purpose |
|------|---------|
| `convert_de_to_gsea.py` | Python script to convert DE results (CSV) to GSEA format (Excel) |
| `run_downstream_analysis.sh` | Bash automation wrapper to run conversion + GSEA pipeline |
| `README.md` | This documentation file |

---

## 🚀 Quick Start

### 1. Single Comparison Conversion

Convert a single comparison result:

```bash
# Navigate to bridge directory
cd bridge

# Convert single comparison
python3 convert_de_to_gsea.py \
  --input ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
  --output ../../RNA-Seq_GO_GSEA_analysis/data/H2O2_vs_Control.xlsx \
  --comparison-name H2O2_vs_Control
```

### 2. Batch Conversion

Convert all comparisons in an experiment:

```bash
python3 convert_de_to_gsea.py \
  --batch \
  --de-output-dir ../output/H2O2_Neuron \
  --gsea-input-dir ../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline
```

### 3. Automated Workflow

Use the wrapper script for full automation:

```bash
# Make script executable (first time only)
chmod +x run_downstream_analysis.sh

# Run for single comparison
./run_downstream_analysis.sh H2O2_vs_Control \
  --experiment-dir H2O2_Neuron \
  --run-all

# Batch mode
./run_downstream_analysis.sh \
  --batch \
  --experiment-dir H2O2_Neuron
```

---

## 📖 Detailed Usage

### `convert_de_to_gsea.py`

**Purpose**: Convert DE analysis results to GSEA-compatible Excel format

**Input Format** (CSV from DE pipeline):
```
Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
Xkr4,15.8,0.52,0.31,1.68,0.093,0.25
Rp1,8.2,-1.35,0.48,-2.81,0.005,0.032
...
```

**Output Format** (Excel for GSEA pipeline):
- **Sheet 1: DE_Results** - Full DE results with standardized columns
- **Sheet 2: Metadata** - Comparison info, gene counts, statistics
- **Sheet 3: Significant_Only** - Filtered results (padj < 0.05)

**Command Line Options**:

```bash
python3 convert_de_to_gsea.py [OPTIONS]

Single File Mode:
  -i, --input           Input CSV file (final_de_results.csv)
  -o, --output          Output Excel file path
  -n, --comparison-name Name of comparison (auto-detected if omitted)
  -s, --source-pipeline Source pipeline: deseq2, edger, or limma (default: deseq2)

Batch Mode:
  -b, --batch           Enable batch conversion mode
  --de-output-dir       DE analysis output directory
  --gsea-input-dir      Target directory for converted files
```

**Examples**:

```bash
# Example 1: Basic conversion
python3 convert_de_to_gsea.py \
  -i ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
  -o ../../RNA-Seq_GO_GSEA_analysis/data/H2O2_vs_Control.xlsx

# Example 2: Specify comparison name
python3 convert_de_to_gsea.py \
  -i ../output/H2O2_Neuron/pairwise/GABA_vs_Control/final_de_results.csv \
  -o ../../RNA-Seq_GO_GSEA_analysis/data/GABA_vs_Control.xlsx \
  -n "GABA_vs_Control"

# Example 3: Batch convert from edgeR results
python3 convert_de_to_gsea.py \
  --batch \
  --source-pipeline edger \
  --de-output-dir ../output/H2O2_Neuron \
  --gsea-input-dir ../../RNA-Seq_GO_GSEA_analysis/data/edger_results
```

---

### `run_downstream_analysis.sh`

**Purpose**: Automated wrapper to convert results and optionally run GSEA pipeline

**Command Line Options**:

```bash
./run_downstream_analysis.sh <comparison_name> [OPTIONS]

Arguments:
  comparison_name       Name of comparison (e.g., H2O2_vs_Control)

Options:
  --experiment-dir      Experiment directory (e.g., H2O2_Neuron)
  --de-output-dir       Custom DE output directory
  --gsea-input-dir      Custom GSEA input directory
  --run-go              Run GO enrichment analysis
  --run-gsea            Run GSEA analysis
  --run-all             Run both GO and GSEA
  --skip-convert        Skip conversion (use existing Excel)
  --batch               Batch mode: convert all comparisons
  -h, --help            Show help message
```

**Examples**:

```bash
# Example 1: Convert only
./run_downstream_analysis.sh H2O2_vs_Control --experiment-dir H2O2_Neuron

# Example 2: Convert and run GO analysis
./run_downstream_analysis.sh H2O2_vs_Control \
  --experiment-dir H2O2_Neuron \
  --run-go

# Example 3: Convert and run both analyses
./run_downstream_analysis.sh GABA_vs_Control \
  --experiment-dir H2O2_Neuron \
  --run-all

# Example 4: Batch convert all comparisons
./run_downstream_analysis.sh \
  --batch \
  --experiment-dir H2O2_Neuron

# Example 5: Custom directories
./run_downstream_analysis.sh H2O2_vs_Control \
  --de-output-dir /path/to/de/output \
  --gsea-input-dir /path/to/gsea/input
```

---

## 🔧 Configuration

### Default Paths

The scripts use these default paths (modify in `run_downstream_analysis.sh` if needed):

```bash
# DE pipeline root
DE_PIPELINE_ROOT = "../RNA-Seq_DE_GO_analysis"

# GSEA pipeline root
GSEA_PIPELINE_ROOT = "../RNA-Seq_GO_GSEA_analysis"

# DE analysis output
DE_OUTPUT_DIR = "${DE_PIPELINE_ROOT}/output"

# GSEA pipeline input
GSEA_INPUT_DIR = "${GSEA_PIPELINE_ROOT}/data/from_de_pipeline"
```

### Column Mapping

The converter automatically handles different DE pipeline outputs:

| Source | Target (GSEA) | Notes |
|--------|---------------|-------|
| `log2FoldChange` | `log2FoldChange` | DESeq2 format (no change) |
| `logFC` | `log2FoldChange` | edgeR format |
| `pvalue` | `pvalue` | All pipelines |
| `padj` | `padj` | DESeq2/edgeR format |
| `FDR` | `padj` | edgeR format |
| `adj.P.Val` | `padj` | limma format |
| `baseMean` | `baseMean` | DESeq2 format |
| `logCPM` | `baseMean` | edgeR format |
| `AveExpr` | `baseMean` | limma format |

---

## 📊 Output Structure

### Directory Layout

After running the bridge scripts:

```
RNA-Seq_GO_GSEA_analysis/
└── data/
    └── from_de_pipeline/           # Created by bridge scripts
        ├── H2O2_vs_Control_DE_results.xlsx
        ├── GABA_vs_Control_DE_results.xlsx
        └── ...
```

### Excel File Structure

Each converted Excel file contains:

1. **DE_Results** sheet:
   - All genes with DE statistics
   - Standardized column names
   - Added columns: `rank_metric`, `regulation`

2. **Metadata** sheet:
   - Comparison name
   - Total gene count
   - Significant gene count (padj < 0.05)
   - Up/down-regulated counts
   - Source pipeline info
   - Conversion timestamp

3. **Significant_Only** sheet:
   - Filtered results (padj < 0.05)
   - Quick reference for downstream analysis

---

## 🧪 Testing

### Test Single Conversion

```bash
# Test with H2O2_vs_Control data
cd bridge
python3 convert_de_to_gsea.py \
  --input ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
  --output test_output.xlsx \
  --comparison-name "H2O2_vs_Control_TEST"

# Check output
ls -lh test_output.xlsx
```

### Test Batch Conversion

```bash
# Test batch mode
python3 convert_de_to_gsea.py \
  --batch \
  --de-output-dir ../output/H2O2_Neuron \
  --gsea-input-dir test_batch_output

# Check results
ls -lh test_batch_output/
```

### Test Wrapper Script

```bash
# Test wrapper (conversion only)
./run_downstream_analysis.sh H2O2_vs_Control \
  --experiment-dir H2O2_Neuron \
  --gsea-input-dir test_wrapper_output

# Verify output
ls -lh test_wrapper_output/
```

---

## 🔍 Troubleshooting

### Issue 1: Python Not Found

**Error**: `python3: command not found`

**Solution**:
```bash
# Check Python installation
python --version
python3 --version

# If using conda
conda activate your_env
```

### Issue 2: Input File Not Found

**Error**: `Input file not found: .../final_de_results.csv`

**Solution**:
- Verify DE analysis completed successfully
- Check experiment directory name matches
- Use absolute paths if relative paths fail

```bash
# Example with absolute path
python3 convert_de_to_gsea.py \
  --input /full/path/to/final_de_results.csv \
  --output /full/path/to/output.xlsx
```

### Issue 3: Missing Columns

**Error**: `KeyError: 'log2FoldChange'`

**Solution**:
- Specify correct source pipeline: `--source-pipeline edger` or `limma`
- Check input CSV has expected columns

```bash
# For edgeR results
python3 convert_de_to_gsea.py \
  --source-pipeline edger \
  --input results.csv \
  --output output.xlsx
```

### Issue 4: Permission Denied (Shell Script)

**Error**: `Permission denied: ./run_downstream_analysis.sh`

**Solution**:
```bash
# Make script executable
chmod +x run_downstream_analysis.sh

# Or run with bash explicitly
bash run_downstream_analysis.sh H2O2_vs_Control --experiment-dir H2O2_Neuron
```

### Issue 5: Excel File Locked

**Error**: `Failed to save Excel: [Errno 13] Permission denied`

**Solution**:
- Close the Excel file if it's open
- Check write permissions in target directory
- Use a different output filename

---

## 🔗 Integration with Pipelines

### Workflow 1: Manual Integration

```bash
# Step 1: Run DE analysis (in RNA-Seq_DE_GO_analysis)
snakemake --configfile config_H2O2_Neuron.yml --cores 4

# Step 2: Convert results
cd bridge
./run_downstream_analysis.sh --batch --experiment-dir H2O2_Neuron

# Step 3: Run GSEA analysis (in RNA-Seq_GO_GSEA_analysis)
cd ../../RNA-Seq_GO_GSEA_analysis
# Use converted files in data/from_de_pipeline/
```

### Workflow 2: Semi-Automated

Add to your DE pipeline Snakefile:

```python
rule convert_to_gsea:
    input:
        "output/{experiment}/pairwise/{comparison}/final_de_results.csv"
    output:
        "../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/{comparison}_DE_results.xlsx"
    shell:
        """
        python3 bridge/convert_de_to_gsea.py \
          --input {input} \
          --output {output} \
          --comparison-name {wildcards.comparison}
        """
```

### Workflow 3: Fully Automated (Advanced)

Create a master Snakefile that orchestrates both pipelines:

```python
# master_workflow/Snakefile
include: "../RNA-Seq_DE_GO_analysis/Snakefile"
include: "../RNA-Seq_GO_GSEA_analysis/workflow/Snakefile"

rule all:
    input:
        # DE results
        expand("de_output/{comparison}/final_de_results.csv", comparison=COMPARISONS),
        # Converted files
        expand("gsea_input/{comparison}_DE_results.xlsx", comparison=COMPARISONS),
        # GSEA results
        expand("gsea_output/{comparison}/gsea_report.html", comparison=COMPARISONS)
```

---

## 📝 Notes

### Compatibility

- **Python**: Requires Python 3.6+
- **Dependencies**: pandas, openpyxl
- **Shell**: Bash 4.0+ (for automation script)
- **OS**: Linux, macOS, Windows (WSL)

### Data Quality

The converter performs these quality checks:
- Removes rows with NA in critical columns (log2FoldChange, pvalue, padj)
- Sorts by significance (padj)
- Validates column names
- Adds ranking metrics for GSEA

### Performance

- **Single file**: ~1-5 seconds per comparison
- **Batch mode**: Processes comparisons sequentially
- **Memory**: Minimal (~100MB for typical datasets)

---

## 🆘 Getting Help

### Check Logs

The scripts provide detailed logging:

```bash
# Redirect to log file
python3 convert_de_to_gsea.py [options] 2>&1 | tee conversion.log

# Check wrapper script output
./run_downstream_analysis.sh [options] 2>&1 | tee workflow.log
```

### Common Questions

**Q: Can I use this with other DE tools (e.g., limma)?**  
A: Yes! Use `--source-pipeline limma` to handle limma-specific column names.

**Q: What if I have custom column names?**  
A: Modify the `standardize_column_names()` function in `convert_de_to_gsea.py`.

**Q: Can I convert non-CSV files?**  
A: The script expects CSV input. Convert other formats to CSV first.

**Q: How do I handle multiple experiments?**  
A: Run batch conversion separately for each experiment directory.

---

## 📚 References

- [RNA-Seq_DE_GO_analysis README](../README.md)
- [RNA-Seq_GO_GSEA_analysis README](../../RNA-Seq_GO_GSEA_analysis/README.md)
- [DESeq2 Documentation](https://bioconductor.org/packages/DESeq2/)
- [GSEA Documentation](https://www.gsea-msigdb.org/)

---

## 📄 License

This bridge layer inherits the license from the parent pipelines.

---

## 👥 Contributing

To improve the bridge scripts:

1. Test with your data
2. Report issues or suggestions
3. Submit improvements

---

**Last Updated**: 2025-12-01  
**Version**: 1.0.0  
**Status**: Production Ready ✅
