# Bridge Directory - Quick Navigation

Welcome to the Pipeline Bridge! This directory connects **RNA-Seq_DE_GO_analysis** (R-based) with **RNA-Seq_GO_GSEA_analysis** (Python-based) using **Snakemake workflows**.

## 📚 Documentation Index

### 🚀 **Start Here**

#### ⭐ **For Snakemake Users** (권장)
- **[SNAKEMAKE_GUIDE.md](SNAKEMAKE_GUIDE.md)** - Snakemake 통합 워크플로우 가이드
  - Snakemake 기반 파이프라인 연결
  - Rule 설명 및 사용법
  - 병렬 처리 및 최적화
  - 디버깅 및 문제 해결

#### For Python Script Users
- **[README.md](README.md)** - Python 스크립트 사용자 가이드
  - 직접 스크립트 실행
  - 세부 옵션 설명
  - 커스터마이징 방법

### 💡 **Learn by Example**
- **[EXAMPLES.sh](EXAMPLES.sh)** - 10 common usage scenarios
  - Single file conversion
  - Batch processing
  - Automated workflows
  - Testing & validation

### 🎨 **Understand the Flow**
- **[WORKFLOW_DIAGRAM.md](WORKFLOW_DIAGRAM.md)** - Visual architecture
  - System overview diagrams
  - Data flow charts
  - Decision trees
  - Performance characteristics

### ✅ **Implementation Details**
- **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** - Technical reference
  - Completion status
  - Features implemented
  - Testing checklist
  - Maintenance notes

## 🛠️ Core Files

### **Snakefile** ⭐
**Purpose**: Snakemake workflow for automated pipeline integration

**Quick Use**:
```bash
# Convert all comparisons
snakemake -s bridge/Snakefile --cores 1

# Specific comparison
snakemake -s bridge/Snakefile --config comparison=H2O2_vs_Control --cores 1

# Check configuration
snakemake -s bridge/Snakefile show_config
```

**Documentation**: [SNAKEMAKE_GUIDE.md](SNAKEMAKE_GUIDE.md)

---

### **convert_de_to_gsea.py**
**Purpose**: Convert DE analysis CSV results to GSEA-compatible Excel format

**Quick Use**:
```bash
# Single file
python3 convert_de_to_gsea.py -i INPUT.csv -o OUTPUT.xlsx

# Batch mode
python3 convert_de_to_gsea.py --batch \
  --de-output-dir ../output/H2O2_Neuron \
  --gsea-input-dir ../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline
```

**Help**: `python3 convert_de_to_gsea.py --help`

---

### **run_downstream_analysis.sh**
**Purpose**: Automated wrapper for conversion + optional GSEA execution

**Quick Use**:
```bash
# Make executable (first time)
chmod +x run_downstream_analysis.sh

# Run conversion
./run_downstream_analysis.sh H2O2_vs_Control --experiment-dir H2O2_Neuron

# Batch mode
./run_downstream_analysis.sh --batch --experiment-dir H2O2_Neuron
```

**Help**: `./run_downstream_analysis.sh --help`

## 🗺️ Navigation Guide

### If you want to...

#### ✨ **Use Snakemake workflow** (권장)
→ Read: [SNAKEMAKE_GUIDE.md](SNAKEMAKE_GUIDE.md)

#### ✨ **Get started immediately**
→ Run: `snakemake -s bridge/Snakefile --cores 1`

#### 📖 **See practical examples**
→ Run: `cat EXAMPLES.sh` or [View EXAMPLES.sh](EXAMPLES.sh)

#### 🔧 **Understand the architecture**
→ Read: [WORKFLOW_DIAGRAM.md](WORKFLOW_DIAGRAM.md)

#### 🐛 **Fix a problem**
→ Check: [Troubleshooting in README.md](README.md#-troubleshooting)

#### 🧪 **Test the scripts**
→ Follow: [Testing section in IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md#-testing-recommendations)

#### 🔬 **Learn technical details**
→ Read: [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)

#### 📊 **See what's included**
→ Check: [Completion Status](IMPLEMENTATION_SUMMARY.md#-completion-status)

## 🎯 Common Workflows

### Workflow 1: Snakemake User (권장)
```
1. Run DE analysis: snakemake --configfile config.yml --cores 4
2. Convert results: snakemake -s bridge/Snakefile --cores 1
3. Check output: ls -lh ../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/
4. Run GSEA analysis with converted files
```

### Workflow 2: First-Time User
```
1. Read README.md Quick Start section
2. Test with single file conversion
3. Verify Excel output
4. Proceed with your analysis
```

### Workflow 2: Batch Processing
```
1. Review EXAMPLES.sh Scenario 2
2. Run batch conversion command
3. Check output directory
4. Use files in GSEA pipeline
```

### Workflow 3: Troubleshooting
```
1. Check error message
2. Read README.md Troubleshooting section
3. Review EXAMPLES.sh for correct syntax
4. Test with simplified command
```

### Workflow 4: Advanced Integration
```
1. Read WORKFLOW_DIAGRAM.md
2. Review Integration section in README.md
3. Implement custom workflow
4. Test and validate
```

## 📋 Quick Command Reference

### Most Common Commands (Snakemake) ⭐

```bash
# 1. Convert all comparisons (recommended)
snakemake -s bridge/Snakefile --cores 1

# 2. Check configuration
snakemake -s bridge/Snakefile show_config

# 3. List detected comparisons
snakemake -s bridge/Snakefile list_comparisons

# 4. Validate inputs
snakemake -s bridge/Snakefile validate_inputs

# 5. Dry-run (preview)
snakemake -s bridge/Snakefile -n

# 6. Specific comparison
snakemake -s bridge/Snakefile --config comparison=H2O2_vs_Control --cores 1

# 7. Parallel processing
snakemake -s bridge/Snakefile --cores 4
```

### Alternative: Python Script Commands

```bash
# 1. Convert single comparison
cd bridge
python3 convert_de_to_gsea.py \
  -i ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
  -o ../../RNA-Seq_GO_GSEA_analysis/data/H2O2_vs_Control.xlsx

# 2. Batch convert all comparisons
python3 convert_de_to_gsea.py --batch \
  --de-output-dir ../output/H2O2_Neuron \
  --gsea-input-dir ../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline

# 3. Automated workflow
./run_downstream_analysis.sh H2O2_vs_Control --experiment-dir H2O2_Neuron

# 4. Get help
python3 convert_de_to_gsea.py --help
./run_downstream_analysis.sh --help
```

## 🔍 File Descriptions

| File | Type | Lines | Purpose |
|------|------|-------|---------|
| `Snakefile` ⭐ | Snakemake | ~400 | **통합 워크플로우 (권장)** |
| `SNAKEMAKE_GUIDE.md` ⭐ | Markdown | ~500 | **Snakemake 사용 가이드** |
| `convert_de_to_gsea.py` | Python | ~600 | Core conversion logic |
| `run_downstream_analysis.sh` | Bash | ~400 | Automation wrapper (대체 방법) |
| `README.md` | Markdown | ~650 | Python script documentation |
| `EXAMPLES.sh` | Bash | ~450 | Usage examples |
| `WORKFLOW_DIAGRAM.md` | Markdown | ~400 | Visual diagrams |
| `IMPLEMENTATION_SUMMARY.md` | Markdown | ~350 | Technical reference |
| `INDEX.md` | Markdown | This file | Navigation guide |

## 🎓 Learning Path

### Beginner (Snakemake 사용자)
1. **SNAKEMAKE_GUIDE.md** - 빠른 시작 섹션
2. Run: `snakemake -s bridge/Snakefile show_config`
3. Run: `snakemake -s bridge/Snakefile -n` (dry-run)
4. Run: `snakemake -s bridge/Snakefile --cores 1`

### Beginner (Python 스크립트 사용자)
1. **README.md** - Start here, read Quick Start
2. **EXAMPLES.sh** - Run Example 1 (single conversion)
3. **README.md** - Review Output Structure section

### Intermediate
1. **EXAMPLES.sh** - Try Scenarios 2-5
2. **README.md** - Study Configuration section
3. **WORKFLOW_DIAGRAM.md** - Understand data flow

### Advanced
1. **IMPLEMENTATION_SUMMARY.md** - Technical details
2. **README.md** - Integration workflows
3. **WORKFLOW_DIAGRAM.md** - Full architecture

## 💡 Tips & Best Practices

### ✅ Do
- Read README.md Quick Start before first use
- Test with single file before batch processing
- Check example commands in EXAMPLES.sh
- Verify output Excel files after conversion
- Use absolute paths to avoid errors

### ⚠️ Don't
- Run batch conversion without testing first
- Skip dependency checking
- Ignore error messages
- Delete original CSV files after conversion
- Modify column names manually in Excel

## 🆘 Getting Help

### Help Resources (in order)

1. **Quick Questions**: Check [README.md Troubleshooting](README.md#-troubleshooting)
2. **Command Syntax**: See [EXAMPLES.sh](EXAMPLES.sh)
3. **Understanding Flow**: Read [WORKFLOW_DIAGRAM.md](WORKFLOW_DIAGRAM.md)
4. **Technical Issues**: Check [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)
5. **Still Stuck**: Run with `--help` flag

### Command-Line Help

```bash
# Python script help
python3 convert_de_to_gsea.py --help

# Bash script help
./run_downstream_analysis.sh --help

# View examples
cat EXAMPLES.sh
```

## 📊 Success Metrics

After using the bridge scripts, you should have:

✅ Excel file(s) in target directory  
✅ Three sheets per file (DE_Results, Metadata, Significant_Only)  
✅ Standardized column names  
✅ Metadata sheet with summary statistics  
✅ No error messages in terminal  
✅ Files ready for GSEA pipeline  

## 🚀 Next Steps

After successful conversion:

1. **Verify Output**
   - Check Excel file exists
   - Open and review sheets
   - Confirm gene counts in Metadata

2. **Proceed to GSEA**
   - Navigate to `RNA-Seq_GO_GSEA_analysis`
   - Open relevant Jupyter notebook
   - Load converted Excel file

3. **Run Advanced Analysis**
   - Execute GO enrichment
   - Perform GSEA
   - Generate visualizations

## 📞 Contact & Support

For questions or issues:
- Check documentation in this directory
- Review error messages carefully
- Test with example data first
- Report persistent issues with error logs

## 🔖 Version Information

**Current Version**: 1.0.0  
**Last Updated**: 2025-12-01  
**Status**: Production Ready ✅  

**Compatibility**:
- Python 3.6+
- Bash 4.0+
- pandas 1.3+
- openpyxl 3.0+

---

**Happy Analyzing! 🧬✨**

For detailed information, start with **[README.md](README.md)**
