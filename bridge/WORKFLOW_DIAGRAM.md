# Pipeline Bridge: Visual Workflow

## 📊 Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                                                                             │
│                     RNA-Seq Analysis Pipeline Ecosystem                     │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘

┌────────────────────────────┐         ┌─────────────────────────────────────┐
│                            │         │                                     │
│  RNA-Seq_DE_GO_analysis    │         │   RNA-Seq_GO_GSEA_analysis          │
│  (R-based)                 │         │   (Python-based)                    │
│                            │         │                                     │
│  ┌──────────────────────┐  │         │  ┌────────────────────────────────┐ │
│  │ DESeq2 / edgeR       │  │         │  │ Advanced GO Enrichment         │ │
│  │ limma-voom           │  │         │  │ - Gene Clustering              │ │
│  │                      │  │         │  │ - Semantic Similarity          │ │
│  │ DE Analysis          │  │         │  │ - Interactive Visualization    │ │
│  └──────────┬───────────┘  │         │  │                                │ │
│             │              │         │  └────────────────────────────────┘ │
│             ▼              │         │                                     │
│  ┌──────────────────────┐  │         │  ┌────────────────────────────────┐ │
│  │ GO Enrichment        │  │         │  │ GSEA Analysis                  │ │
│  │ KEGG Pathway         │  │         │  │ - Rank-based Enrichment        │ │
│  │ Visualization        │  │         │  │ - Leading Edge Analysis        │ │
│  └──────────┬───────────┘  │         │  │ - Enrichment Score Plots       │ │
│             │              │         │  └────────────────────────────────┘ │
│             ▼              │         │                                     │
│  ┌──────────────────────┐  │         │                                     │
│  │ Results:             │  │         │                                     │
│  │ - CSV files          │◄─┼─────────┤ Input Requirements:                 │
│  │ - PNG plots          │  │   🔗    │ - Excel format                     │
│  │ - Excel summaries    │  │ BRIDGE  │ - Standardized columns              │
│  └──────────────────────┘  │  LAYER  │ - Metadata sheet                    │
│                            │         │                                     │
└────────────────────────────┘         └─────────────────────────────────────┘

                                  ▲
                                  │
                    ┌─────────────┴─────────────┐
                    │                           │
                    │   Bridge Scripts          │
                    │   ───────────────         │
                    │   1. convert_de_to_gsea.py│
                    │   2. run_downstream.sh    │
                    │                           │
                    └───────────────────────────┘
```

## 🔄 Data Flow Diagram

```
Step 1: DE Analysis (RNA-Seq_DE_GO_analysis)
═══════════════════════════════════════════════

  Raw Counts + Metadata
         │
         ▼
  ┌──────────────┐
  │   DESeq2     │
  │   edgeR      │ ──► Differential Expression
  │   limma-voom │
  └──────┬───────┘
         │
         ▼
  ┌──────────────────────────────────────────┐
  │  final_de_results.csv                    │
  ├──────────────────────────────────────────┤
  │  Gene | baseMean | log2FC | padj | ...   │
  │  TP53 | 1234.5   | 2.3    | 0.001| ...   │
  │  BRCA1| 567.8    | -1.5   | 0.023| ...   │
  └──────────────────────────────────────────┘


Step 2: Bridge Conversion
══════════════════════════

  final_de_results.csv
         │
         ▼
  ┌────────────────────────┐
  │  convert_de_to_gsea.py │
  │                        │
  │  • Load CSV            │
  │  • Standardize columns │
  │  • Add metadata        │
  │  • Filter quality      │
  │  • Sort by padj        │
  └──────────┬─────────────┘
             │
             ▼
  ┌──────────────────────────────────────┐
  │  H2O2_vs_Control_DE_results.xlsx     │
  ├──────────────────────────────────────┤
  │  Sheet 1: DE_Results (all genes)     │
  │  Sheet 2: Metadata (summary)         │
  │  Sheet 3: Significant_Only (filtered)│
  └──────────────────────────────────────┘


Step 3: Advanced Analysis (RNA-Seq_GO_GSEA_analysis)
═════════════════════════════════════════════════════

  H2O2_vs_Control_DE_results.xlsx
         │
         ▼
  ┌──────────────────┐
  │  GO Pipeline     │  ──► Advanced GO visualization
  │  GSEA Pipeline   │  ──► Gene Set Enrichment Analysis
  └──────────────────┘
```

## 🚦 Workflow Execution Paths

### Path 1: Manual Single Conversion
```
User Action                    Bridge Script                 Output
───────────                    ─────────────                 ──────
Run DE pipeline                                              final_de_results.csv
    │
    ├─► Run converter          convert_de_to_gsea.py ────► Excel file
    │
    └─► Open GSEA notebook     (manual)                  ──► GSEA results
```

### Path 2: Automated Batch Workflow
```
User Action                    Bridge Script                 Output
───────────                    ─────────────                 ──────
Run DE pipeline                                              Multiple CSV files
    │
    └─► Run automation         run_downstream_analysis.sh
            │
            ├─► Check deps     (validation)
            │
            ├─► Batch convert  convert_de_to_gsea.py ────► Multiple Excel files
            │
            └─► (Optional)     Trigger GSEA              ──► GSEA results
```

### Path 3: Integrated Snakemake (Future)
```
snakemake all
    │
    ├─► DE analysis rules      (Existing Snakefile)      ──► CSV results
    │
    ├─► Bridge rule            convert_de_to_gsea.py ────► Excel files
    │
    └─► GSEA rules             (GSEA Snakefile)          ──► Final results
```

## 📁 File Transformation Detail

### Input File Structure
```csv
# final_de_results.csv (DESeq2 output)
Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,symbol
ENSG00000141510,1234.5,2.34,0.45,5.2,2.3e-07,1.5e-05,TP53
ENSG00000012048,567.8,-1.52,0.38,-4.0,6.4e-05,2.1e-03,BRCA1
...
```

### Output File Structure
```
H2O2_vs_Control_DE_results.xlsx

┌─────────────────────────────────────────────────────────┐
│ Sheet 1: DE_Results                                     │
├─────────────────────────────────────────────────────────┤
│ GeneID  │ log2FoldChange │ pvalue    │ padj      │ ...  │
│ TP53    │ 2.34           │ 2.3e-07   │ 1.5e-05   │ ...  │
│ BRCA1   │ -1.52          │ 6.4e-05   │ 2.1e-03   │ ...  │
└─────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────┐
│ Sheet 2: Metadata                                       │
├─────────────────────────────────────────────────────────┤
│ Parameter                     │ Value                   │
│ Comparison                    │ H2O2_vs_Control         │
│ Total Genes                   │ 15234                   │
│ Significant Genes (padj<0.05) │ 1234                    │
│ Up-regulated                  │ 567                     │
│ Down-regulated                │ 667                     │
│ Source Pipeline               │ RNA-Seq_DE_GO_analysis  │
│ Conversion Date               │ 2025-12-01 14:30:15     │
└─────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────┐
│ Sheet 3: Significant_Only (padj < 0.05)                 │
├─────────────────────────────────────────────────────────┤
│ GeneID  │ log2FoldChange │ pvalue    │ padj      │ ...  │
│ TP53    │ 2.34           │ 2.3e-07   │ 1.5e-05   │ ...  │
│ BRCA1   │ -1.52          │ 6.4e-05   │ 2.1e-03   │ ...  │
│ ...     │ ...            │ ...       │ < 0.05    │ ...  │
└─────────────────────────────────────────────────────────┘
```

## 🎯 Decision Tree: Which Path to Use?

```
Do you need to run GSEA analysis on DE results?
│
├─ Yes ─► Is it a single comparison?
│   │
│   ├─ Yes ─► Use: python3 convert_de_to_gsea.py -i INPUT -o OUTPUT
│   │
│   └─ No (multiple) ─► Use: python3 convert_de_to_gsea.py --batch
│
└─ No ─► Continue with standard DE_GO_analysis pipeline only
```

## 🔧 Component Interaction Map

```
┌─────────────────────────────────────────────────────────────────┐
│                        Bridge Components                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────────────┐         ┌─────────────────────────────┐   │
│  │ User Interface   │         │ Core Functions              │   │
│  ├──────────────────┤         ├─────────────────────────────┤   │
│  │ CLI Arguments    │────────►│ load_de_results()           │   │
│  │ --input          │         │ standardize_column_names()  │   │
│  │ --output         │         │ add_required_columns()      │   │
│  │ --batch          │         │ filter_and_sort()           │   │
│  │ --comparison     │         │ save_to_excel()             │   │
│  └──────────────────┘         └─────────────────────────────┘   │
│                                                                 │
│  ┌──────────────────┐         ┌─────────────────────────────┐   │
│  │ Automation       │         │ Output Generation           │   │
│  ├──────────────────┤         ├─────────────────────────────┤   │
│  │ Dependency Check │────────►│ Multi-sheet Excel           │   │
│  │ Path Validation  │         │ Formatted headers           │   │
│  │ Batch Processing │         │ Auto column width           │   │
│  │ Error Handling   │         │ Freeze panes                │   │
│  └──────────────────┘         └─────────────────────────────┘   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## 📊 Performance Characteristics

```
Operation          | Time        | Memory    | Notes
───────────────────┼─────────────┼───────────┼──────────────────────
Single conversion  | 1-5 sec     | ~100 MB   | Typical dataset
Batch (5 files)    | 5-25 sec    | ~100 MB   | Sequential processing
Batch (20 files)   | 20-100 sec  | ~150 MB   | Sequential processing
Large file (50k)   | 10-30 sec   | ~500 MB   | Single comparison
```

## 🎨 User Experience Flow

```
1. Scientist completes DE analysis
   ├─► Has results: final_de_results.csv
   └─► Wants: Advanced GSEA analysis

2. Navigate to bridge directory
   └─► cd bridge

3. Choose workflow:
   ├─► Single: python3 convert_de_to_gsea.py -i INPUT -o OUTPUT
   └─► Batch:  ./run_downstream_analysis.sh --batch --experiment-dir EXPDIR

4. Observe progress:
   ├─► Colored terminal output
   ├─► Progress indicators
   └─► Success/error messages

5. Verify results:
   ├─► Check Excel file created
   ├─► Review metadata sheet
   └─► Confirm gene counts

6. Proceed to GSEA:
   └─► Open RNA-Seq_GO_GSEA_analysis notebooks
```

## 🔐 Error Handling Flow

```
Input Validation
    │
    ├─► File exists? ──No──► Error: File not found
    │                          └─► Exit with message
    └─► Yes
          │
          ├─► Valid CSV? ──No──► Error: Invalid format
          │                        └─► Exit with message
          └─► Yes
                │
                ├─► Required columns? ──No──► Error: Missing columns
                │                               └─► Exit with message
                └─► Yes
                      │
                      └─► Process successfully
                            │
                            ├─► NA values? ──Yes──► Remove & log
                            │
                            └─► Save output
                                  │
                                  ├─► Permission denied? ──Yes──► Error message
                                  │
                                  └─► Success!
```

## 📈 Integration Roadmap

```
Current State (v1.0)
    ├─► Manual conversion ✓
    ├─► Batch processing ✓
    └─► Basic automation ✓

Near Future (v1.1)
    ├─► Snakemake integration
    ├─► Automated GSEA execution
    └─► Quality reports

Future Enhancements (v2.0)
    ├─► GUI wrapper
    ├─► Real-time monitoring
    ├─► Cloud integration
    └─► API endpoints
```

---

**Legend**:
- `┌─┐ └─┘` : Boxes/containers
- `│ ├ └ ─` : Connections
- `▼ ►` : Data flow direction
- `─►` : Process flow
- `🔗` : Integration point
- `✓` : Completed feature
