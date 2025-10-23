# RNA-seq 데이터 차발현 유전자 탐색 및 GO/KEGG분석 파이프라인 

## 📔 프로젝트 개요

본 프로젝트는 RNA-seq 카운트 데이터를 사용하여 차등 발현 유전자(DEG) 분석, 유전자 기능(GO) 및 경로(KEGG) 농축 분석을 수행하는 유연하고 재현성 높은 R 기반 파이프라인입니다.

이 파이프라인의 가장 큰 특징은 `config.yml` 설정 파일 하나로 `DESeq2`, `edgeR`, `limma-voom` 세 가지의 주요 DE 분석 알고리즘을 자유롭게 선택할 수 있다는 점입니다. 또한 분석 종(Human/Mouse), 실험 디자인, 분석 대상 유전자 그룹, 시각화 옵션, 결과 출력 형식 등 모든 단계를 중앙에서 제어할 수 있도록 설계되었습니다.

airway 예제 데이터를 사용하여 전체 분석 과정을 즉시 재현할 수 있습니다.

### ✨ 주요 특징
- 다중 분석 알고리즘 지원: `DESeq2`, `edgeR`, `limma-voom` 중 원하는 DE 분석 방법을 `config` 파일에서 선택하여 실행 가능합니다.

- 중앙화된 설정: `config.yml` 파일을 통해 모든 파라미터를 중앙에서 관리하여 코드 수정 없이 분석을 제어합니다.

- 유연한 분석 옵션:
    * 종 선택: `Human` 또는 `Mouse` 데이터 분석을 쉽게 전환할 수 있습니다.
    * 타겟 유전자 선택: 전체 DEG, 상향(Up-regulated), 하향(Down-regulated) 조절 유전자 그룹에 대한 기능 분석을 선택적으로 수행합니다.
    * GO Ontology 선택: `BP`, `CC`, `MF` 중 원하는 GO 분석을 자유롭게 조합할 수 있습니다.

- 자동화된 환경 구축: 분석에 필요한 R 패키지를 자동으로 확인하고 설치하여 높은 재현성을 보장합니다.

- 체계적인 결과 관리:
    * 실행 시마다 타임스탬프가 찍힌 폴더에 모든 결과물이 저장됩니다.
    * 분석 옵션에 따라 결과 파일 이름이 동적으로 생성됩니다 (예: go_enrichment_up_BP.csv).
    * CSV와 함께 사용이 편리한 Excel(.xlsx) 파일로 결과 export를 지원합니다.

- 고품질 시각화: config.yml을 통해 플롯의 색상, 폰트 크기, 라벨링 등 세부 속성을 조정하여 논문에 바로 사용할 수 있는 수준의 그래프를 생성합니다.

---

## 🚀 시작하기

### 사전 요구사항

-   [R](https://www.r-project.org/) (버전 4.2 이상 권장)
-   [RStudio](https://posit.co/download/rstudio-desktop/) 또는 [Visual Studio Code](https://code.visualstudio.com/) (+ R 확장 프로그램)
-   Jupyter Notebook 환경에서 R을 사용하기 위한 [IRkernel](https://irkernel.github.io/)

### 설치 및 실행

1.  **프로젝트 복제 또는 다운로드**
    ```bash
    git clone https://github.com/ibs-CMG-NGS/RNA-Seq_DE_GO_analysis
    cd RNA-Seq_DE_GO_analysis
    ```

2.  **`config.yml` 설정**
    프로젝트의 `config.yml` 파일을 열어 자신의 데이터와 분석 목적에 맞게 파라미터를 수정합니다. (아래 config.yml 설명 참고) 특히 `de_analysis` 섹션의 `method`를 확인하세요. 

3.  **파이프라인 실행**

    **방법 1: Jupyter Notebook (GUI 환경)**
    
    `notebooks/analysis_pipeline.ipynb` 파일을 열고, 첫 번째 셀부터 순서대로 실행(`Shift + Enter`)하세요.
    -   **최초 실행 시:** `Step 0`에서 분석에 필요한 모든 R 패키지를 자동으로 설치합니다.
    -   **분석 완료 후:** `output/` 폴더 아래에 `YYYY-MM-DD_HH-MM-SS` 형식의 폴더가 생성되고, 그 안에 모든 결과물(CSV, PNG)이 저장됩니다.

    **방법 2: 커맨드 라인 (CLI) - Snakemake 등 파이프라인 도구 호환**
    
    리눅스 서버나 HPC 환경에서 Snakemake 등의 워크플로우 관리 도구와 함께 사용할 수 있습니다.
    
    ```bash
    # 전체 파이프라인 실행 (의존성 설치 포함)
    Rscript run_pipeline.R --install-deps
    
    # 기본 실행 (config.yml 사용)
    Rscript run_pipeline.R
    
    # 또는 shell 스크립트 사용
    ./run_pipeline.sh
    
    # 특정 config 파일 지정
    Rscript run_pipeline.R --config my_config.yml
    
    # 출력 디렉토리 지정
    Rscript run_pipeline.R --output results/my_analysis
    
    # 특정 분석 단계만 실행
    Rscript run_pipeline.R --step de           # DE 분석만
    Rscript run_pipeline.R --step plots        # 플롯 생성만
    Rscript run_pipeline.R --step enrichment   # Enrichment 분석만
    Rscript run_pipeline.R --step go_plots     # GO 플롯만
    
    # Verbose 모드로 실행 (디버깅용)
    Rscript run_pipeline.R --verbose
    
    # 도움말 보기
    Rscript run_pipeline.R --help
    ```
    
    **CLI 옵션:**
    - `--config, -c`: 설정 파일 경로 (기본값: config.yml)
    - `--output, -o`: 출력 디렉토리 경로 (기본값: output/YYYY-MM-DD_HH-MM-SS)
    - `--step, -s`: 실행할 분석 단계 (all, de, plots, enrichment, go_plots)
    - `--install-deps`: R 패키지 자동 설치
    - `--verbose, -v`: 상세 로그 출력
    - `--help`: 도움말 표시

### Snakemake 워크플로우 통합

프로젝트에 포함된 `Snakefile`을 사용하여 Snakemake로 파이프라인을 실행할 수 있습니다:

```bash
# 전체 파이프라인 실행 (4 코어 사용)
snakemake --cores 4

# Dry-run (실제 실행 없이 계획만 확인)
snakemake --cores 4 --dry-run

# 특정 단계까지만 실행
snakemake --cores 4 de_analysis        # DE 분석만
snakemake --cores 4 enrichment         # Enrichment까지

# DAG 시각화
snakemake --dag | dot -Tpng > workflow.png

# 결과 정리
snakemake clean
```

Snakemake는 자동으로 의존성을 관리하고, 병렬 처리를 수행하며, 실패 시 재시작 기능을 제공합니다.

### Nextflow 워크플로우 통합

Nextflow를 선호하는 경우, 프로젝트에 포함된 `nextflow_pipeline.nf`를 사용할 수 있습니다:

```bash
# 전체 파이프라인 실행
nextflow run nextflow_pipeline.nf

# 커스텀 config 사용
nextflow run nextflow_pipeline.nf --config my_config.yml

# 출력 디렉토리 지정
nextflow run nextflow_pipeline.nf --outdir results/exp1

# 의존성 설치 포함
nextflow run nextflow_pipeline.nf --install_deps true
```

---

## ⚙️ `config.yml` 상세 설명

이 파일은 분석 파이프라인의 모든 것을 제어하는 두뇌 역할을 합니다.

```yaml
# -------------------------
# Project Settings
# -------------------------
# 분석할 종(species)을 선택합니다: "human" 또는 "mouse"
# 이 설정에 따라 아래 organism_db와 kegg_code가 자동으로 결정됩니다.
species: "human"

# -------------------------
# File Paths
# -------------------------
count_data_path: "data/raw/airway_scaledcounts.csv"
metadata_path: "data/raw/airway_metadata.csv"
output_dir: "output"

# -------------------------
# Analysis Parameters
# -------------------------
# DE 분석 관련 설정
de_analysis:
  # 사용할 DGE 분석 방법을 선택: "DESeq2", "edgeR", "limma-voom"
  method: "DESeq2"
  # 메타데이터의 컬럼을 이용한 실험 디자인
  design_formula: "~ dex"
  # 유의성 기준 (Adjusted p-value)
  padj_cutoff: 0.05
  # Log2 Fold Change 절대값 기준
  log2fc_cutoff: 1.0

# enrichment analysis 관련 설정
enrichment:
  # 분석할 유전자 목록: "total", "up", "down" 중에서 원하는 항목만 남기세요.
  gene_lists: ["total", "up", "down"]
  # 실행할 Gene Ontology 목록: "BP", "CC", "MF" 중에서 원하는 항목만 남기세요.
  go_ontologies: ["BP", "CC", "MF"]
  pvalue_cutoff: 0.05
  qvalue_cutoff: 0.2
# -------------------------
# Annotation Databases (자동 설정)
# -------------------------
# species 설정에 따라 결정되는 값들이므로 직접 수정할 필요가 없습니다.
databases:
  human:
    organism_db: "org.Hs.eg.db"
    kegg_code: "hsa"
  mouse:
    organism_db: "org.Mm.eg.db"
    kegg_code: "mmu"

# -------------------------
# Plot Aesthetics (시각화 속성)
# -------------------------
plot_aesthetics:
  # Volcano Plot 속성
  volcano:
    title: "Volcano Plot"
    up_color: "#FF5733"      # 상향 조절 유전자 색상
    down_color: "#3375FF"    # 하향 조절 유전자 색상
    base_color: "grey"       # 그 외 유전자 색상
    point_size: 2.5          # 점 크기
    label_top_n: 0          # 상위 N개 유전자에 라벨 표시 (0이면 표시 안함)
    base_font_size: 14       # 기본 폰트 크기

  # Dot Plot (GO & KEGG) 속성
  dotplot:
    show_n_categories: 15    # 표시할 상위 카테고리 개수
    font_size: 12            # 폰트 크기
    low_color: "#FFC300"     # p-value가 낮을 때의 색상 (진한 색)
    high_color: "#C70039"    # p-value가 높을 때의 색상 (연한 색)

# -------------------------
# GO Bar Plot Settings
# -------------------------
go_barplot:
  # 플롯에 표시할 상위 GO Term 개수
  top_n: 10
  # 종합 플롯에 포함할 Ontology 목록: "BP", "CC", "MF" 중 선택
  namespaces: ["BP", "CC", "MF"]
  # 각 Ontology에 대한 막대 색상
  colors:
    BP: "#E57373"  # Biological Process
    CC: "#64B5F6"  # Cellular Component
    MF: "#81C784"  # Molecular FunctionS

# -------------------------
# Export Options
# -------------------------
export:
  # true로 설정하면 .xlsx 엑셀 파일도 함께 생성합니다.
  export_to_excel: true
```

---

## 🔬 스크립트 및 결과물 설명

### `01_run_de_analysis.R`

-   **목적:** 차발현 유전자(Differential Gene Expression, DGE) 분석을 수행하고, 유전자 심볼을 Annotation합니다.
-   **입력:** `data/raw/` 폴더의 카운트 및 메타 데이터
-   **출력:** 
    - `final_de_results.csv`: 정규화된 카운트 + 유전자별 DGE 분석 결과 (log2FoldChange, p-value, padj 등)
    - `final_de_results.xlsx`: 위의 결과를 엑셀 파일

## `02_generate_plots.R`

-   **목적:** DGE 분석 결과를 바탕으로 주요 시각화 자료를 생성합니다.
-   **입력:** `final_de_results.csv`, 원본 카운트 데이터
-   **출력:**
    - `pca_plot.png`: 샘플 간의 관계를 보여주는 PCA 플롯
    - `volcano_plot.png`: 유의미한 DEG를 한눈에 보여주는 Volcano 플롯

## `03_enrichment_analysis.R`

-   **목적:** config 파일 설정에 따라 지정된 유전자 그룹과 Ontology 조합에 대한 모든 기능 농축 분석(GO & KEGG)을 자동으로 수행합니다.
-   **입력:** `final_de_results.csv`
-   **출력:**
    - `go_enrichment_{up|down|total}_{BP|CC|MF}.csv`, `go_dotplot_{up|down|total}_{BP|CC|MF}.png`: Gene Ontology(GO) 분석 결과 및 시각화, 차발현 유전자 속성과 
    - `kegg_enrichment_{up|down|total}.csv`, `kegg_dotplot_{up|down|total}.png`: KEGG Pathway 분석 결과 및 시각화

## `04_generate_go_plots.R`

-   **목적** config 파일 설정과 분석 결과를 바탕으로 GO 분석 결과를 bar plot으로 생성합니다. 
-   **입력** `go_enrichment_{up|down|total}_{BP|CC|MF}.csv`
-   **출력** 
    - `go_barplot_{up|down|total}.png`: 3x1의 세로 레이아웃을 갖는 bar plot. 위에서부터 아래로 `BP`, `CC`, `MF` 순으로 정렬 

---

## 📁 폴더 구조

```
RNA-Seq_DE_GO_analysis/
├── .here                       # 프로젝트 루트를 지정하는 표지 파일
├── config.yml                  # 모든 경로 및 파라미터 설정
├── data/
│   ├── raw/                    # 원본 데이터 (counts, metadata)
│   └── processed/              # (현재 사용 안함, 필요시 가공 데이터 저장)
├── notebooks/
│   └── analysis_pipeline.ipynb # 전체 파이프라인을 실행하는 메인 노트북
├── output/                     # 모든 분석 결과물이 저장되는 폴더
└── src/
├── analysis/                   # 핵심 분석 단계별 스크립트
│   ├── 01_run_deseq2.R
│   ├── 02_generate_plots.R
│   │── 03_enrichment_analysis.R
│   └── 04_generate_go_plots.R
└── utils/                      # 유틸리티 스크립트
    ├── install_dependencies.R
    └── load_data.R
```

---

## DE 분석 알고리즘 비교: DESeq2 vs edgeR vs limma-voom

- 이 패키지는 세 가지의 서로 다른 통계적 접근 방식을 사용하는 DE 분석 알고리즘을 지원합니다. 
- 각 방법의 특징을 이해하면 당신의 데이터와 연구 목적에 가장 적합한 도구를 선택할 수 있습니다.

| 구분 | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) | [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) | [limma-voom](https://bioconductor.org/packages/release/bioc/html/limma.html) |
|:------|:--------|:--------|:-------------|
| **핵심 통계 모델** | 음이항 분포 (Negative Binomial) <br> 엄격한 분산 추정 | 음이항 분포 (Negative Binomial) <br> 유연하고 빠른 분산 추정 | 선형 모델 (Linear Model) + 정규분포 <br> 로그 변환된 카운트 데이터에 적용 |
| **장점** | - 높은 신뢰도: 위양성(False Positive) 억제에 효과적이라 가장 보수적이고 안정적인 결과 제공 <br> - 편리한 사용법: 함수들이 직관적이고 체계적 <br> - 가장 큰 커뮤니티: 관련 자료를 찾기 쉬움 | - 가장 빠른 속도: 분석 속도가 매우 빠름 <br> - 높은 유연성: 복잡한 실험 디자인(GLM)에 강함 <br> - 균형 잡힌 성능: 속도와 정확도의 균형이 우수 | - 많은 샘플 수에 강력: 그룹당 샘플 수가 많아질수록(6개 이상) 성능 향상 <br> - 검증된 통계 기법: 오랫동안 검증된 선형 모델과 베이지안 통계 사용 |
| **단점** | - 느린 속도: 샘플 수가 많아지면 느려짐 <br> - 보수적 경향: 실제 유의미한 유전자를 놓칠 수 있음 (위음성) | - 덜 보수적: DESeq2보다 위양성이 약간 더 높을 수 있음 | - 적은 샘플 수에 불리: 그룹당 샘플 수가 적을 때(3개 미만) 통계적 힘 약화 |
| **추천 상황** | 표준 분석, 신뢰도가 중요할 때, 샘플 수가 적을 때 (그룹당 3~5개) | 빠른 분석이 필요하거나 복잡한 실험 디자인, 균형 잡힌 결과를 원할 때 | 샘플 수가 충분히 많을 때 (그룹당 6개 이상) |

---

## 💡 주요 개념 (FAQ)

### DE 분석 및 데이터 관련
- Normalized Counts와 FPKM의 차이점은 무엇인가요?
    -둘 다 시퀀싱 깊이를 보정하지만, FPKM은 유전자 길이까지 추가로 보정합니다. DE 분석은 동일 유전자를 샘플 간에 비교하는 것이므로, 변하지 않는 값인 유전자 길이를 보정할 필요가 없습니다. 따라서 DESeq2, edgeR 등이 사용하는 `normalized_counts`(또는 CPM) 방식이 현대 DE 분석의 표준입니다.

- Ensembl ID는 있는데 왜 유전자 심볼(Gene Symbol)은 비어있나요?
    - 오류가 아니며, 주로 생물학적인 이유 때문입니다. Non-coding Genes (단백질 미생성 유전자), 아직 기능이 밝혀지지 않은 Novel Genes (신규 유전자) 등은 공식 유전자 심볼이 없는 경우가 많습니다.

### 기능 농축 분석 (Enrichment Analysis) 관련
- pvalue_cutoff와 qvalue_cutoff는 무엇인가요?
    - P-value Cutoff는 "이 GO Term이 농축된 것이 우연일 확률"에 대한 느슨한 1차 필터입니다. (pvalue_cutoff: 0.05 -> 우연일 확률이 5% 미만인 후보들을 일단 선별)
    - Q-value Cutoff는 수천 개의 GO Term을 동시에 검정할 때 발생하는 통계적 오류(위양성)를 보정한 엄격한 최종 필터입니다. "유의미하다고 선언한 결과 중 최대 20%는 위양성일 수 있음을 감수하겠다"는 의미로, **FDR(False Discovery Rate)**과 같은 개념입니다.

- Dot plot의 Count와 GeneRatio는 무엇을 의미하나요?
    - Count (점의 크기): 분석에 사용한 내 유전자 목록 중, 특정 GO Term에 포함되는 유전자의 개수입니다. 점이 클수록 더 많은 유전자가 그 기능에 관여한다는 의미입니다.
    - GeneRatio (x축 위치): Count를 분석에 사용한 전체 유전자 개수로 나눈 값(비율)입니다. 이 값이 높을수록(오른쪽), 해당 기능이 내 유전자 목록 전체에서 차지하는 생물학적 중요도 또는 비중이 크다는 의미입니다.

- GeneRatio와 Enrichment Score는 비슷한 개념인가요?
    - 아닙니다, 두 개념은 다른 분석 방식에서 나옵니다.
    - GeneRatio: 우리가 사용하는 ORA(Over-Representation Analysis) 방식의 결과로, 미리 선별된 유의미한 유전자 목록 내에서의 비율을 나타냅니다.
    - Enrichment Score: GSEA(Gene Set Enrichment Analysis) 방식의 결과로, 전체 유전자의 발현 순위 안에서 특정 유전자 그룹의 방향성 있는 쏠림 현상을 측정하는 더 복잡한 통계 값입니다.

- GeneRatio도 Cutoff 기준이 있나요?


    - 아니요, 없습니다. GeneRatio는 P-value처럼 통계적 유의성을 판단하는 기준이 아니라, 영향력의 크기를 나타내는 척도입니다. 먼저 padj < 0.05 기준으로 통계적으로 유의미한 Term들을 걸러낸 후, 그중에서 GeneRatio가 높은 순서대로 결과를 해석하며 중요도의 우선순위를 매기는 데 사용합니다.


