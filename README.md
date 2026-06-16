# RNA-seq 데이터 차발현 유전자 탐색 및 Gene Ontology Analysis 파이프라인 

## 📔 프로젝트 개요

본 프로젝트는 RNA-seq 카운트 데이터를 사용하여 차등 발현 유전자(DEGs)를 식별하고, 유전자 기능(GO) 및 경로(KEGG) 농축 분석을 수행하는 **유연하고 재현성 높은 R 기반 분석 패키지**입니다.

중앙 설정 파일(`config.yml`)을 통해 `DESeq2`, `edgeR`, `limma-voom` 세 가지 주요 DE 분석 알고리즘, 분석 종(Human/Mouse), 실험 디자인, 분석 대상 유전자 그룹, 시각화 옵션, 결과 출력 형식 등 모든 단계를 제어할 수 있습니다.
 
airway 예제 데이터를 사용하여 전체 분석 과정을 즉시 재현할 수 있습니다.

### ✨ 주요 특징
* **다중 DE 알고리즘 지원**: `config.yml`에서 `DESeq2`, `edgeR`, `limma-voom` 선택 가능.

* **중앙화된 설정**: 모든 파라미터와 파일 경로는 `config.yml`에서 관리.
* **유연한 분석 옵션**: 종 선택, 타겟 유전자 그룹(Up/Down/Total), GO Ontology(BP/CC/MF) 조합 분석.
* **자동화된 환경 관리**: `Conda` 환경 파일을 통해 R 및 모든 패키지 의존성 관리.
* **체계적인 결과 관리**: 고정된 출력 폴더(`config.yml`에서 지정)에 모든 결과물과 사용된 설정 파일(`config_used.yml`) 저장.
* **고품질 시각화**: PCA, Volcano plot, GO Dot/Bar plot 등의 세부 속성 제어 가능.
* **다양한 실행 옵션**: 재현성 높은 파이프라인 실행을 위한 **Snakemake**(권장)와 단계별 실행 및 결과 탐색을 위한 **Jupyter Notebook** 지원.
* **Google Drive 자동 업로드**: `rclone`을 통해 파이프라인 완료 후 결과 폴더를 자동으로 Google Drive에 백업. `config.yml`에서 한 줄로 활성화 가능.

---

## 🚀 시작하기: 실행 방법 선택
이 파이프라인은 두 가지 방식으로 실행할 수 있습니다.

1. Snakemake (권장): Linux/macOS 워크스테이션 환경에서 전체 파이프라인을 안정적이고 재현성 있게 실행하는 데 가장 적합합니다. 의존성 관리, 병렬 처리, 부분 재실행 등 강력한 기능을 제공합니다.

2. Jupyter Notebook: Windows 환경 사용자나, 파이프라인의 각 단계를 개별적으로 실행하며 중간 결과를 확인하고 싶을 때 유용합니다.

### 프로젝트 복제 또는 다운로드
방법1과 방법2를 실행하기에 앞서 github repository를 local에 복제하거나 다운로드합니다. 

```bash
#프로젝트 폴더를 다운로드할 폴더로 working directory를 이동 후 아래 명령어를 실행
git clone https://github.com/ibs-CMG-NGS/RNA-Seq_DE_GO_analysis

cd RNA-Seq_DE_GO_analysis #프로젝트 루트로 이동
```

### 🐍 방법 1: Snakemake로 실행하기 (권장)
Snakemake는 Python 기반의 워크플로우 관리 시스템으로, 규칙(rule) 기반으로 파일 간의 의존성을 정의하고 필요한 작업만 자동으로 실행해 줍니다.

#### 1. 환경 설정 (최초 1회)
Snakemake는 실행 환경과 분석 환경을 분리하여 관리하는 것이 가장 좋습니다.

a) R 분석 환경 생성: 프로젝트 루트 폴더에서 다음 명령어를 실행하여 `environment.yml`에 정의된 R 및 관련 패키지 환경(`rna-seq-de-go-analysis`)을 생성합니다.

```bash
conda env create -f environment.yml
```

b) Snakemake 실행 환경 생성: Snakemake 자체를 실행하기 위한 최소한의 환경(`snakemake_env`)을 생성합니다. 프로젝트 루트에 아래 내용으로 `snakemake_environment.yml` 파일을 만듭니다.

```yaml
# snakemake_environment.yml
name: snakemake_env
channels:
  - conda-forge
  - bioconda
dependencies:
  - snakemake-minimal >=7.0 # Snakemake 실행에 필요한 최소 패키지
  - python >=3.8
  - pyyaml # Snakemake가 config.yml을 읽기 위해 필요
  # (선택) DAG 시각화를 위한 graphviz 추가 가능
  # - graphviz
```

그리고 다음 명령어로 환경을 생성합니다.

```bash
conda env create -f snakemake_environment.yml
```

#### 2. 설정 파일 준비 (`configs/` 폴더)

**새로운 설정 파일 구조:**
- `configs/template/config.yml`: Git에 추적되는 템플릿 (수정하지 마세요!)
- `configs/config_*.yml`: 사용자별 설정 파일 (Git에서 자동 무시됨)

a) 템플릿에서 새 설정 파일 생성:

```bash
cp configs/template/config.yml configs/config_my_experiment.yml
```

b) `configs/config_my_experiment.yml` 파일을 열어 다음 항목들을 수정:
- 데이터 경로 (`count_data_path`, `metadata_path`)
- 출력 폴더 (`output_dir`)
- DE 분석 방법 (`de_analysis.method`)
- 비교할 그룹 쌍 (`de_analysis.pairwise_comparisons`)
- 종 (`species`)
- 기타 분석 파라미터

c) `Snakefile`의 첫 부분에서 `CONFIG_FILE` 변수를 수정:

```python
CONFIG_FILE = "configs/config_my_experiment.yml"
```

자세한 설정 방법은 `configs/README.md`를 참조하세요.

#### 3. 파이프라인 실행
a) Snakemake 환경 활성화:

```bash
conda activate snakemake_env
```

b) Snakemake 실행:

```bash
# 전체 파이프라인 실행 (예: 4개 코어 사용)
snakemake --cores 4 --use-conda

# --- 유용한 Snakemake 명령어 ---

# 실행 계획 미리보기 (Dry-run, 실제 실행 안 함)
snakemake --cores 4 --use-conda --dry-run
# 또는 단축 명령어
snakemake -np --use-conda

# 특정 규칙(단계)까지만 실행 (예: DE 분석까지만)
snakemake --cores 4 --use-conda --until run_de_analysis

# 특정 파일 생성 (예: Volcano plot만 다시 생성)
snakemake --cores 4 --use-conda --force results/volcano_plot.png

# 파이프라인 구조(DAG) 이미지로 보기 (graphviz 설치 필요)
# snakemake --dag | dot -Tpng > pipeline_dag.png

# 결과 폴더 및 로그 초기화
snakemake --cores 4 --use-conda clean
```
Snakemake는 `Snakefile`에 정의된 규칙에 따라 `config.yml`을 읽고, 각 R 스크립트를 실행할 때 자동으로 `rna-seq-de-go-analysis` 환경을 활성화하여 분석을 수행합니다. 결과는 `config.yml`의 `output_dir`에 지정된 폴더에 저장됩니다.

#### � 논문용 GO Summary Table 자동 생성

파이프라인을 실행하면 각 pairwise 비교마다 **`final_go_results.xlsx`** 파일이 자동으로 생성됩니다. 이 파일은:
- ✅ 모든 GO enrichment 결과를 하나의 Excel 파일로 통합
- ✅ Gene set별 (UP/DOWN/TOTAL), Ontology별 (BP/CC/MF) 워크시트 자동 구성
- ✅ 상위 유의한 GO term 요약 시트 포함
- ✅ 분석 파라미터 메타데이터 자동 기록
- ✅ 전문적인 서식 (색상 테마, 열 너비, freeze pane) 적용
- ✅ 논문 supplementary material로 바로 제출 가능
- ✅ 개별 연구자에게 전달하기 편리

**출력 위치**: `output/{비교이름}/pairwise/{비교군}_vs_{기준군}/final_go_results.xlsx`

**예시**: `output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_go_results.xlsx`

#### � Config 파일 변경하기

다른 데이터셋을 분석하거나 설정을 변경하려면 **`Snakefile`의 6번째 줄**만 수정하면 됩니다:

```python
# Snakefile 상단 (6번째 줄)
CONFIG_FILE = "config_H2O2_Neuron.yml"  # ← 여기만 변경!
```

예시:
```python
# 다른 프로젝트 분석 시
CONFIG_FILE = "config_Shank2.yml"

# 파라미터 테스트 시
CONFIG_FILE = "config_test_prefilter10.yml"
```

**장점**:
- ✅ 한 곳에서만 수정 → 실수 방지
- ✅ 모든 규칙(rule)에 자동 반영
- ✅ 이전 방식처럼 여러 곳을 수정할 필요 없음

**이전 방식 (번거로움)**:
```python
# 예전에는 각 rule마다 일일이 수정해야 했음
config_file = "config_H2O2_Neuron.yml"  # Rule 1
config_file = "config_H2O2_Neuron.yml"  # Rule 2
config_file = "config_H2O2_Neuron.yml"  # Rule 3
# ... (총 7군데 수정 필요)
```

#### 💡 Snakemake 활용 예시: 여러 조건 실행
여러 config 파일로 순차적으로 분석하려면:

```bash
# 방법 1: CONFIG_FILE 변수만 변경하고 실행 (권장)
# Snakefile 6번째 줄을 수정하고
snakemake --cores 4 --use-conda

# 방법 2: 명령줄에서 config 파일 지정 (고급 사용자)
snakemake --cores 4 --use-conda --configfile config_A.yml
snakemake --cores 4 --use-conda --configfile config_B.yml
```

### 📓 방법 2: Jupyter Notebook으로 실행하기
이 방법은 각 분석 단계를 직접 실행하고 중간 결과를 확인하는 데 유용합니다.

#### 1. 환경 설정 (최초 1회)

a) R 분석 환경 생성: Snakemake 방식과 동일하게, 프로젝트 루트 폴더에서 다음 명령어를 실행하여 `environment.yml`에 정의된 R 분석 환경(`rna-seq-de-go-analysis`)을 생성합니다.

```bash
conda env create -f environment.yml
```

b) Jupyter Notebook 및 R 커널 설치: 생성된 R 분석 환경에 Jupyter Notebook과 IRkernel을 설치합니다.

```bash
conda activate rna-seq-de-go-analysis
conda install jupyter notebook -c conda-forge
Rscript -e 'install.packages("IRkernel"); IRkernel::installspec(user = FALSE)'
```
(user = FALSE는 시스템 전체에 커널을 등록하여 Jupyter에서 찾기 쉽게 합니다.)

#### 2. 설정 (`config.yml` 수정)
`config.yml` 파일을 열어 분석할 데이터 경로, DE 분석 방법, 종, 출력 폴더 등 모든 파라미터를 수정합니다. Jupyter 방식에서는 R 스크립트가 타임스탬프 기반의 출력 폴더를 생성할 수 있으므로 `output_dir` 설정은 선택사항입니다. (현재 R 스크립트는 `config.yml`의 `output_dir`을 사용하도록 되어 있으므로, Snakemake와 동일하게 고정된 폴더 이름을 지정하는 것이 좋습니다.)

#### 3. 파이프라인 실행
a) R 분석 환경 활성화:

```bash
conda activate rna-seq-de-go-analysis
```

b) Jupyter Notebook 실행:

```bash
jupyter notebook
```
c) 노트북 파일 열기 및 실행: 웹 브라우저에서 notebooks/analysis_pipeline.ipynb 파일을 열고, 첫 번째 셀부터 순서대로 실행(Shift + Enter)합니다. 각 셀은 src/analysis/ 폴더의 R 스크립트를 호출하여 해당 분석 단계를 수행합니다.

주의: Jupyter 방식에서는 install_dependencies.R 스크립트를 사용하지 않습니다. 모든 패키지는 Conda 환경 생성 시 설치됩니다. 노트북의 Step 0 셀은 삭제하거나 주석 처리해야 합니다.

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

# Gene ID 타입: "ENSEMBL", "ENTREZID", "SYMBOL" 등
# Raw count 파일의 gene ID 형식을 지정합니다.
# 생략 시 ID 패턴을 통해 자동 감지됩니다.
# - ENSEMBL: ENSG00000... (human) 또는 ENSMUSG00000... (mouse)
# - ENTREZID: 숫자로만 구성된 ID (예: 1234, 5678)
# - SYMBOL: 유전자 심볼 (예: TP53, GAPDH)
gene_id_type: "ENSEMBL"

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
  
  # 고급 옵션
  advanced_options:
    # Pre-filtering: 전체 샘플 합계가 이 값 미만인 유전자는 분석 전 제거
    # 권장값: 10 (0이면 필터링 미적용)
    prefilter_threshold: 10

# enrichment analysis 관련 설정
enrichment:
  # 분석할 유전자 목록: "total", "up", "down" 중에서 원하는 항목만 남기세요.
  gene_lists: ["total", "up", "down"]
  # 실행할 Gene Ontology 목록: "BP", "CC", "MF" 중에서 원하는 항목만 남기세요.
  go_ontologies: ["BP", "CC", "MF"]
  pvalue_cutoff: 0.05
  qvalue_cutoff: 0.2
  
  # GO 유전자 세트 크기 제한 (clusterProfiler enrichGO 파라미터)
  min_gs_size: 10    # 최소 유전자 세트 크기 (너무 작으면 통계적으로 불안정)
  max_gs_size: 500   # 최대 유전자 세트 크기 (너무 크면 너무 일반적이어서 의미 약함)
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

-   **목적:** `config` 설정에 따라 `DESeq2`, `edgeR`, `limma-voom` 중 하나로 차발현 유전자 분석(Defferential expression analysis)을 수행. 유전자 심볼 Annotation, 정규화 카운트와 결과 병합.
-   **입력:** `data/raw/` 폴더의 카운트 및 메타 데이터
-   **출력:** 
    - `final_de_results.csv`: 정규화된 카운트 + 유전자별 DGE 분석 결과 (log2FoldChange, p-value, padj 등)
    - `final_de_results.xlsx`: 위의 결과를 엑셀 파일 형식으로 저장.
    - `config_used.yml`: 분석에 사용된 `config` 파일의 복사본 

## `02_generate_plots.R`

-   **목적:** 분석 방법에 맞는 방식(Vst 또는 LogCPM)으로 PCA 플롯 생성. `config` 설정에 맞춘 Volcano 플롯 생성.
-   **입력:** `final_de_results.csv`, 원본 카운트 데이터
-   **출력:**
    - `pca_plot.png`: 샘플 간의 관계를 보여주는 PCA 플롯
    - `volcano_plot.png`: 유의미한 DEG를 한눈에 보여주는 Volcano 플롯

## `03_enrichment_analysis.R`

-   **목적:** `config` 설정(Up/Down/Total, BP/CC/MF)에 따라 GO/KEGG 농축 분석 수행. 4가지 지표(`FoldEnrichment`/`GeneRatio`, `Count`, `p.adjust`)를 시각화하는 Dot Plot(버블 차트) 생성.
-   **입력:** `final_de_results.csv`
-   **출력:**
    - `go_enrichment_{up|down|total}_{BP|CC|MF}.csv`, `go_dotplot_{up|down|total}_{BP|CC|MF}.png`: Gene Ontology(GO) 분석 결과 및 시각화
    - `kegg_enrichment_{up|down|total}.csv`, `kegg_dotplot_{up|down|total}.png`: KEGG Pathway 분석 결과 및 시각화

## `04_generate_go_plots.R`

-   **목적** `config` 설정에 따라 GO 결과를 취합하여 Ontology별(BP, CC, MF) 3x1 Subplot 형태의 Bar plot 생성. 
-   **입력** `go_enrichment_{up|down|total}_{BP|CC|MF}.csv`
-   **출력** 
    - `go_barplot_{up|down|total}.png`: 3x1의 세로 레이아웃을 갖는 bar plot. 위에서부터 아래로 `BP`, `CC`, `MF` 순으로 정렬 

## `05_generate_go_table.R` ✨ 

-   **목적** GO enrichment 분석 결과를 논문 supplementary material 형식의 통합 Excel 파일로 생성. 개별 연구자 전달 및 논문 제출용으로 최적화된 형식.
-   **입력** 모든 `go_enrichment_{up|down|total}_{BP|CC|MF}.csv` 파일
-   **출력** 
    - `final_go_results.xlsx`: 다중 워크시트로 구성된 Excel 파일
      - **All_Results**: 모든 GO 결과 통합
      - **UP_regulated**, **DOWN_regulated**, **TOTAL_regulated**: Gene set별 결과
      - **Biological_Process**, **Cellular_Component**, **Molecular_Function**: Ontology별 결과
      - **Top_Terms**: 각 카테고리의 상위 20개 가장 유의한 GO term
      - **Analysis_Info**: 분석 파라미터 및 메타데이터
    - 전문적인 서식 적용 (헤더 스타일, 자동 열 너비, freeze pane 등)

---

## 📁 폴더 구조

```
RNA-Seq_DE_GO_analysis/
├── .here                       # 프로젝트 루트를 지정하는 표지 파일
├── config.yml                  # ★ 모든 분석 설정
├── environment.yml             # ★ R 분석 환경 정의
├── snakemake_environment.yml   # ★ (권장) Snakemake 실행 환경 정의
├── data/
│   ├── raw/                    # 원본 데이터 (counts, metadata)
│   └── processed/              # (현재 사용 안함, 필요시 가공 데이터 저장)
├── notebooks/
│   └── analysis_pipeline.ipynb # Jupyter 실행용 노트북
├── output/                     # 분석 결과 저장 폴더
├── src/
│   ├── analysis/                   # 핵심 분석 단계별 스크립트
│   │   ├── 01_run_de_analysis.R
│   │   ├── 02_generate_plots.R
│   │   ├── 03_enrichment_analysis.R
│   │   ├── 04_generate_go_plots.R
│   │   └── 05_generate_go_table.R     # ✨ 논문용 GO 통합 Excel 생성
│   └── utils/                      # 유틸리티 스크립트
│        └── load_data.R
├── Snakefile                   # ★ (권장) Snakemake 파이프라인 정의
├── bridge/                     # 🔗 Pipeline 연결 스크립트
│   ├── convert_de_to_gsea.py      # DE 결과 → GSEA 포맷 변환
│   ├── run_downstream_analysis.sh # 자동화 wrapper 스크립트
│   ├── README.md                  # 상세 사용 설명서
│   └── EXAMPLES.sh                # 사용 예제 모음
└── README.md                   # 프로젝트 설명 및 사용 방법
```

---

## 🔗 Pipeline 연결: Advanced GSEA 분석으로 확장하기

이 파이프라인으로 DE 분석과 기본 GO enrichment를 완료한 후, 더 심화된 GO 분석이나 GSEA(Gene Set Enrichment Analysis)가 필요하다면 **RNA-Seq_GO_GSEA_analysis** 파이프라인과 연결할 수 있습니다.

### Bridge Layer 개요

`bridge/` 디렉토리에는 두 파이프라인을 연결하는 통합 워크플로우가 포함되어 있습니다:

```
RNA-Seq_DE_GO_analysis  →  Bridge Snakefile  →  RNA-Seq_GO_GSEA_analysis
   (R-based, DESeq2)       (자동 변환 워크플로우)    (Python-based, Advanced)
```

### 빠른 시작

#### 🎯 방법 1: Snakemake 통합 워크플로우 (권장) ⭐

```bash
# 모든 비교군 결과를 GSEA 포맷으로 자동 변환
snakemake -s bridge/Snakefile --cores 1

# 특정 비교군만 변환
snakemake -s bridge/Snakefile \
  --config comparison=H2O2_vs_Control experiment=H2O2_Neuron \
  --cores 1

# 병렬 처리 (여러 비교군 동시 변환)
snakemake -s bridge/Snakefile --cores 4
```

**장점**: 
- ✅ 의존성 자동 관리
- ✅ 병렬 처리 지원
- ✅ 변경된 파일만 재처리
- ✅ 로그 자동 생성
- ✅ 이미 익숙한 Snakemake 환경

#### 방법 2: Python 스크립트 직접 실행

```bash
cd bridge

# 단일 비교
python3 convert_de_to_gsea.py \
  --input ../output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \
  --output ../../RNA-Seq_GO_GSEA_analysis/data/H2O2_vs_Control.xlsx

# 일괄 변환
python3 convert_de_to_gsea.py \
  --batch \
  --de-output-dir ../output/H2O2_Neuron \
  --gsea-input-dir ../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline
```

### 주요 기능

- **Snakemake 통합**: 두 파이프라인 모두 Snakemake 기반으로 자연스럽게 연결
- **자동 포맷 변환**: CSV → Excel 변환 및 컬럼명 표준화
- **병렬 처리**: 여러 비교군 동시 처리 가능
- **의존성 추적**: 파일 변경 시 필요한 부분만 재실행
- **품질 검증**: NA 값 제거, 통계적 유의성 기준 정렬
- **메타데이터 포함**: 분석 정보와 요약 통계 자동 생성

### 상세 문서

자세한 사용법, 옵션 설명, 문제 해결 방법은 다음을 참조하세요:

- **⭐ Snakemake 가이드** (권장): [`bridge/SNAKEMAKE_GUIDE.md`](bridge/SNAKEMAKE_GUIDE.md)
- **Python 스크립트 가이드**: [`bridge/README.md`](bridge/README.md)
- **사용 예제 모음**: [`bridge/EXAMPLES.sh`](bridge/EXAMPLES.sh)
- **워크플로우 다이어그램**: [`bridge/WORKFLOW_DIAGRAM.md`](bridge/WORKFLOW_DIAGRAM.md)

---

## ☁️ Google Drive 자동 업로드

파이프라인이 완료되면 `rclone`을 통해 결과 폴더 전체를 Google Drive에 자동으로 업로드할 수 있습니다. 이미 GDrive에 동일한 파일이 있으면 건너뛰므로, 분석을 재실행해도 안전하게 사용할 수 있습니다.

### 사전 준비 (최초 1회)

#### 1. rclone 설치

```bash
# Ubuntu/Debian
sudo apt install rclone

# 또는 공식 설치 스크립트
curl https://rclone.org/install.sh | sudo bash
```

#### 2. Google Drive remote 설정

```bash
rclone config
```

대화형 설정 진행:
1. `n` (New remote)
2. 이름 입력 (예: `my_gdrive`)
3. 스토리지 타입: `drive` (Google Drive)
4. Client ID / Secret: 비워두고 Enter (기본값 사용)
5. Scope: `1` (전체 드라이브 접근)
6. 브라우저 인증 완료

설정된 remote 이름 확인:
```bash
rclone listremotes
# 출력 예: my_gdrive:
```

### config.yml 설정

`configs/config_my_experiment.yml`에 아래 블록을 추가합니다:

```yaml
# -------------------------
# Upload Options
# -------------------------
upload:
  # true로 설정하면 파이프라인 완료 후 Google Drive에 자동 업로드
  gdrive: true

  # rclone remote 이름 (rclone listremotes 로 확인)
  remote: "my_gdrive:"

  # Google Drive 내 업로드할 상위 폴더 이름
  # 실제 업로드 경로: {remote}{dest_folder}/{output_dir 폴더명}/
  # 예: my_gdrive:RNA-Seq_Results/mouse-h2o2-neuron-2026/
  dest_folder: "RNA-Seq_Results"
```

`gdrive: false`로 설정하면 업로드 단계를 건너뜁니다.

### 동작 방식

- **트리거**: `summary_report.html`과 `methods_section.md`가 모두 생성된 후 — 즉 모든 분석이 완료된 시점에 실행됩니다.
- **업로드 경로**: `{remote}{dest_folder}/{output_dir 폴더명}/`
  - 예: `output_dir: output/mouse-h2o2-neuron-2026` → `my_gdrive:RNA-Seq_Results/mouse-h2o2-neuron-2026/`
- **중복 처리**: `rclone copy`를 사용하므로 동일한 파일은 건너뜁니다. 증분 업로드가 가능합니다.
- **완료 표시**: `output/{프로젝트}/.gdrive_upload_done.flag` 파일 생성 (Snakemake 재실행 시 업로드 중복 방지)
- **로그**: `output/{프로젝트}/logs/10_upload_to_gdrive.log`

### 업로드 확인

```bash
# GDrive에 업로드된 파일 목록 확인
rclone ls my_gdrive:RNA-Seq_Results/mouse-h2o2-neuron-2026/

# 업로드 로그 확인
cat output/mouse-h2o2-neuron-2026/logs/10_upload_to_gdrive.log
```

### 문제 해결

| 증상 | 원인 | 해결 방법 |
|:-----|:-----|:----------|
| `rclone: command not found` | rclone 미설치 | `sudo apt install rclone` 실행 |
| `Failed to create file system` | remote 이름 오류 | `rclone listremotes`로 정확한 이름 확인 후 config 수정 |
| `Token expired` | 인증 만료 | `rclone config reconnect my_gdrive:` 실행 |
| 업로드 단계가 실행되지 않음 | `gdrive: false` 설정 | config에서 `gdrive: true`로 변경 |

---

## DE 분석 알고리즘 비교: DESeq2 vs edgeR vs limma-voom

- 이 패키지는 세 가지의 서로 다른 통계적 접근 방식을 사용하는 DE 분석 알고리즘을 지원합니다. 
- 각 방법의 특징을 이해하면 당신의 데이터와 연구 목적에 가장 적합한 도구를 선택할 수 있습니다.

| 구분 | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) | [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) | [limma-voom](https://bioconductor.org/packages/release/bioc/html/limma.html) |
|:------|:--------|:--------|:-------------|
| **핵심 통계 모델** | 음이항 분포 (Negative Binomial) <br> 엄격한 분산 추정 | 음이항 분포 (Negative Binomial) <br> 유연하고 빠른 분산 추정 | 선형 모델 (Linear Model) <br> voom을 통한 분산 정규화 |
| **장점** | - 높은 신뢰도: 위양성(False Positive) 억제에 효과적이라 가장 보수적이고 안정적인 결과 제공 <br> - 편리한 사용법: 함수들이 직관적이고 체계화되어 있음 <br> - 자동화된 정규화: 샘플 간 크기 인자를 자동 계산 | - 빠른 속도: 대용량 데이터셋에서도 신속한 처리 <br> - 유연성: 복잡한 실험 디자인에 대응 가능 <br> - 메모리 효율: 큰 데이터셋에서도 메모리 사용량이 적음 | - 초고속: 세 방법 중 가장 빠른 속도 <br> - 복잡한 디자인 지원: 다중 요인, 반복 측정 등 복잡한 실험 디자인에 강력함 <br> - 유연한 모델링: 선형 모델 프레임워크로 다양한 통계 분석 가능 |
| **단점** | - 느린 속도: 샘플 수가 많아지면 느려짐 <br> - 보수적 경향: 실제 유의미한 유전자를 놓칠 수 있음 (위음성) | - 덜 보수적: DESeq2보다 위양성 가능성이 약간 높을 수 있음 <br> - 샘플 수가 적을 때: 분산 추정이 불안정할 수 있음 | - 전처리 필요: voom 변환이 필수적 <br> - 음이항 분포 미사용: RNA-seq의 과분산을 직접 모델링하지 않음 <br> - 극단적 저발현 유전자: 처리가 까다로울 수 있음 |
| **추천 상황** | 표준 분석, 신뢰도가 중요할 때, 샘플 수가 적을 때 (그룹당 3~5개) | 빠른 분석이 필요하거나 복잡한 실험 디자인, 균형 잡힌 결과를 원할 때 | 매우 큰 데이터셋, 복잡한 실험 디자인(다중 요인, 반복 측정 등), 속도가 중요한 탐색적 분석 |

### 참고문헌

- **DESeq2**: Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550. [https://doi.org/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

- **edgeR**: Robinson, M.D., McCarthy, D.J., Smyth, G.K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1):139-140. [https://doi.org/10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)

- **limma**: Ritchie, M.E., Phipson, B., Wu, D., et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7):e47. [https://doi.org/10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)

---

## 🔬 GO/KEGG 분석 알고리즘: clusterProfiler

이 파이프라인은 GO(Gene Ontology) 및 KEGG(Kyoto Encyclopedia of Genes and Genomes) enrichment 분석에 [**clusterProfiler**](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) 패키지를 사용합니다.

### clusterProfiler란?

**clusterProfiler**는 생물학적 테마를 찾아내기 위한 R/Bioconductor의 표준 도구로, 차등 발현된 유전자 목록이 특정 생물학적 기능이나 경로에 통계적으로 유의하게 농축되어 있는지 분석합니다.

### 핵심 기능

| 분석 유형 | 함수 | 설명 |
|:---------|:-----|:-----|
| **GO Enrichment** | `enrichGO()` | Gene Ontology 용어 농축 분석 <br> - Biological Process (BP): 생물학적 과정 <br> - Cellular Component (CC): 세포 위치 <br> - Molecular Function (MF): 분자 기능 |
| **KEGG Pathway** | `enrichKEGG()` | KEGG 대사 경로 및 신호전달 경로 농축 분석 |
| **시각화** | `dotplot()`, `barplot()` | 농축된 기능/경로의 시각화 |

### 분석 방법론: ORA (Over-Representation Analysis)

이 파이프라인은 **ORA(과대표현 분석)** 방식을 사용합니다:

1. **입력**: 통계적으로 유의한 차등 발현 유전자 목록 (예: padj < 0.05, |log2FC| > 1)
2. **백그라운드**: 분석에 사용된 전체 유전자 세트
3. **통계 검정**: Fisher's exact test 또는 hypergeometric test
   - "이 유전자 목록에 특정 GO term이 우연보다 많이 포함되어 있는가?"
4. **보정**: 다중 검정 보정 (Benjamini-Hochberg FDR)

**예시**:
```
전체 유전자: 20,000개
DEG (Up-regulated): 500개
특정 GO term에 속하는 전체 유전자: 300개
DEG 중 해당 GO term에 속하는 유전자: 50개

→ Fisher's exact test로 이 농축이 우연인지 평가
→ p-value < 0.05, q-value < 0.2이면 유의한 농축으로 판단
```

### 주요 파라미터 (config.yml)

```yaml
enrichment:
  pvalue_cutoff: 0.05      # 1차 필터: 우연일 확률 < 5%
  qvalue_cutoff: 0.2       # 최종 필터: FDR 보정 후 < 20%
  min_gs_size: 10          # 최소 유전자 세트 크기
  max_gs_size: 500         # 최대 유전자 세트 크기
```

- **pvalue_cutoff**: 농축의 통계적 유의성 1차 판단 기준
- **qvalue_cutoff**: 다중 검정 보정 후 최종 유의성 기준 (더 엄격)
- **min_gs_size/max_gs_size**: 너무 작거나 큰 유전자 세트 제외
  - 너무 작으면(< 10): 통계적으로 불안정
  - 너무 크면(> 500): 너무 일반적이어서 생물학적 의미가 약함

### clusterProfiler의 장점

✅ **포괄적인 기능**: GO, KEGG, GSEA 등 다양한 분석 지원  
✅ **최신 데이터베이스**: Annotation 패키지를 통한 자동 업데이트  
✅ **강력한 시각화**: publication-quality 그래프 자동 생성  
✅ **표준화된 워크플로우**: RNA-seq 연구의 사실상 표준  
✅ **다중 종 지원**: Human, Mouse 등 다양한 모델 생물

### ORA vs GSEA 비교

| 특징 | ORA (본 파이프라인) | GSEA |
|:-----|:-------------------|:-----|
| **입력** | 선별된 유의 유전자 목록 | 전체 유전자의 발현 순위 |
| **통계 방법** | Fisher's exact test | Enrichment score + permutation |
| **장점** | 간단하고 직관적, 빠름 | 약한 신호도 감지 가능 |
| **단점** | Cutoff에 민감함 | 계산 복잡, 해석 어려움 |

본 파이프라인은 **ORA** 방식을 사용하여 명확하고 해석하기 쉬운 결과를 제공합니다. GSEA 분석이 필요한 경우 `bridge/` 디렉토리의 변환 도구를 사용하여 결과를 GSEA 형식으로 변환할 수 있습니다.

### 참고문헌

- **clusterProfiler**: Yu, G., Wang, L.G., Han, Y., He, Q.Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS: A Journal of Integrative Biology*, 16(5):284-287. [https://doi.org/10.1089/omi.2011.0118](https://doi.org/10.1089/omi.2011.0118)

- **clusterProfiler 4.0 (업데이트)**: Wu, T., Hu, E., Xu, S., et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*, 2(3):100141. [https://doi.org/10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141)

- **Gene Ontology**: The Gene Ontology Consortium (2021). The Gene Ontology resource: enriching a GOld mine. *Nucleic Acids Research*, 49(D1):D325-D334. [https://doi.org/10.1093/nar/gkaa1113](https://doi.org/10.1093/nar/gkaa1113)

- **KEGG**: Kanehisa, M., Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research*, 28(1):27-30. [https://doi.org/10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)

---

## 🧬 Gene ID 타입 설정

이 파이프라인은 다양한 Gene ID 형식을 지원하며, `config.yml`의 `gene_id_type` 설정으로 입력 데이터의 ID 형식을 지정할 수 있습니다.

### 지원되는 Gene ID 타입

| ID 타입 | 설명 | 예시 | 사용 시기 |
|:--------|:-----|:-----|:---------|
| **ENSEMBL** | Ensembl 유전자 ID | `ENSG00000141510` (human)<br>`ENSMUSG00000051951` (mouse) | RNA-seq 정량화 도구(STAR, Salmon 등)의 기본 출력 |
| **ENTREZID** | NCBI Entrez 유전자 ID | `7157`, `5594` | 오래된 마이크로어레이 데이터<br>NCBI 기반 분석 |
| **SYMBOL** | 유전자 심볼 | `TP53`, `GAPDH` | 읽기 쉬운 결과 제시용<br>문헌 기반 분석 |

### 설정 방법

#### 1. **명시적 지정 (추천)** ⭐

`config.yml`에서 직접 지정:

```yaml
# config.yml
species: "mouse"
gene_id_type: "ENSEMBL"  # 또는 "ENTREZID", "SYMBOL"
```

**장점**: 명확하고 오류 가능성이 낮음

#### 2. **자동 감지 (Fallback)**

`gene_id_type`을 생략하면 첫 번째 유전자 ID의 패턴을 분석하여 자동 감지:

```yaml
# config.yml
species: "mouse"
# gene_id_type: 생략 시 자동 감지
```

**자동 감지 규칙**:
- `ENSMUSG00000...` → `ENSEMBL` (마우스)
- `ENSG00000...` → `ENSEMBL` (사람)
- 숫자로만 구성 (예: `12345`) → `ENTREZID`
- 대문자로 시작하는 영숫자 (예: `TP53`) → `SYMBOL`

**주의**: 자동 감지는 편리하지만, 명시적 지정이 더 안전합니다.

### 내부 처리 방식

GO/KEGG enrichment 분석은 **Entrez ID**를 필요로 합니다. 파이프라인은 자동으로 변환을 수행합니다:

1. **ENTREZID**: 변환 불필요 → 바로 사용
2. **ENSEMBL**: `AnnotationDbi::mapIds()`로 Entrez ID 변환
3. **SYMBOL**: `AnnotationDbi::mapIds()`로 Entrez ID 변환

**변환 예시**:
```r
# ENSEMBL → ENTREZID
ENSMUSG00000051951 → 18999  # Gapdh 유전자

# 일부 유전자는 변환 실패 (NA) 가능 → 자동 제외
```

### 데이터 준비 가이드

#### ✅ Raw count 파일 형식

```csv
Geneid,Sample1,Sample2,Sample3
ENSMUSG00000051951,1234,5678,9012
ENSMUSG00000021803,234,567,890
```

**중요**: 
- 첫 번째 열이 Gene ID
- `gene_id_type`과 일치하는 ID 형식 사용
- 헤더 행 필수

#### ⚠️ 주의사항

1. **일관성**: 모든 유전자가 동일한 ID 타입이어야 함
2. **버전 번호**: ENSEMBL ID의 버전 번호는 제거 권장
   - ❌ `ENSMUSG00000051951.4`
   - ✅ `ENSMUSG00000051951`
3. **변환율**: ENSEMBL → ENTREZID 변환 시 일부 유전자는 매핑 실패 가능 (보통 80-95% 성공)

### 문제 해결

**증상**: "No valid Entrez IDs after conversion" 오류

**해결 방법**:
1. `gene_id_type`이 실제 데이터와 일치하는지 확인
2. Raw count 파일의 첫 몇 줄 확인:
   ```bash
   head -5 data/raw/your_counts.csv
   ```
3. Organism database가 올바른지 확인 (`species` 설정)
4. ENSEMBL ID 버전 번호 제거 시도

---

## 🧹 Pre-filtering: 저발현 유전자 제거

RNA-seq 분석에서 모든 샘플에서 발현되지 않거나 매우 낮게 발현되는 유전자를 사전에 제거하는 것은 **표준적이고 권장되는 절차**입니다. 이 파이프라인은 `config.yml`의 `prefilter_threshold` 설정으로 이를 자동화합니다.

### ❓ 왜 저발현 유전자를 제거해야 할까요?

#### 1. **통계적 검정력 향상** ⭐ (가장 중요)
- RNA-seq 분석에서는 수천~수만 개의 유전자를 동시에 검정합니다
- p-value를 FDR(False Discovery Rate)로 보정할 때, 검정하는 유전자 수가 많을수록 보정이 더 엄격해집니다
- 정보가 없는 유전자를 제거하면 다중 검정 부담이 줄어들어, **실제로 유의한 유전자를 더 잘 찾을 수 있습니다**

**예시**:
```
제거 전: 30,000개 유전자 검정 → FDR 보정 매우 엄격 → 유의한 유전자 500개 발견
제거 후: 15,000개 유전자 검정 → FDR 보정 덜 엄격 → 유의한 유전자 800개 발견
```

#### 2. **분산 추정의 정확도 향상**
- DESeq2와 edgeR은 유전자별 분산을 추정할 때 전체 유전자의 패턴을 사용합니다
- 발현되지 않는 유전자는 분산 추정에 노이즈만 추가하여 정확도를 떨어뜨립니다
- 저발현 유전자 제거 → 더 정확한 분산 추정 → 더 신뢰할 수 있는 통계 결과

#### 3. **계산 효율성**
- 분석 시간 단축 (특히 대규모 데이터셋)
- 메모리 사용량 감소

#### 4. **생물학적 의미**
- 모든 샘플에서 발현이 0이거나 극히 낮은 유전자는:
  - 해당 조직/세포에서 발현되지 않는 유전자
  - 기술적 노이즈일 가능성이 높음
- 이런 유전자는 비교 분석에서 어떤 정보도 제공하지 않습니다

### 🔧 설정 방법

`config.yml` 파일에서 `prefilter_threshold` 값을 설정합니다:

```yaml
de_analysis:
  advanced_options:
    # 전체 샘플의 raw count 합계가 이 값 미만인 유전자 제거
    prefilter_threshold: 10  # 권장값
```

### 📊 권장 필터링 기준

| 설정값 | 의미 | 추천 상황 |
|:-------|:-----|:----------|
| **0** | 필터링 없음 | ❌ 비추천 (다중 검정 부담 증가) |
| **10** | 전체 샘플 합계 < 10 제거 | ✅ **대부분의 경우 권장** (DESeq2 공식 권장) |
| **5** | 전체 샘플 합계 < 5 제거 | 샘플 수가 매우 적을 때 (< 6개) |
| **20** | 전체 샘플 합계 < 20 제거 | 샘플 수가 많고 (> 20개) 보수적 분석 원할 때 |

### 💡 필터링의 실제 효과

**예시 데이터**: 9개 샘플, 57,012개 유전자

```csv
Geneid,Sample1,Sample2,Sample3,...,Sample9
ENSMUSG00000104478,0,0,0,...,0          # 합계 = 0 → 제거됨
ENSMUSG00000104385,0,0,0,...,0          # 합계 = 0 → 제거됨
ENSMUSG00000086053,0,0,2,...,0          # 합계 = 2 → 제거됨 (< 10)
ENSMUSG00000102135,2,12,8,...,0         # 합계 = 58 → 유지됨
ENSMUSG00000051285,3280,3033,3139,...,3940  # 합계 = 29,769 → 유지됨
```

**결과** (`prefilter_threshold: 10` 사용 시):
- 제거: ~30,000개 유전자 (모든 샘플에서 0이거나 극히 낮은 발현)
- 유지: ~27,000개 유전자
- **효과**: FDR 보정이 덜 엄격해져서 더 많은 진짜 DEG를 발견할 수 있음

### ⚠️ 주의사항

1. **너무 엄격한 필터링은 피하세요**
   - `prefilter_threshold`를 너무 높게 설정하면 (예: 100) 실제로 중요한 저발현 유전자를 놓칠 수 있습니다
   - 일부 전사인자, 신호전달 분자 등은 낮게 발현되지만 생물학적으로 중요할 수 있습니다

2. **필터링은 DE 분석 전에 적용됩니다**
   - 최종 결과 파일에는 필터링을 통과한 유전자만 포함됩니다
   - 제거된 유전자는 결과에 나타나지 않습니다

3. **그룹별 발현 패턴은 고려되지 않습니다**
   - 현재 필터링은 "전체 샘플의 합계"만 봅니다
   - 예: Control에서는 0이지만 Treatment에서 높게 발현되는 유전자도 합계가 threshold 이상이면 유지됩니다
   - 이는 의도된 동작으로, DE 분석 전 단계에서는 그룹별 차이를 예단하지 않습니다

### 🔬 더 정교한 필터링 기준 (참고)

본 파이프라인은 단순 합계 기준을 사용하지만, 다른 연구에서는 이런 기준도 사용됩니다:

- **CPM 기반**: 최소 N개 샘플에서 CPM(Counts Per Million) > X
  ```r
  # 예: 최소 3개 샘플에서 CPM > 1
  keep <- rowSums(cpm(counts) > 1) >= 3
  ```
- **Absolute count 기반**: 최소 N개 샘플에서 count > X
  ```r
  # 예: 최소 3개 샘플에서 count > 10
  keep <- rowSums(counts > 10) >= 3
  ```

본 파이프라인의 단순 합계 기준은 사용이 쉽고 대부분의 경우 충분히 효과적입니다.

### 📚 참고 문헌

- DESeq2 공식 vignette: "Pre-filtering the dataset"
  - 권장: 최소 10 reads 이상
  - 링크: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering

- edgeR User's Guide: "Filtering"
  - 권장: CPM 기반 필터링 또는 최소 count 기준

---

## 🔄 Pairwise 비교 시 정규화 전략

다중 그룹 실험에서 pairwise 비교를 수행할 때, 정규화 방법에 따라 결과가 달라질 수 있습니다. 이 파이프라인은 두 가지 정규화 전략을 지원하며, `config.yml`의 `advanced_options.pairwise_normalization` 설정으로 선택할 수 있습니다.

### 1. **Subset Normalization** (`"subset"` - 기본값)

**방법**: 각 pairwise 비교마다 해당 두 그룹의 샘플만으로 정규화를 수행합니다.

**장점**:
- ✅ **독립적인 비교**: 각 비교가 다른 그룹의 영향을 받지 않습니다
- ✅ **보수적인 결과**: 더 엄격한 기준으로 차등 발현 유전자를 선별합니다
- ✅ **명확한 해석**: 두 그룹 간의 순수한 차이만 반영됩니다

**단점**:
- ⚠️ **일관성 부족**: 비교마다 다른 정규화 기준이 적용되어 결과를 직접 비교하기 어렵습니다
- ⚠️ **전체 맥락 손실**: 실험 전체의 생물학적 맥락이 반영되지 않습니다

**추천 상황**:
- 각 비교를 독립적으로 해석하고 싶을 때
- 다른 그룹의 극단적인 발현 패턴이 특정 비교에 영향을 주는 것을 피하고 싶을 때
- 보수적인 DEG 선별이 중요할 때

### 2. **Global Normalization** (`"global"`)

**방법**: 전체 샘플을 사용하여 size factors/normalization factors를 계산한 후, pairwise 비교 시에도 이 값을 유지합니다.

**장점**:
- ✅ **일관된 정규화**: 모든 비교에서 동일한 정규화 기준이 적용됩니다
- ✅ **전체 맥락 유지**: 실험 전체의 생물학적 맥락이 반영됩니다
- ✅ **비교 가능성**: 서로 다른 pairwise 비교 결과를 직접 비교할 수 있습니다
- ✅ **통계적 안정성**: 더 많은 샘플을 사용하여 정규화 계수를 추정하므로 더 안정적입니다

**단점**:
- ⚠️ **다른 그룹의 영향**: 비교에 포함되지 않은 그룹의 발현 패턴도 정규화에 영향을 줍니다
- ⚠️ **민감도 증가**: Subset 방식보다 더 많은 DEG를 검출할 수 있어, 위양성 가능성이 약간 높아질 수 있습니다

**추천 상황**:
- 여러 pairwise 비교 결과를 종합적으로 분석하고 싶을 때
- 시계열 실험이나 용량 반응 실험처럼 전체 실험 맥락이 중요할 때
- 샘플 수가 적어서 정규화 계수 추정의 안정성이 중요할 때

### 설정 방법

`config.yml` 파일의 `advanced_options` 섹션에서 설정합니다:

```yaml
de_analysis:
  advanced_options:
    # Pairwise 비교 시 정규화 방법 선택
    pairwise_normalization: "subset"  # 또는 "global"
```

### 실제 적용 예시

**시나리오**: Control, Low dose, High dose 세 그룹이 있고, "Low vs Control"과 "High vs Control" 두 가지 비교를 수행한다고 가정합니다.

- **Subset 방식**: 
  - "Low vs Control" 비교 시: Low와 Control 샘플만으로 정규화
  - "High vs Control" 비교 시: High와 Control 샘플만으로 정규화
  - → 각 비교의 normalized count 값이 서로 다른 기준으로 계산됨

- **Global 방식**: 
  - Control, Low, High 모든 샘플을 사용하여 정규화
  - "Low vs Control", "High vs Control" 모두 동일한 normalized count 값 사용
  - → 두 비교 결과를 직접 비교 가능 (예: "Low에서는 유의하지 않지만 High에서는 유의한 유전자" 식별 가능)

---

## 💡 주요 개념 (FAQ)

### Pre-filtering 및 데이터 품질 관련
- **모든 샘플에서 raw count가 0인 유전자는 삭제해도 되나요?**
    - **네, 반드시 삭제하는 것이 좋습니다.** 오히려 삭제하지 않으면 다중 검정 보정이 과도하게 엄격해져서 실제로 유의한 유전자를 놓칠 수 있습니다. `prefilter_threshold: 10` 설정을 권장합니다. 자세한 내용은 위의 "🧹 Pre-filtering: 저발현 유전자 제거" 섹션을 참고하세요.

- **prefilter_threshold를 0으로 설정하면 어떻게 되나요?**
    - 필터링이 적용되지 않아 모든 유전자가 분석에 포함됩니다. 이는 권장하지 않으며, 분석 시간이 길어지고 다중 검정 부담으로 인해 통계적 검정력이 떨어집니다.

- **prefilter_threshold를 얼마로 설정해야 하나요?**
    - **대부분의 경우 10을 권장합니다** (DESeq2 공식 권장 기준). 샘플 수가 매우 적다면(<6개) 5로 낮추고, 샘플 수가 매우 많다면(>20개) 20으로 높일 수 있습니다.

### DE 분석 및 데이터 관련
- **Normalized Counts와 FPKM의 차이점은 무엇인가요?**
    - 둘 다 시퀀싱 깊이를 보정하지만, FPKM은 유전자 길이까지 추가로 보정합니다. DE 분석은 동일 유전자를 샘플 간에 비교하는 것이므로, 변하지 않는 유전자 길이를 굳이 보정할 필요가 없습니다. 따라서 DESeq2, edgeR 등은 유전자 길이 보정 없이 라이브러리 크기(시퀀싱 깊이)만 보정한 normalized counts를 사용합니다.
 
- **Ensembl ID는 있는데 왜 유전자 심볼(Gene Symbol)은 비어있나요?**
    - 오류가 아니며, 주로 생물학적인 이유 때문입니다. Non-coding Genes (단백질 미생성 유전자), 아직 기능이 밝혀지지 않은 Novel Genes (신규 유전자) 등은 공식적인 Gene Symbol이 없을 수 있습니다. 이런 경우 Ensembl ID만 존재하며, 추후 연구가 진행되면 심볼이 부여됩니다.

### 기능 농축 분석 (Enrichment Analysis) 관련
- **pvalue_cutoff와 qvalue_cutoff는 무엇인가요?**
    - **P-value Cutoff**: "이 GO Term이 농축된 것이 우연일 확률"에 대한 느슨한 1차 필터입니다. (pvalue_cutoff: 0.05 → 우연일 확률이 5% 미만인 후보들을 일단 선별)
    - **Q-value Cutoff**: 수천 개의 GO Term을 동시에 검정할 때 발생하는 통계적 오류(위양성)를 보정한 엄격한 최종 필터입니다. "유의미하다고 선언한 결과 중 실제로 틀릴 확률"을 제어합니다. (qvalue_cutoff: 0.2 → 선언된 결과 중 최대 20%가 위양성일 수 있음을 허용)
 
- **Dot plot의 Count와 GeneRatio는 무엇을 의미하나요?**
    - **Count (점의 크기)**: 분석에 사용한 내 유전자 목록 중, 특정 GO Term에 포함되는 유전자의 개수입니다. 점이 클수록 더 많은 유전자가 그 기능에 관여합니다.
    - **GeneRatio (x축 위치)**: Count를 분석에 사용한 전체 유전자 개수로 나눈 값(비율)입니다. 이 값이 높을수록(오른쪽), 해당 기능이 내 유전자 목록 전체에서 차지하는 비중이 큽니다.
 
- **GeneRatio와 Enrichment Score는 비슷한 개념인가요?**
    - 아닙니다, 두 개념은 다른 분석 방식에서 나옵니다.
    - **GeneRatio**: 우리가 사용하는 ORA(Over-Representation Analysis) 방식의 결과로, 미리 선별된 유의미한 유전자 목록 내에서의 비율을 나타냅니다.
    - **Enrichment Score**: GSEA(Gene Set Enrichment Analysis) 방식의 결과로, 전체 유전자의 발현 순위 안에서 특정 유전자 그룹의 방향성 있는 쏠림 현상을 측정하는 지표입니다. 이 파이프라인에서는 ORA 방식을 사용하므로 Enrichment Score는 계산되지 않습니다.
 
- **GeneRatio도 Cutoff 기준이 있나요?**
    - 아니요, 없습니다. GeneRatio는 P-value처럼 통계적 유의성을 판단하는 기준이 아니라, 영향력의 크기를 나타내는 척도입니다. 먼저 padj < 0.05 기준으로 유의미한 GO Term을 선별한 후, GeneRatio를 참고하여 "얼마나 많은 유전자가 관여하는지" 해석에 활용합니다.


