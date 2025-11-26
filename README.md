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

#### 2. 설정 (`config.yml` 수정)
`config.yml` 파일을 열어 분석할 데이터 경로(`count_data_path`, `metadata_path`), 사용할 DE 분석 방법(`de_analysis.method`), 종(`species`), 출력 폴더(`output_dir` - 고정된 이름 사용, 예: **"results"**) 등 모든 파라미터를 사용자의 환경에 맞게 수정합니다.

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

#### 💡 Snakemake 활용 예시: 여러 조건 실행
만약 여러 데이터셋이나 파라미터 조합으로 분석을 반복하고 싶다면, `config.yml` 파일을 여러 개 만들고(`config_A.yml`, `config_B.yml`) Snakemake 실행 시 --configfile 옵션으로 지정하면 됩니다.

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
│   │   │── 03_enrichment_analysis.R
│   │   └── 04_generate_go_plots.R
│   └── utils/                      # 유틸리티 스크립트
│        └── load_data.R
├── Snakefile                   # ★ (권장) Snakemake 파이프라인 정의
└── README.md                   # 프로젝트 설명 및 사용 방법
```

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
 
---

## � Pairwise 비교 시 정규화 전략

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

## �💡 주요 개념 (FAQ)

### DE 분석 및 데이터 관련
- Normalized Counts와 FPKM의 차이점은 무엇인가요?
    -둘 다 시퀀싱 깊이를 보정하지만, FPKM은 유전자 길이까지 추가로 보정합니다. DE 분석은 동일 유전자를 샘플 간에 비교하는 것이므로, 변하지 않는 �[...]
 
- Ensembl ID는 있는데 왜 유전자 심볼(Gene Symbol)은 비어있나요?
    - 오류가 아니며, 주로 생물학적인 이유 때문입니다. Non-coding Genes (단백질 미생성 유전자), 아직 기능이 밝혀지지 않은 Novel Genes (신규 유전자) 등은 �[...]

### 기능 농축 분석 (Enrichment Analysis) 관련
- pvalue_cutoff와 qvalue_cutoff는 무엇인가요?
    - P-value Cutoff는 "이 GO Term이 농축된 것이 우연일 확률"에 대한 느슨한 1차 필터입니다. (pvalue_cutoff: 0.05 -> 우연일 확률이 5% 미만인 후보들을 일단 선별[...]
    - Q-value Cutoff는 수천 개의 GO Term을 동시에 검정할 때 발생하는 통계적 오류(위양성)를 보정한 엄격한 최종 필터입니다. "유의미하다고 선언한 결과 ��[...]
 
- Dot plot의 Count와 GeneRatio는 무엇을 의미하나요?
    - Count (점의 크기): 분석에 사용한 내 유전자 목록 중, 특정 GO Term에 포함되는 유전자의 개수입니다. 점이 클수록 더 많은 유전자가 그 기능에 관여한[...]
    - GeneRatio (x축 위치): Count를 분석에 사용한 전체 유전자 개수로 나눈 값(비율)입니다. 이 값이 높을수록(오른쪽), 해당 기능이 내 유전자 목록 전체에�[...]
 
- GeneRatio와 Enrichment Score는 비슷한 개념인가요?
    - 아닙니다, 두 개념은 다른 분석 방식에서 나옵니다.
    - GeneRatio: 우리가 사용하는 ORA(Over-Representation Analysis) 방식의 결과로, 미리 선별된 유의미한 유전자 목록 내에서의 비율을 나타냅니다.
    - Enrichment Score: GSEA(Gene Set Enrichment Analysis) 방식의 결과로, 전체 유전자의 발현 순위 안에서 특정 유전자 그룹의 방향성 있는 쏠림 현상을 측정하는 ��[...]
 
- GeneRatio도 Cutoff 기준이 있나요?
    - 아니요, 없습니다. GeneRatio는 P-value처럼 통계적 유의성을 판단하는 기준이 아니라, 영향력의 크기를 나타내는 척도입니다. 먼저 padj < 0.05 기준으로 [...]


