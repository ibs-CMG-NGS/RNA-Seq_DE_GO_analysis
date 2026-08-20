# RNA-seq 데이터 차등발현 유전자 탐색 및 Gene Ontology Analysis 파이프라인

## 📔 프로젝트 개요

본 프로젝트는 RNA-seq 카운트 데이터를 사용하여 차등 발현 유전자(DEGs)를 식별하고, 유전자 기능(GO) 및 경로(KEGG) 농축 분석을 수행하는 **유연하고 재현성 높은 R/Snakemake 기반 분석 파이프라인**입니다.

중앙 설정 파일(`configs/config_*.yml`)을 통해 `DESeq2`/`edgeR`/`limma-voom` DE 분석, 다중 그룹 pairwise 비교, omnibus test, time-series 분석, coexpression module 분석, GO term clustering, QC 리포트, CMG-SeqViewer 연동, Google Drive 자동 업로드까지 전체 워크플로우를 제어합니다.

airway 예제 데이터를 사용하여 전체 분석 과정을 즉시 재현할 수 있습니다.

### ✨ 주요 특징
* **다중 DE 알고리즘 지원**: `DESeq2`, `edgeR`, `limma-voom` 중 선택.
* **다중 그룹 pairwise 비교**: 그룹이 3개 이상이어도 `pairwise_comparisons` 목록만으로 원하는 비교쌍을 모두 자동 실행.
* **Omnibus test + Multi-group export**: 전체 그룹 간 차이를 보는 LRT/F-test와, 그 결과를 CMG-SeqViewer MULTI_GROUP 뷰용으로 내보내는 기능.
* **Time-series 분석 (maSigPro)**: 적은 replicate/time point에서도 안정적으로 동작하는 count 기반 회귀로 발현 궤적을 검정하고, 유의 유전자를 패턴별로 클러스터링. `include_groups`로 한 프로젝트 안에 시간축이 다른 여러 시계열(예: acute/chronic)도 구성 가능.
* **Coexpression Module 분석**: Omnibus test로 걸러진 유의 유전자 서브셋에 대해 `DEGreport::degPatterns`로 발현 패턴 기반 클러스터링. `include_groups`로 특정 그룹 조합만 대상으로 제한 가능.
* **GO term 심화 분석 3종**:
  - **Term clustering**: Jaccard 유전자 중복도 기반 `pairwise_termsim()` + `treeplot()`으로 유의 GO term을 클러스터링하고, CMG-SeqViewer가 바로 임포트할 수 있는 `final_go_clustered_results.xlsx`로 내보냄.
  - **GO Slim rollup**: `gofilter()`로 GO DAG 상위 level만 남겨 up/down 대칭 bar chart(`go_slim_overview_*.png`) 생성 — 실험 전체 방향을 빠르게 파악.
  - **rrvgo 의미론적 축약**: GO DAG 의미 거리 기반으로 유의 term에 `parentTerm` 라벨을 부여하고 treemap/scatter plot 생성.
* **QC 리포트 자동 생성**: 샘플 전체(Global) 및 비교쌍별(Pairwise) QC 플롯 + HTML 리포트.
* **CMG-SeqViewer 연동**: DE/GO/multi-group/time-series/coexpression 결과를 parquet + staging JSON으로 자동 export, 프로젝트 전체를 `metadata.json`으로 집계.
* **Methods section 자동 생성**: 사용된 통계 방법/파라미터를 논문 Methods 초안 형태로 자동 기록.
* **논문용 GO Summary Table 자동 생성**: `final_go_results.xlsx` — Gene Set(UP/DOWN/TOTAL) × Ontology(BP/CC/MF)별 워크시트로 자동 구성.
* **Google Drive 자동 업로드**: `rclone`을 통해 파이프라인 완료 후 결과 폴더를 자동으로 백업 (증분 업로드, 중복 방지).
* **배치 실행 스크립트**: `run_batch.sh`로 여러 프로젝트(config)를 한 번에 순차 실행 가능.

---

## 🚀 시작하기

### 1. 저장소 클론

```bash
# 프로젝트를 다운로드할 폴더로 이동 후
git clone https://github.com/ibs-CMG-NGS/RNA-Seq_DE_GO_analysis
cd RNA-Seq_DE_GO_analysis
```

#### 특정 브랜치를 사용하려면

```bash
# 원격 저장소의 브랜치 목록 확인 (clone 전에도 조회 가능)
git ls-remote --heads https://github.com/ibs-CMG-NGS/RNA-Seq_DE_GO_analysis

# 방법 A: clone할 때 바로 브랜치 지정
git clone -b <branch-name> https://github.com/ibs-CMG-NGS/RNA-Seq_DE_GO_analysis
cd RNA-Seq_DE_GO_analysis

# 방법 B: 기본 브랜치(main)로 clone한 뒤 전환
git clone https://github.com/ibs-CMG-NGS/RNA-Seq_DE_GO_analysis
cd RNA-Seq_DE_GO_analysis
git branch -a              # 로컬 + 원격 브랜치 전체 목록
git checkout <branch-name> # 예: git checkout add-multiple-comparisons
```

### 2. 환경 설정 (최초 1회)

이 파이프라인은 **Snakemake, rclone(Google Drive 업로드), R 및 모든 분석 패키지를 하나의 conda 환경(`rna-seq-de-go-analysis`)에서 함께 관리**합니다. Snakemake용 별도 환경을 만들 필요가 없습니다.

```bash
conda env create -f environment.yml
conda activate rna-seq-de-go-analysis
```

이 환경 안에 `snakemake`, `rclone`, `R` 및 `DESeq2`/`edgeR`/`limma`/`clusterProfiler`/`DEGreport`/`maSigPro`/`rrvgo` 등 분석에 필요한 모든 R 패키지가 함께 설치됩니다.

### 3. 설정 파일 준비 (`configs/` 폴더)

- `configs/template/config.yml`: Git에 추적되는 템플릿 (직접 수정하지 말고 복사해서 사용)
- `configs/config_*.yml`: 프로젝트별 설정 파일 (Git에서 자동 무시됨, `.gitignore` 참고)

```bash
cp configs/template/config.yml configs/config_my_experiment.yml
```

`configs/config_my_experiment.yml`을 열어 최소한 아래 항목을 수정합니다:
- `species`, `gene_id_type`
- `count_data_path`, `metadata_path`, `output_dir`
- `de_analysis.design_formula`, `de_analysis.group_variable`
- `de_analysis.pairwise_comparisons`

나머지 옵션(omnibus test, time-series, coexpression modules, GO term clustering, QC, GDrive 업로드 등)은 아래 [⚙️ config.yml 상세 설명](#️-configyml-상세-설명)을 참고하세요.

### 4. 파이프라인 실행

Snakemake는 `--configfile`로 어떤 config를 쓸지 그때그때 지정합니다. **Snakefile을 직접 수정할 필요가 없습니다.**

```bash
# 단일 프로젝트 실행
snakemake --configfile configs/config_my_experiment.yml --cores 4 --use-conda

# 실행 계획만 미리 보기 (Dry-run)
snakemake --configfile configs/config_my_experiment.yml --cores 4 --use-conda -n

# 특정 파일만 다시 생성 (예: 특정 pair의 volcano plot만)
snakemake --configfile configs/config_my_experiment.yml --cores 4 --use-conda \
  output/my_experiment/pairwise/TreatmentA_vs_Control/volcano_plot.png

# config가 바뀐 뒤 이미 완료된 작업까지 강제로 다시 실행하고 싶을 때
snakemake --configfile configs/config_my_experiment.yml --cores 4 --use-conda --forceall
```

> `--use-conda`를 붙이면 Snakemake가 각 rule 실행 시 `rna-seq-de-go-analysis` 환경을 자동으로 활성화합니다.
> `config_file`이 거의 모든 rule의 input으로 선언돼 있어서, config 내용을 수정하면(같은 파일이어도) 다음 실행 시 관련 단계가 자동으로 다시 계산됩니다 — 이미 완료된 프로젝트의 config를 고칠 때는 전체 재계산 시간을 감안하세요.

#### 여러 프로젝트를 한 번에 실행하려면: `run_batch.sh`

`configs/config_*.yml`로 등록된 프로젝트들을 순차적으로 한 번에 실행할 수 있습니다.

```bash
bash run_batch.sh                              # 모든 프로젝트 순차 실행
bash run_batch.sh --list                       # 실행 대상 프로젝트 목록만 확인
bash run_batch.sh --dry-run                    # 전체 실행 계획만 확인
bash run_batch.sh --include mouse-h2o2,kkj      # 이름에 해당 문자열 포함된 것만 실행
bash run_batch.sh --exclude 2026-hiy            # 해당 프로젝트 제외
bash run_batch.sh --force                       # 조건 변경 후 완료된 작업도 강제 재실행
bash run_batch.sh --cores 8                     # 프로젝트당 사용할 코어 수 지정
```

`count_data_path`가 존재하지 않는 프로젝트는 자동으로 건너뛰고, 결과는 `logs/batch/batch_<timestamp>.log`(전체 요약) + 프로젝트별 로그로 남습니다.

---

## 🔬 파이프라인 단계 및 스크립트 설명

전체 흐름은 `Snakefile`에 rule 단위로 정의되어 있습니다. 아래는 각 단계의 스크립트와 역할입니다.

### 1단계 — 차등발현 분석 (DE)
| 스크립트 | 역할 |
|---|---|
| `01a_run_omnibus_test.R` | (선택) 전체 그룹 간 차이를 보는 omnibus test(LRT/F-test) 실행 → `omnibus_test_results.csv` |
| `01b_run_pairwise_de.R` | `pairwise_comparisons`에 정의된 각 비교쌍에 대해 DESeq2/edgeR/limma-voom 실행 → `final_de_results.csv/xlsx` |
| `01c_run_masigpro_timeseries.R` | (선택) maSigPro 기반 time-series 회귀 + 발현 패턴 클러스터링 → `time_series/` |
| `09_export_multi_group.R` | (선택) omnibus 통계 + normalized count를 CMG-SeqViewer MULTI_GROUP 포맷으로 export |
| `10_run_coexpression_modules.R` | (선택) omnibus 유의 유전자 서브셋에 대해 `DEGreport::degPatterns` 클러스터링 → `coexpression_modules/` |

### 2단계 — QC 및 시각화
| 스크립트 | 역할 |
|---|---|
| `02_generate_plots.R` | Global PCA plot, pairwise volcano plot 생성 |
| `02a_generate_qc_plots.R` / `02c_generate_global_qc_report.R` | 전체 샘플 대상 QC 플롯 + HTML 리포트 |
| `02b_generate_pairwise_qc_plots.R` / `02d_generate_pairwise_qc_report.R` | 비교쌍별 QC 플롯(MA plot, heatmap 등) + HTML 리포트 |

### 3단계 — GO/KEGG Enrichment 및 심화 분석
| 스크립트 | 역할 |
|---|---|
| `03_enrichment_analysis.R` | `clusterProfiler::enrichGO()`/`enrichKEGG()` 실행 + dotplot. Gene Set(up/down)에 대해 term clustering(Jaccard)·GO Slim rollup·rrvgo 의미 축약도 이 안에서 함께 수행 |
| `04_generate_go_plots.R` | Ontology별(BP/CC/MF) 통합 GO bar plot 생성 |
| `05_generate_go_table.R` | `final_go_results.xlsx` 생성 (Gene Set × Ontology 워크시트) |
| `05b_generate_clustered_go_table.R` | `03`의 term clustering 결과를 CMG-SeqViewer Clustered-GO 포맷(`final_go_clustered_results.xlsx`)으로 집계 (전역 유일 `cluster_id` 부여) |
| `05c_generate_go_slim_overview.R` | `03`의 GO Slim 결과를 up/down 합쳐 ontology별 대칭 bar chart(`go_slim_overview_*.png`)로 집계 |
| `11_run_group_enrichment.R` | time-series 클러스터 / coexpression module처럼 그룹 컬럼으로 나뉜 유전자 리스트에 대해 그룹별 GO/KEGG enrichment 수행 |

### 4단계 — 리포트 및 배포
| 스크립트 | 역할 |
|---|---|
| `06_export_seqviewer.R` / `06b_aggregate_seqviewer.R` | pair별 CMG-SeqViewer parquet/staging 생성 후 프로젝트 전체 `metadata.json`으로 집계 |
| `07_generate_summary_report.R` | 전체 결과를 요약한 HTML 리포트 생성 |
| `08_generate_methods_section.R` | 사용된 통계 방법/파라미터를 논문 Methods 초안(Markdown)으로 자동 기록 |
| (Snakefile 내 `upload_to_gdrive` rule) | `rclone`으로 결과 폴더 전체를 Google Drive에 업로드 |

각 프로젝트는 `output/{output_dir}/config_used.yml`(및 각 하위 분석 폴더의 동일 파일)로 그 실행에 실제 사용된 config의 스냅샷을 함께 저장합니다.

---

## ⚙️ config.yml 상세 설명

전체 옵션과 기본값·주석은 **`configs/template/config.yml`이 항상 최신 소스입니다.** 여기서는 섹션별 개요만 정리합니다.

```yaml
species: "human"                 # "human" / "mouse" / "rat"
gene_id_type: "ENSEMBL"          # "ENSEMBL" / "ENTREZID" / "SYMBOL" (생략 시 자동 감지)

count_data_path: "data/raw/..."
metadata_path: "data/raw/..."
output_dir: "output/my_experiment"

de_analysis:
  method: "DESeq2"                       # "DESeq2" / "edgeR" / "limma-voom"
  design_formula: "~ group"
  group_variable: "group"
  run_omnibus_test: true                 # 전체 그룹 간 차이 검정(LRT/F-test)
  multi_group_export: {...}              # omnibus 결과를 CMG-SeqViewer MULTI_GROUP으로 export
  time_series: {...}                     # maSigPro 기반 time-series 분석 (include_groups 지원)
  coexpression_modules: {...}            # degPatterns 기반 coexpression module 분석 (include_groups 지원)
  pairwise_comparisons: [[비교군, 기준군], ...]
  padj_cutoff: 0.05
  log2fc_cutoff: 1.0
  advanced_options: {...}                # pre-filtering, 정규화 전략, VST, PCA 옵션

enrichment:
  gene_lists: ["total", "up", "down"]
  go_ontologies: ["BP", "CC", "MF"]
  pvalue_cutoff / qvalue_cutoff / min_gs_size / max_gs_size / min_gene_count / plot_top_n
  term_cluster: {...}                    # Jaccard 기반 GO term clustering (CMG Clustered-GO 포맷 export)
  go_slim: {...}                         # GO Slim rollup (up/down 대칭 bar chart)
  rrvgo: {...}                           # 의미론적 축약 (treemap/scatter, parentTerm 라벨)

qc_plots: {...}
databases: {...}                         # species별 organism DB / KEGG code (자동)
plot_aesthetics: {...}
go_barplot: {...}

export:
  export_to_excel: true
  seqviewer: true                        # CMG-SeqViewer parquet + staging JSON 생성

upload:
  gdrive: true
  remote: "gd_pargilbong:"
  dest_folder: "CMG_folder/Data/01_RNA-SEQ"
```

### 그룹이 여러 개일 때 서브셋 분석하기 (`include_groups`)

`time_series`와 `coexpression_modules` 둘 다 `include_groups`를 지원합니다 — metadata 전체 그룹 중 일부만 골라서 그 서브셋으로 시계열/coexpression 분석을 할 수 있습니다. 단, Snakemake 기본 슬롯은 프로젝트당 하나만 지원하므로, **같은 종류의 서브셋 분석을 두 개 이상** 돌리려면(예: acute 시계열 + chronic 시계열) 파생 config(`variant_label` 설정 포함)를 만들어 스크립트를 직접 실행해야 합니다. 실제 예시는 `docs/plan_time_series_coexpression_modules.md`를 참고하세요.

---

## 📊 주요 결과물

| 파일 | 생성 스크립트 | 설명 |
|---|---|---|
| `output/{project}/pairwise/{pair}/final_de_results.xlsx` | `01b` | 비교쌍별 DE 결과 |
| `output/{project}/pairwise/{pair}/final_go_results.xlsx` | `05` | GO/KEGG enrichment 결과, Gene Set(UP/DOWN/TOTAL) × Ontology(BP/CC/MF) 워크시트 |
| `output/{project}/pairwise/{pair}/final_go_clustered_results.xlsx` | `05b` | Jaccard 기반 클러스터링된 GO 결과 — CMG-SeqViewer가 바로 임포트 가능한 포맷(`cluster_id` 포함) |
| `output/{project}/pairwise/{pair}/go_slim_overview_{BP,CC,MF}.png` | `05c` | GO Slim 기반 up/down 대칭 bar chart |
| `output/{project}/pairwise/{pair}/go_rrvgo_*` | `03` | rrvgo treemap/scatter plot + 의미 축약 CSV |
| `output/{project}/multi_group_result.csv` | `09` | 전체 그룹 omnibus 통계 + normalized count |
| `output/{project}/time_series/`, `coexpression_modules/` | `01c`, `10` | 시계열/coexpression 분석 결과 + 클러스터별 GO/KEGG |
| `output/{project}/summary_report.html` | `07` | 전체 결과 요약 리포트 |
| `output/{project}/methods_section.md` | `08` | 논문 Methods 초안 |
| `output/{project}/seqviewer/metadata.json` + `datasets/*.parquet` | `06`, `06b` | CMG-SeqViewer 임포트용 데이터셋 |

`final_go_results.xlsx`의 워크시트는 `{GENESET}_{ONTOLOGY}` 형식(예: `UP_BP`, `DOWN_CC`)으로 구성되며, 별도의 통합 "All_Results"/"Top_Terms" 시트는 생성하지 않습니다(개별 시트로 충분히 탐색 가능하도록 설계).

---

## 📁 폴더 구조

```
RNA-Seq_DE_GO_analysis/
├── environment.yml             # ★ R + Snakemake + rclone 통합 conda 환경 정의
├── Snakefile                   # ★ 파이프라인 규칙 정의 (--configfile로 config 지정)
├── run_batch.sh                # 여러 프로젝트 순차 실행 스크립트
├── configs/
│   ├── template/config.yml     # ★ Git 추적 템플릿 (직접 수정 금지, 복사해서 사용)
│   └── config_*.yml            # 프로젝트별 설정 파일 (Git에서 자동 무시)
├── data/
│   └── raw/                    # 원본 카운트/메타데이터
├── output/                     # 프로젝트별 분석 결과 저장 폴더
├── src/
│   └── analysis/                   # 단계별 R 스크립트 (위 "파이프라인 단계" 표 참고)
├── bridge/                     # 🔗 RNA-Seq_GO_GSEA_analysis 파이프라인 연동용 변환 스크립트
├── notebooks/
│   └── analysis_pipeline.ipynb # (참고용, Snakemake 사용을 권장)
└── logs/batch/                 # run_batch.sh 실행 로그
```

---

## 🔗 Pipeline 연결: Advanced GSEA 분석으로 확장하기

이 파이프라인으로 DE 분석과 기본 GO enrichment를 완료한 후, 더 심화된 GSEA(Gene Set Enrichment Analysis)가 필요하다면 `bridge/` 디렉토리의 변환 도구로 **RNA-Seq_GO_GSEA_analysis** 파이프라인과 연결할 수 있습니다. 자세한 내용은 `bridge/README.md`를 참고하세요.

```bash
cd bridge
python3 convert_de_to_gsea.py \
  --input ../output/{project}/pairwise/{pair}/final_de_results.csv \
  --output ../../RNA-Seq_GO_GSEA_analysis/data/{pair}.xlsx
```

---

## ☁️ Google Drive 자동 업로드

파이프라인이 완료되면(`summary_report.html` + `methods_section.md` 생성 시점) `rclone`을 통해 결과 폴더 전체를 Google Drive에 자동으로 업로드할 수 있습니다. `rclone copy`를 사용하므로 동일한 파일은 건너뛰고, 재실행해도 안전합니다.

### 사전 준비 (최초 1회)

`rclone`은 `environment.yml`에 포함되어 있어 conda 환경 생성 시 함께 설치됩니다. Google Drive remote만 별도로 설정하면 됩니다.

```bash
conda activate rna-seq-de-go-analysis
rclone config
```

대화형 설정: `n`(새 remote) → 이름 입력(예: `my_gdrive`) → 스토리지 타입 `drive` → Client ID/Secret은 비워두고 Enter → Scope `1` → 브라우저 인증.

```bash
rclone listremotes   # 예: my_gdrive:
```

### config 설정

```yaml
upload:
  gdrive: true
  remote: "my_gdrive:"
  # 실제 업로드 경로: {remote}{dest_folder}/{output_dir 폴더명}/
  dest_folder: "RNA-Seq_Results"
```

`gdrive: false`로 두면 업로드 단계를 건너뜁니다. 완료 표시는 `output/{project}/.gdrive_upload_done.flag`, 로그는 `output/{project}/logs/10_upload_to_gdrive.log`에 남습니다.

---

## DE 분석 알고리즘 비교: DESeq2 vs edgeR vs limma-voom

| 구분 | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) | [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) | [limma-voom](https://bioconductor.org/packages/release/bioc/html/limma.html) |
|:------|:--------|:--------|:-------------|
| **핵심 통계 모델** | 음이항 분포, 엄격한 분산 추정 | 음이항 분포, 유연하고 빠른 분산 추정 | 선형 모델 + voom 분산 정규화 |
| **장점** | 위양성 억제에 효과적, 직관적인 사용법, 자동 정규화 | 빠른 속도, 복잡한 디자인 대응, 메모리 효율 | 초고속, 복잡한 디자인(다중 요인/반복측정)에 강함 |
| **단점** | 샘플이 많으면 느려짐, 다소 보수적 | DESeq2보다 위양성 가능성 약간 높음 | voom 전처리 필수, 극단적 저발현 유전자 처리 까다로움 |
| **추천 상황** | 표준 분석, 신뢰도 중요, 그룹당 3~5개 샘플 | 빠른 분석, 복잡한 디자인 | 매우 큰 데이터셋, 다중 요인 디자인 |

### 참고문헌
- **DESeq2**: Love, M.I., Huber, W., Anders, S. (2014). *Genome Biology*, 15:550. https://doi.org/10.1186/s13059-014-0550-8
- **edgeR**: Robinson, M.D., McCarthy, D.J., Smyth, G.K. (2010). *Bioinformatics*, 26(1):139-140. https://doi.org/10.1093/bioinformatics/btp616
- **limma**: Ritchie, M.E., Phipson, B., Wu, D., et al. (2015). *Nucleic Acids Research*, 43(7):e47. https://doi.org/10.1093/nar/gkv007

---

## 🔬 GO/KEGG 분석: clusterProfiler + 심화 분석

GO/KEGG enrichment는 [**clusterProfiler**](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)의 `enrichGO()`/`enrichKEGG()`를 사용한 **ORA(Over-Representation Analysis)** 방식입니다:

1. **입력**: 통계적으로 유의한 차등 발현 유전자 목록 (예: padj < 0.05, |log2FC| > 1)
2. **통계 검정**: Fisher's exact test — "이 유전자 목록에 특정 GO term이 우연보다 많이 포함되어 있는가?"
3. **보정**: Benjamini-Hochberg FDR

### 유의 term이 너무 많을 때: 3가지 심화 분석

FDR/fold-enrichment로 걸러도 수백 개 term이 남는 경우가 흔해서, 이 파이프라인은 세 가지 보조 분석을 함께 제공합니다(모두 `enrichment` 섹션에서 개별 on/off):

| 방법 | 기반 | 산출물 | 용도 |
|---|---|---|---|
| **Term clustering** (`term_cluster`) | 유전자 중복도(Jaccard) | `final_go_clustered_results.xlsx`, treeplot | CMG-SeqViewer에서 클러스터 단위로 탐색 |
| **GO Slim** (`go_slim`) | GO DAG level | `go_slim_overview_*.png` | 실험 전체 방향을 bar chart 한 장으로 파악 |
| **rrvgo** (`rrvgo`) | GO DAG 의미 거리 | treemap/scatter plot | 의미상 유사한 term을 사람이 읽기 쉬운 상위 개념으로 묶어서 조망 |

### ORA vs GSEA

| 특징 | ORA (본 파이프라인) | GSEA |
|:-----|:-------------------|:-----|
| **입력** | 선별된 유의 유전자 목록 | 전체 유전자의 발현 순위 |
| **장점** | 간단하고 직관적, 빠름 | 약한 신호도 감지 가능 |
| **단점** | Cutoff에 민감함 | 계산 복잡, 해석 어려움 |

GSEA가 필요하면 위 [Pipeline 연결](#-pipeline-연결-advanced-gsea-분석으로-확장하기) 섹션의 `bridge/` 도구를 사용하세요.

### 참고문헌
- **clusterProfiler**: Yu, G., Wang, L.G., Han, Y., He, Q.Y. (2012). *OMICS*, 16(5):284-287.
- **clusterProfiler 4.0**: Wu, T., Hu, E., Xu, S., et al. (2021). *The Innovation*, 2(3):100141.
- **Gene Ontology**: The Gene Ontology Consortium (2021). *Nucleic Acids Research*, 49(D1):D325-D334.
- **KEGG**: Kanehisa, M., Goto, S. (2000). *Nucleic Acids Research*, 28(1):27-30.
- **rrvgo**: Sayols, S. (2023). Bioconductor. https://bioconductor.org/packages/rrvgo/

---

## 🧬 Gene ID 타입 설정

`gene_id_type`으로 입력 데이터의 ID 형식을 지정합니다.

| ID 타입 | 설명 | 예시 |
|:--------|:-----|:-----|
| **ENSEMBL** | Ensembl 유전자 ID | `ENSG00000141510`(human), `ENSMUSG00000051951`(mouse) |
| **ENTREZID** | NCBI Entrez 유전자 ID | `7157`, `5594` |
| **SYMBOL** | 유전자 심볼 | `TP53`, `GAPDH` |

생략 시 첫 번째 gene ID 패턴으로 자동 감지되지만(`ENSMUSG...`→mouse ENSEMBL, `ENSG...`→human ENSEMBL, 숫자만→ENTREZID, 그 외→SYMBOL), 명시적 지정을 권장합니다. GO/KEGG 분석은 내부적으로 Entrez ID가 필요하므로 ENSEMBL/SYMBOL 입력 시 자동 변환됩니다(변환 실패 유전자는 자동 제외, 보통 80~95% 성공).

**Raw count 파일 형식 예시**:
```csv
Geneid,Sample1,Sample2,Sample3
ENSMUSG00000051951,1234,5678,9012
ENSMUSG00000021803,234,567,890
```
버전 번호(`ENSMUSG00000051951.4`)는 제거하는 것을 권장합니다.

---

## 🧹 Pre-filtering: 저발현 유전자 제거

`de_analysis.advanced_options`의 두 값으로 필터링 기준을 정합니다:

```yaml
de_analysis:
  advanced_options:
    prefilter_min_count: 1        # 이 count 이상인 것을 "발현됨"으로 침
    prefilter_min_samples: "auto" # 최소 몇 개 샘플에서 발현돼야 유지할지
                                   # "auto" → 비교 그룹 중 가장 작은 그룹 크기를 자동 사용
```

즉 "`prefilter_min_count` 이상 발현된 샘플이 `prefilter_min_samples`개 이상인 유전자만 유지"하는 방식입니다. 전체 샘플 합계 기준이 아니라 **샘플 수 기준**이라 그룹 크기가 다른 다중 그룹 비교에서도 생물학적으로 더 타당합니다.

### 왜 저발현 유전자를 제거해야 할까요?
1. **통계적 검정력 향상**: 검정하는 유전자 수가 줄면 FDR 보정 부담이 줄어 실제 유의 유전자를 더 잘 찾을 수 있음
2. **분산 추정 정확도 향상**: 노이즈성 저발현 유전자가 분산 추정에 끼는 영향 제거
3. **계산 효율성**: 분석 시간·메모리 절감

### 주의사항
- 너무 엄격하면(예: min_samples를 그룹 크기보다 크게) 생물학적으로 중요한 저발현 유전자(전사인자 등)를 놓칠 수 있습니다.
- 필터링은 DE 분석 전 단계이며, 그룹별 발현 차이를 미리 예단하지 않고 "발현 여부"만 봅니다.

### 참고 문헌
- DESeq2 공식 vignette: "Pre-filtering the dataset" — https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering
- edgeR User's Guide: "Filtering"

---

## 🔄 Pairwise 비교 시 정규화 전략

다중 그룹 실험에서 pairwise 비교 시 정규화 방법에 따라 결과가 달라질 수 있습니다. `de_analysis.advanced_options.pairwise_normalization`으로 선택합니다.

### 1. `"subset"` — 각 비교마다 해당 두 그룹의 샘플만으로 정규화
- ✅ 독립적인 비교, 더 보수적인 결과, 명확한 해석
- ⚠️ 비교마다 정규화 기준이 달라 서로 직접 비교하기 어려움, 전체 맥락 손실

### 2. `"global"` — 전체 샘플로 정규화한 뒤 pairwise 비교 (템플릿 기본값)
- ✅ 모든 비교에 일관된 정규화, 전체 실험 맥락 유지, 서로 다른 비교 결과를 직접 비교 가능, 더 안정적인 정규화 계수 추정
- ⚠️ 비교에 포함되지 않은 그룹의 발현 패턴도 정규화에 영향, subset 대비 민감도가 약간 높아질 수 있음

**추천**: 여러 pairwise 비교를 종합적으로 보거나(시계열·용량반응 실험 등), 그룹당 샘플 수가 적어 정규화 안정성이 중요할 때는 `"global"`을, 각 비교를 완전히 독립적으로 해석하고 싶을 때는 `"subset"`을 사용하세요.

```yaml
de_analysis:
  advanced_options:
    pairwise_normalization: "global"  # 또는 "subset"
```

---

## 💡 주요 개념 (FAQ)

**모든 샘플에서 raw count가 0인 유전자는 삭제해도 되나요?**
네, 권장합니다. 삭제하지 않으면 다중 검정 보정이 과도하게 엄격해져 실제 유의 유전자를 놓칠 수 있습니다. `prefilter_min_count: 1`, `prefilter_min_samples: "auto"`가 기본값입니다.

**Normalized Counts와 FPKM의 차이는?**
DE 분석은 동일 유전자를 샘플 간 비교하므로 유전자 길이 보정이 불필요합니다. DESeq2/edgeR은 라이브러리 크기(시퀀싱 깊이)만 보정한 normalized count를 사용합니다.

**Ensembl ID는 있는데 Gene Symbol이 비어있는 이유는?**
오류가 아닙니다. Non-coding gene, 아직 명명되지 않은 novel gene 등은 공식 Symbol이 없을 수 있습니다.

**pvalue_cutoff와 qvalue_cutoff의 차이는?**
- **P-value cutoff**: "이 GO term 농축이 우연일 확률"에 대한 느슨한 1차 필터
- **Q-value cutoff**: 다중 검정 보정(FDR) 후 엄격한 최종 필터

**Dot plot의 Count와 GeneRatio는?**
- **Count**(점 크기): 내 유전자 목록 중 해당 GO term에 포함되는 유전자 개수
- **GeneRatio**(x축): Count / 분석에 사용한 전체 유전자 수

**GeneRatio와 GSEA의 Enrichment Score는 같은 개념인가요?**
아닙니다. GeneRatio는 ORA 방식(본 파이프라인)의 지표이고, Enrichment Score는 GSEA 방식의 지표입니다. 이 파이프라인은 ORA만 사용하므로 Enrichment Score는 계산되지 않습니다.

**GeneRatio도 유의성 cutoff가 있나요?**
아니요. GeneRatio는 유의성이 아니라 영향력의 크기를 나타내는 척도입니다. `padj`로 먼저 유의한 term을 고르고, GeneRatio는 그 안에서 "얼마나 많은 유전자가 관여하는지" 해석에 참고합니다.
