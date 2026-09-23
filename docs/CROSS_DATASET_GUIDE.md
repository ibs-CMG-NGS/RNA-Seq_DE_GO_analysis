# Cross-Dataset 비교 가이드

이미 파이프라인을 완주한 프로젝트 2개 이상을 사후에 비교해서 insight를 뽑는
독립 도구다. 한 프로젝트 안의 pairwise 조건 비교(`enrichment.cross_condition`,
`12_run_cross_condition_comparison.R`)와는 다른 축 — 여기서는 **프로젝트(데이터셋)
자체**가 비교 단위다. 같은 종의 RNA-seq 프로젝트 2개, 다른 종의 RNA-seq
프로젝트, 심지어 RNA-seq 프로젝트와 [atac-seq-da-analysis](../../atac-seq-da-analysis)
저장소의 ATAC-seq 프로젝트까지도 하나의 config로 함께 비교할 수 있다.

## 왜 별도 도구인가

- Snakemake `configfile:`은 프로젝트 1개짜리 `config` 변수 구조라 여러
  `output_dir`을 동시에 참조하는 이 작업을 자연스럽게 표현하지 못한다.
- "이미 끝난 프로젝트 2개를 사람이 판단해서 비교"하는 성격이라 파일 변경 감지
  기반 자동 재실행(Snakemake의 핵심 가치)이 필요 없다.
- 그래서 `configs/cross_dataset_*.yaml`이라는 별도 config 파일 + 독립 스크립트
  2개(`Rscript ... <config>`)로 직접 실행한다. 어떤 프로젝트의 `rule all`에도
  편입되지 않는다.

## 스크립트 2개

| 스크립트 | 비교 대상 | 원형(포팅 출처) |
|---|---|---|
| `18_run_cross_dataset_go_comparison.R` | pairwise DE/DA 비교(조건 vs Control)의 GO/KEGG term — common/flip/exclusive/mixed + rrvgo 의미론적 축약 + heatmap/dotplot/UpSet/pair 산점도 | atac-seq-da-analysis `20_run_cross_species_go_comparison.R` |
| `19_run_cluster_cross_dataset_comparison.R` | time-series 클러스터 / coexpression 모듈의 GO term-set Jaccard 유사도(다대다 매칭) | atac-seq-da-analysis `28_run_cluster_cross_species_comparison.R`(N-데이터셋으로 일반화) |

## config 스키마

`configs/template/cross_dataset_template.yaml`이 주석 포함 전체 템플릿이다.
핵심 필드:

```yaml
output_dir: "output/cross_dataset_LABEL_A-vs-LABEL_B"

datasets:
  - label: "mouse_rna"                                # groups/flips/exclusives의 축
    config: "configs/config_2026-06-hiy-mouse-rna.yml" # 그 프로젝트의 1개짜리 config
    assay: "rna"                                       # "rna"(기본값) | "atac"
    pair_map:                                           # (선택) 정식 이름 -> 실제 폴더명
      "Acute_1D_vs_Control": "D1_vs_Control"
  - label: "human_rna"
    config: "configs/config_2026-06-hiy-human-rna.yml"
    assay: "rna"

pairs:                                                  # 18번만 사용(19번은 불필요)
  - "Acute_1D_vs_Control"
  - "Acute_3D_vs_Control"

go_ontologies: [BP, KEGG]
fdr_cutoff: 0.05
fold_enrichment_cutoff: 2.0
```

### `assay`: rna | atac

RNA-Seq_DE_GO_analysis와 atac-seq-da-analysis는 GO enrichment 표 자체(컬럼)는
동일하지만 그 표를 어디서 찾는지·geneID가 어떤 포맷인지가 다르다(직접 대조
확인, [`docs/atac_pipeline_alignment_request.md`](atac_pipeline_alignment_request.md)에
차이 전체 정리 및 장기적으로 이 차이를 없애자는 정합 제안). `assay` 필드로
데이터셋마다 다음 프리셋 중 하나를 선택한다(18/19번 스크립트 안에 내장):

| | `assay: rna`(기본값) | `assay: atac` |
|---|---|---|
| pairwise GO/KEGG | `pairwise/{pair}/enrichment/go_enrichment_{dir}_{ont}.csv` | `pairwise/{pair}/go_enrichment_{dir}_{ont}.csv`(서브폴더 없음) |
| geneID 포맷 | Entrez ID (자동으로 SYMBOL 매핑 후 비교) | Gene SYMBOL (그대로 사용) |
| time-series 클러스터 GO | `time_series[_{variant}]/go_termcluster_cluster{N}_BP.csv` | `time_series/{variant}/go_enrichment/masigpro_cluster_{k}_BP.csv` |
| coexpression 모듈 GO | `coexpression_modules[_{variant}]/go_termcluster_module{N}_BP.csv` | `coexpression_modules/[{variant}/]go_enrichment/module_{id}_BP.csv` |

atac-seq-da-analysis 프로젝트를 참조할 때는 `config:`에 **그 레포의 절대경로**를
쓴다(두 레포가 `~/ngs-pipeline/` 하위 형제 디렉토리라는 전제) — 상대경로를 써도
동작은 하지만(각 데이터셋 config 파일이 위치한 레포 루트 기준으로 자동 resolve),
가독성을 위해 절대경로를 권장한다.

```yaml
  - label: "mouse_atac"
    config: "/home/ngs/ngs-pipeline/atac-seq-da-analysis/configs/config_2025-ljh-hiy-mouse.yaml"
    assay: "atac"
```

### `pair_map`

같은 비교(예: "급성 1일차 vs 대조군")인데 프로젝트마다 pairwise 폴더명이 다를 수
있다 — `pairs:`에는 항상 정식(canonical) 이름을 쓰고, 데이터셋별로 실제 폴더명이
다르면 `pair_map`으로 매핑한다. Condition id(그래프 축 라벨, common/flip 판정
키)는 항상 정식 이름을 쓰고 파일 경로 생성에만 매핑된 이름이 쓰인다.

## Worked Example 3종 (실제 검증에 쓰인 config, 그대로 재실행 가능)

### 1. RNA vs RNA (다른 종) — `configs/cross_dataset_hiy-mouse-vs-human-rna.yaml`

```bash
Rscript src/analysis/18_run_cross_dataset_go_comparison.R configs/cross_dataset_hiy-mouse-vs-human-rna.yaml
Rscript src/analysis/19_run_cluster_cross_dataset_comparison.R configs/cross_dataset_hiy-mouse-vs-human-rna.yaml
```

`flip_DOWN_mouse_to_UP_human_BP.csv`의 최상위 term이 `chromosome segregation`,
`nuclear division` 등 세포주기 유전자 — mouse는 억제·human은 활성화되는 방향으로
갈린다는, H2O2 자극에 대한 종간 반응 차이의 핵심 신호.

### 2. ATAC vs ATAC (다른 종) — `configs/cross_dataset_hiy-mouse-vs-human-atac.yaml`

```bash
Rscript src/analysis/19_run_cluster_cross_dataset_comparison.R configs/cross_dataset_hiy-mouse-vs-human-atac.yaml
```

atac-seq-da-analysis 자신의 `20_run_cross_species_go_comparison.R`/
`28_run_cluster_cross_species_comparison.R`이 이미 만들어둔
`output/cross_species_hiy-mouse-vs-human/`의 Jaccard 행렬과 값이 정확히
일치함을 확인(포팅 검증).

### 3. RNA vs ATAC (같은 종) — `configs/cross_dataset_hiy-mouse-rna-vs-atac.yaml`

```bash
Rscript src/analysis/18_run_cross_dataset_go_comparison.R configs/cross_dataset_hiy-mouse-rna-vs-atac.yaml
Rscript src/analysis/19_run_cluster_cross_dataset_comparison.R configs/cross_dataset_hiy-mouse-rna-vs-atac.yaml
```

지금까지 어느 레포에도 없었던 새로운 종류의 비교 — "이 유전자가 전사체 수준에서
움직이는 시점/패턴과 그 근처 크로마틴이 열리고 닫히는 시점/패턴이 얼마나
겹치는가"를 GO term-set Jaccard로 정량화한다.

## 산출물

`18_run_cross_dataset_go_comparison.R` (output_dir 바로 아래):

| 파일 | 내용 |
|---|---|
| `common_{up,down}_{loose,strict}_{ont}.csv` | 여러 데이터셋에서 공통으로 같은 방향인 term. loose=최소 데이터셋 수 충족, strict=반대 방향 등장 없음 |
| `common_{up,down}_strict_rrvgo_{ont}.csv` | strict 결과에 rrvgo 의미론적 축약 적용(GO만, KEGG 제외) |
| `flip_{UP,DOWN}_{A}_to_{B}_{ont}.csv` | 데이터셋 A에서 UP/DOWN, B에서 반대 방향인 term + gene Jaccard(low_overlap 플래그 포함) |
| `exclusive_{label}_{ont}.csv` | 특정 데이터셋(들)에만 있고 나머지엔 전혀 없는 term |
| `mixed_{label}_{ont}.csv` | 같은 데이터셋 안에서 조건마다 방향이 혼재하는 term |
| `pair_scatter_{pair}_{A}_vs_{B}_{ont}.png`/`_data.csv` | pair별 두 데이터셋 직접 산점도(concordant/discordant/데이터셋-특이 색 구분) |
| `cross_dataset_{heatmap,dotplot}_{ont}.png` | 전체 후보 term heatmap/dot plot |
| `upset_{UP,DOWN}_{ont}.png`/`_membership_*.csv` | (dataset::pair)별 유의 term 중첩 구조 |
| `final_cross_dataset_go_results.xlsx` | 위 카테고리 전부를 시트별로 취합 |
| `condition_count_log.txt` | 조건×방향별 유의 term 수 로그(sparsity 경고 포함) |

`19_run_cluster_cross_dataset_comparison.R` (데이터셋 쌍마다):

| 파일 | 내용 |
|---|---|
| `ts_cluster_{A}_vs_{B}_jaccard_matrix.csv` / `coexpr_module_{A}_vs_{B}_jaccard_matrix.csv` | 모든 클러스터/모듈 조합의 GO term-set Jaccard 전체 행렬 |
| `..._best_match.csv` | 클러스터/모듈별 최고 매칭 상대 + 공유 term 수 |
| `..._jaccard_heatmap.png` | 위 행렬 히트맵 |
| `cluster_cross_dataset_provenance.csv` | 소스 파일 패턴·처리 방식 기록 |

## 트러블슈팅

- **모든 조건/데이터셋에서 `0 significant terms`가 찍히면 십중팔구 경로 문제다.**
  `assay` 값이 그 데이터셋의 실제 구조와 맞는지, `pair_map`이 필요한데 빠지지
  않았는지부터 확인. `condition_count_log.txt`(18번) 또는 콘솔 로그(19번)에서
  어느 (데이터셋, pair/클러스터) 조합이 비어있는지 바로 확인 가능.
- **다른 레포(ATAC) 데이터셋을 참조했는데 경로가 안 맞으면**: 그 데이터셋
  `config:`에 적힌 project config의 `output_dir`이 상대경로인 경우, "그 config
  파일이 들어있는 레포 루트" 기준으로 자동 resolve된다(예:
  `atac-seq-da-analysis/configs/config_X.yaml` 안의 `output_dir: "output/Y"`는
  `atac-seq-da-analysis/output/Y`로 풀린다) — 이 규칙과 다르게 배치된 프로젝트라면
  `config:`에 절대경로를 직접 쓸 것.
- **`flip_*`/`common_*` 결과가 의심스러울 정도로 적/많으면**: `sparsity_warn_threshold`
  미만인 (데이터셋, 조건, 방향) 조합이 있는지 로그를 확인 — 유의 term이 원래
  적은 비교가 섞이면 common/flip 판정이 왜곡될 수 있다는 경고.
- **gene-level Jaccard(`flip_*.csv`의 `Gene Jaccard` 컬럼)는 정식 ortholog 매핑이
  아니다** — 대소문자 정규화(`toupper()`)한 심볼 문자열 일치일 뿐이므로 참고용
  지표로만 사용할 것.
