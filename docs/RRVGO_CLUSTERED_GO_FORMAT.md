# RRVGO Clustered-GO Format (`final_go_rrvgo_clustered_results.xlsx`)

각 pairwise 비교 폴더(`output/{project}/pairwise/{pair}/`)에 `05d_generate_rrvgo_clustered_go_table.R`이 생성하는 파일의 포맷 스펙입니다. 같은 폴더의 `final_go_clustered_results.xlsx`(`docs/UPSTREAM_CLUSTERED_GO_FORMAT.md`, cmg-seqviewer 레포 소유)와 **의도적으로 같은 컬럼 계약을 공유**하되, 클러스터링 알고리즘이 다르다는 것을 구분할 수 있도록 컬럼 2개를 추가했습니다.

이 문서는 cmg-seqviewer 쪽 임포터가 아직 구현되지 않은 상태에서, 이 저장소가 먼저 확정한 산출물 포맷을 기록해두는 용도입니다. 임포터 작업 시 그대로 참고하면 됩니다.

## 두 클러스터드-GO 파일의 차이

| | `final_go_clustered_results.xlsx` | `final_go_rrvgo_clustered_results.xlsx` |
|---|---|---|
| 생성 스크립트 | `05b_generate_clustered_go_table.R` | `05d_generate_rrvgo_clustered_go_table.R` |
| 클러스터링 기준 | Term 간 유전자 중복도 (Jaccard, `clusterProfiler::pairwise_termsim(method="JC")`) | GO DAG 의미 유사도 (`rrvgo::calculateSimMatrix(method="Rel")` + `reduceSimMatrix()`) |
| "묶임"의 의미 | 유의 유전자 집합이 실제로 겹치는 term끼리 | GO 트리 상에서 개념적으로 가까운 term끼리 (유전자가 안 겹쳐도 묶일 수 있음) |
| 활성화 조건 | `config.enrichment.term_cluster.enabled` (기본 `true`) | `config.enrichment.rrvgo.enabled` (기본 `true`) |
| 대상 gene set | up / down (total 제외) | up / down (total 제외) |
| 최소 term 수 | 3개 이상 | 2개 이상 |

같은 GO ID가 두 파일 모두에 등장하더라도, `cluster_id`로 묶인 동료 term은 서로 다를 수 있습니다 — 별개의 분석 결과이지 하나를 다른 하나로 대체하는 관계가 아닙니다.

## 시트 구성

`final_go_clustered_results.xlsx`와 동일하게, gene set × ontology 조합마다 시트 하나:
`UP_BP`, `UP_CC`, `UP_MF`, `DOWN_BP`, `DOWN_CC`, `DOWN_MF` (최대 6개, `enrichment.go_ontologies` 설정에 따라 줄어들 수 있음). 유의 term이 하나도 없는 프로젝트는 `No Results` 단일 시트의 placeholder 파일이 생성됩니다.

## 컬럼 스펙

| 컬럼명 | 타입 | 설명 |
|---|---|---|
| `Gene Set` | string | `UP` 또는 `DOWN` |
| `Ontology` | string | `BP` / `CC` / `MF` |
| `GO ID` | string | GO term ID (예: `GO:0022836`) |
| `GO Term` | string | GO term 설명 |
| `Gene Ratio` | string | `n/N` 형식 (clusterProfiler 원본 그대로) |
| `Background Ratio` | string | `n/N` 형식 |
| `P-value` | double | raw p-value |
| `Adjusted P-value` | double | BH-FDR 보정값. **컬럼명은 항상 `Adjusted P-value`이며 `padj`가 아님** (기존 스펙과 동일 규칙) |
| `Gene Count` | int | 해당 term에 속한 유의 유전자 개수 |
| `Gene Symbols` | string | `/`로 구분된 gene symbol 목록 (Entrez ID → Symbol 변환) |
| `cluster_id` | string | 그룹 내(gene set × ontology) 전역 유일 zero-padded 3자리 문자열(`"001"`, `"002"`, ...). 최종 클러스터 크기가 1인 term은 `"Singleton"` |
| `Representative Term` | string | **(rrvgo 전용 추가 컬럼)** 해당 클러스터를 대표하는 term 이름 (`rrvgo::reduceSimMatrix()`의 `parentTerm`) — score(기본 `-log10(padj)`)가 가장 높은 term이 대표로 선택됨 |
| `Algorithm` | string | **(rrvgo 전용 추가 컬럼)** 고정값 `"Semantic (rrvgo)"`. 두 클러스터드-GO 파일을 한 화면에서 다룰 경우 클러스터링 방식을 구분하는 용도 |

## `cluster_id` 부여 규칙

`final_go_clustered_results.xlsx`와 동일한 규칙입니다:
- 파일 전체(모든 시트)에 걸쳐 유일한 3자리 zero-padded 번호(`"001"`, `"002"`, ...)
- 클러스터 크기가 1인 term(다른 어떤 term과도 안 묶인 경우)은 번호 대신 `"Singleton"`
- 정렬 순서: `cluster_id`(Singleton은 맨 뒤) → `Adjusted P-value` 오름차순

## 데이터 소스 (구현 참고용)

`go_rrvgo_{gene_set}_{ont}.csv`(rrvgo의 `reducedTerms` 그대로, `go/cluster/parent/score/size/term/parentTerm/...` 컬럼)에는 `Gene Ratio`/`Background Ratio`/`P-value`/`Gene Count`/`Gene Symbols`가 없습니다. `05d`는 이를 같은 폴더의 `go_enrichment_{gene_set}_{ont}.csv`(원본 `enrichResult` 전체 덤프, 모든 gene set에 대해 항상 생성됨)와 GO ID 기준으로 join해서 채웁니다.
