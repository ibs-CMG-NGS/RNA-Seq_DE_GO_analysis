# maSigPro 시계열 분석 + Omnibus-필터 기반 Coexpression Module 분석 도입

> **구현 현황 (2026-07-22)**: Phase 0~8 전체 구현 및 GSE205486 실데이터 end-to-end 테스트 완료.
> 테스트 중 발견한 실제 버그: maSigPro의 dispersion(theta) 추정을 시간 다항식(연속형) design으로 하면
> 잔차 자유도 부족으로 dispersion이 비정상적으로 크게 잡혀(theta≈0.27) GLM이 `NA/NaN/Inf` 에러로
> 죽는 문제가 있어, Time×Series 실험 셀을 나타내는 factor 모델로 dispersion을 재추정하도록 수정함
> (theta≈1.27로 안정화). 상세 내용은 `src/analysis/01c_run_masigpro_timeseries.R`의
> "Dispersion(theta) 외부 추정" 주석 참고. WGCNA/CEMiTool 정식 도입은 계획대로 미착수 상태(Future Phase).

## Context

현재 파이프라인은 그룹 간 차이를 보는 omnibus test(LRT/F-test)와 pairwise DE만 지원하고, (1) 시간 경과에 따른 발현 궤적(trajectory)을 모델링하는 시계열 분석과 (2) 유전자들이 함께 움직이는 상관 구조(coexpression module)를 찾는 분석이 없다. 대화에서 확인한 사용자 상황(그룹당 n=3, 적은 time point, 3그룹 이상 비교)에 맞춰 두 가지를 추가한다:

1. **maSigPro** — count 기반 회귀로 적은 replicate/time point에서도 안정적으로 시계열 유의 유전자를 찾고 패턴별로 클러스터링까지 해주는 도구.
2. **Coexpression module 분석** — omnibus test로 걸러진 유의 유전자 서브셋에 한해 발현 패턴 기반 클러스터링(**DEGreport::degPatterns**)을 수행. 전체 transcriptome 규모의 WGCNA/CEMiTool은 지금 샘플 규모(n=3 × 그룹)에서 불안정하므로 **의도적으로 이번 범위에서 제외**하고 향후 별도 단계로 보류한다.

두 모듈 모두 기존 파이프라인 컨벤션(Snakemake rule 패턴, R 스크립트 하우스 스타일, config 섹션 구조)을 그대로 따르고, **결과를 CMG-SeqViewer로 반입 가능한 형태(parquet + staging JSON)로도 내보낸다** — 이 파이프라인의 다른 모든 산출물(`final_de_results`, `final_go_results`, `multi_group_result`)이 이미 그렇게 하고 있고, 사용자가 실제로 CMG-SeqViewer에서 탐색/시각화를 하기 때문에 신기능도 동일한 경로로 소비 가능해야 한다.

## 핵심 설계 결정 (조사로 확정된 사항)

- **메타데이터 스키마 확장 필요**: 현재 예시 메타데이터(`GSE142445`, `GSE205486` 등)는 `sample_id, condition, replicate`만 있고 숫자형 `time` 컬럼이 없다. `condition` 문자열("4hr_ipsilateral", "Lesion_D1")을 정규식으로 자동 파싱하는 대신, **명시적 `time`(숫자) + 선택적 `series`(범주형, 예: 부위/유전자형) 컬럼을 config에서 지정하도록 요구**한다 — 기존 `group_variable` 방식과 동일한 철학(암묵적 추론보다 명시적 선언).
- **maSigPro는 count 모드로 사용**: `library(MASS)` 명시적 로드 필요(안 하면 `negative.binomial()` not found 에러), `theta`(분산 파라미터)는 자동 추정되지 않으므로 `edgeR::estimateGLMCommonDisp()`로 외부 추정 후 전달. `degree`는 `n_timepoints - 1`을 넘을 수 없고, series별로 **최소 3개의 서로 다른 time point**가 있어야 함 — 스크립트에서 사전 검증 후 명확한 에러로 중단.
- **Coexpression은 hand-rolled hclust 대신 `DEGreport::degPatterns` 채택**: replicate를 그룹 평균으로 정리한 뒤 Kendall correlation + `cluster::diana`로 클러스터링해 n=3 노이즈에 더 강하고, 범주형 그룹(시간이 아니어도 됨)에 바로 적용 가능. 단, 기본 `minc=15`(최소 클러스터 크기)가 유의 유전자가 수십 개 수준일 때 조용히 클러스터를 병합/폐기하므로 **config에서 명시적으로 낮춰 노출**해야 함(예: 3~5).
- **CMG-SeqViewer 연동은 기존 컨벤션 그대로 확장**: `dataset_type` 필드는 파이프라인마다 자유롭게 새 값을 정의해서 써온 전례가 있음(RNA-seq: `differential_expression`/`go_analysis`/`multi_group`, ATAC-seq: `differential_accessibility` — sibling 저장소 `atac-seq-da-analysis/ATAC_PIPELINE_SEQVIEWER_PACKAGING.md` 확인). 앱 쪽 코드 수정 없이 이 저장소 안에서 새 `dataset_type`("time_series", "coexpression_module")을 추가하는 것이 기존 관행과 일치. 단, 같은 문서에서 지적된 실수(조건/조직/tags 필드를 비워둠)를 반복하지 않도록 메타데이터 필드를 처음부터 꼼꼼히 채운다.
- WGCNA/CEMiTool은 이번 계획에 설계 언급만 하고 **구현하지 않음** (Future Phase로 분리, 트리거 조건 명시).

## 재사용할 기존 코드/패턴

- `src/utils/load_data.R`의 `create_de_object()` — config 기반 `DESeqDataSet` 생성 (VST 변환 전 단계까지).
- `01b_run_pairwise_de.R`의 하우스 스타일: positional args(`commandArgs`), `%||%` null-coalescing 헬퍼, `stop(paste(...))` 방식 에러 처리, `cat()` 기반 로깅, 끝에 `config_used.yml` 복사, `openxlsx`로 xlsx export(`export.export_to_excel` 플래그로 게이팅).
- `02b_generate_pairwise_qc_plots.R:164-192`의 heatmap 패턴: `t(scale(t(mat)))` z-score, `colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)`, `png()`/`dev.off()` 래핑 — coexpression 모듈 heatmap에 그대로 재사용.
- `09_export_multi_group.R` 전체 — **가장 중요한 참고 대상**. omnibus 통계 + normalized count를 합쳐 CSV/parquet/staging JSON을 "한 스크립트 안에서" 만드는 패턴이 이번 두 모듈(시리즈/omnibus 결과처럼 pair 단위가 아닌 프로젝트 단위 산출물)과 정확히 같은 모양이라, 그대로 복제해서 쓴다: `make_alias_slug()`, `new_uuid()`, `write_parquet_dataset()` 헬퍼, `seqviewer/datasets/`·`seqviewer/staging/` 디렉토리 생성, `enabled=false`면 `cat(...); quit(save="no", status=0)`로 조용히 스킵.
- `06_export_seqviewer.R`의 staging entry 필드 스키마(`dataset_id, alias, original_filename, dataset_type, experiment_condition, organism, cell_type, tissue, timepoint, row_count, gene_count, significant_genes, import_date, file_path, notes, tags`) — 두 신규 모듈의 staging entry도 이 필드셋을 그대로 채운다(특히 `timepoint`는 지금까지 항상 빈 문자열이었는데, time_series 모듈에서는 실제로 채울 수 있는 첫 사례가 됨).
- `06b_aggregate_seqviewer.R`은 `staging/*_entries.json`을 전부 글롭(glob)해서 병합하므로 **수정 불필요** — 단, Snakemake DAG는 명시적 `input`으로만 의존성을 인식하므로 `Snakefile`의 `aggregate_seqviewer` rule에 새 staging JSON 경로를 `input`으로 추가해야 함(글롭만으로는 실행 순서가 보장되지 않음).

---

## Phase 0 — 설계 확정 (config 스키마 최종화)
**내용**: 아래 config 블록을 `configs/template/config.yml`에 확정 (seqviewer export 옵션 포함).
```yaml
de_analysis:
  time_series:
    enabled: false
    time_variable: "time"       # 메타데이터의 숫자형 time 컬럼명 (신규 컬럼 요구)
    series_variable: null       # 여러 시리즈 비교 시 범주형 컬럼명 (null=단일 시리즈)
    degree: "auto"              # "auto" = min(n_timepoints-1, 2), 정수 지정 가능
    q_value: 0.05
    rsq_cutoff: 0.6
    pattern_k: 6                # see.genes 클러스터 개수
    export_seqviewer: true      # parquet + staging JSON 생성 여부

  coexpression_modules:
    enabled: false
    padj_cutoff: null           # null이면 de_analysis.padj_cutoff 재사용
    min_genes: 10                # 이보다 적으면 스킵(경고 후 종료)
    min_cluster_size: 5          # degPatterns minc 오버라이드
    export_seqviewer: true       # parquet + staging JSON 생성 여부
```
**비용**: 1~2시간. **난이도**: 낮음 — 이번 대화에서 대부분 확정됨, 문서화만 남음.

## Phase 1 — 환경 셋업
**내용**: `environment.yml`에 `bioconductor-masigpro`, `bioconductor-degreport` 추가. 겸사겸사 **기존에 누락돼 있던 `r-arrow`, `r-jsonlite`도 추가**(`09_export_multi_group.R`/`06_export_seqviewer.R`가 이미 `library(arrow)`/`library(jsonlite)`를 쓰는데 `environment.yml`에는 선언이 안 되어 있어 실제 설치된 환경과 파일이 drift된 상태 — 이번에 같이 바로잡음). conda 환경 재구축 후 `library(maSigPro)`, `library(DEGreport)`, `library(MASS)`, `library(arrow)`, `library(jsonlite)` 로딩 확인.
**비용**: 실작업 1~2시간 + conda dependency solve/설치 대기(환경에 따라 수십 분~1시간). **난이도**: 낮음 — 패키지 추가 자체는 단순하지만 conda solve 시간이 변수.

## Phase 2 — maSigPro 시계열 모듈 구현
**내용**: `src/analysis/01c_run_masigpro_timeseries.R` 신규 작성 + `Snakefile`에 `run_masigpro_timeseries` rule 추가(`time_series.enabled`일 때만 `rule all`에 편입).
- 메타데이터에서 `time_variable`(+ `series_variable`) 존재 검증, series별 고유 time point ≥3 검증 (미달 시 어떤 series가 문제인지 명시하고 `stop()`)
- `edesign` 구성: `Time`, `Replicates`(동일 Time×Series 조합당 동일 값), series별 0/1 더미 컬럼
- `edgeR::estimateGLMCommonDisp()`로 `theta` 외부 추정 → `p.vector(..., counts=TRUE, family=negative.binomial(theta), theta=theta)` → `T.fit()` → `get.siggenes()` → `see.genes()`로 클러스터 패턴 산출
- 출력: 유의 유전자 CSV(gene_id, p-value, R², cluster_id), 패턴 플롯 PNG, `config_used.yml` 복사, xlsx export(하우스 스타일 준수)
**비용**: 6~10시간. **난이도**: 높음 — count 모드 API가 까다로움(MASS 명시적 로드 누락, theta 미추정, edesign Replicates 의미 오해가 흔한 실수 포인트로 확인됨). 실제 R 세션에서 반복 디버깅 필요.

## Phase 3 — maSigPro 결과 리포트 통합
**내용**: `07_generate_summary_report.R`, `08_generate_methods_section.R`에 `time_series.enabled` 조건부 섹션 추가(maSigPro 인용: Conesa et al. 2006 포함).
**비용**: 3~4시간. **난이도**: 중간 — `multi_group_export`가 이미 만들어둔 "조건부 섹션 포함" 패턴을 따라가면 되지만, 두 리포트 스크립트 내부 구조를 파악하는 시간이 필요.

## Phase 4 — maSigPro → CMG-SeqViewer 연동
**내용**: `01c_run_masigpro_timeseries.R` 안에 `export_seqviewer: true`일 때 실행되는 export 블록 추가 (`09_export_multi_group.R` 패턴 그대로 복제).
- 유의 유전자 × (gene_symbol, p-value, R², cluster_id, time로 정렬된 샘플별 normalized/VST 값) 통합 데이터프레임 구성
- `write_parquet_dataset()`으로 `seqviewer/datasets/*.parquet` 저장
- staging entry: `dataset_type="time_series"`, `experiment_condition`에 series 설명, **`timepoint` 필드에 실제 time 값 목록 기록**(기존 파이프라인들이 늘 비워뒀던 필드를 처음으로 채우는 사례), `tags=c("time_series","maSigPro", ...)`, `notes`에 degree/Q-value 명시
- `seqviewer/staging/time_series_entries.json`으로 저장
- `Snakefile`: `rule all`에 조건부 output 추가 + **`aggregate_seqviewer` rule의 `input`에 이 staging JSON 경로를 명시적으로 추가**(글롭만으로는 DAG 의존성이 안 잡힘 — `06b_aggregate_seqviewer.R` 자체는 수정 불필요)
**비용**: 3~4시간. **난이도**: 중간 — 패턴 복제라 로직 자체는 쉽지만, Snakemake DAG에 새 staging 파일을 두 군데(`rule all`, `aggregate_seqviewer` input)에 빠짐없이 연결해야 하는 실수 포인트가 있음.

## Phase 5 — Coexpression module 분석 구현
**내용**: `src/analysis/10_run_coexpression_modules.R` 신규 작성 + `Snakefile`에 `run_coexpression_modules` rule 추가(`omnibus_test_results.csv`에 의존, `coexpression_modules.enabled`일 때만 편입).
- omnibus CSV에서 `padj_cutoff` 기준 유의 유전자 필터 → `min_genes` 미달 시 경고 후 정상 종료(스킵)
- `create_de_object()` → `vst()` → 유의 유전자로 서브셋
- `DEGreport::degPatterns(ma=..., metadata=meta, time=group_variable, minc=min_cluster_size)` 호출
- 출력: gene→module 배정 CSV, degPatterns 기본 제공 패턴 플롯, xlsx export, `config_used.yml` 복사
**비용**: 3~5시간. **난이도**: 중간 — `degPatterns` 채택으로 난이도는 낮아졌지만 `minc` 기본값(15)을 명시적으로 낮추지 않으면 조용히 클러스터가 사라지는 리스크를 코드/문서 양쪽에서 방어해야 함.

## Phase 6 — Coexpression 결과 리포트 통합
**내용**: 기존 heatmap 패턴(z-score + RdBu + png/dev.off) 재사용한 module heatmap 추가, `07_generate_summary_report.R`/`08_generate_methods_section.R`에 조건부 섹션 반영(DEGreport 인용 포함).
**비용**: 2~3시간. **난이도**: 낮음 — `degPatterns`가 자체 ggplot 패턴 플롯을 제공해 추가 시각화 부담이 적음.

## Phase 7 — Coexpression → CMG-SeqViewer 연동
**내용**: `10_run_coexpression_modules.R` 안에 `export_seqviewer: true`일 때 실행되는 export 블록 추가 (Phase 4와 동일 패턴).
- gene_id, gene_symbol, module_id, omnibus 통계(padj 등), 그룹별로 정렬된 샘플별 normalized/VST 값 통합
- staging entry: `dataset_type="coexpression_module"`, `tags=c("coexpression","module", 그룹명...)`, `notes`에 `minc`/상관 방법 명시
- `seqviewer/staging/coexpression_modules_entries.json`으로 저장
- `Snakefile`: Phase 4와 동일하게 `rule all` + `aggregate_seqviewer` input 양쪽에 연결
**비용**: 2~3시간. **난이도**: 낮음 — Phase 4에서 만든 패턴을 재사용.

## Phase 8 — 통합 테스트 & 문서화
**내용**: 실제 시계열성 데이터(`GSE142445-stroke-timecourse-mouse`, `GSE205486-spinalcord-mouse`)의 메타데이터 **사본**에 `time`(필요시 `series`) 컬럼을 추가해 end-to-end 테스트 config 구성 후 두 모듈 모두 실행 검증. `configs/template/config.yml` 주석 정비, README/pipeline-integration 문서 갱신. **`seqviewer/metadata.json`에 새 `dataset_type` 항목이 정상 병합되는지, parquet를 실제로 로드해서 컬럼/값이 맞는지 확인** (CMG-SeqViewer 앱 자체 UI 테스트는 이 저장소 범위 밖이므로 데이터 계약까지만 검증).
**비용**: 4~6시간. **난이도**: 중간 — n=3 규모에서 나오는 경고/엣지케이스(예: 특정 series의 time point 부족) 대응이 실제로 발생할 가능성이 높음.

**총 예상 비용**: 약 25~37시간 (약 4~6 작업일 상당).

---

## Future Phase (보류) — WGCNA / CEMiTool 정식 도입

지금은 구현하지 않고 아래 트리거 조건 중 하나가 충족되면 별도 계획으로 재검토:
- 전체 샘플 수가 충분히 늘어나(대략 15~20개 이상) 안정적인 네트워크 구축이 가능해질 때
- 유의 유전자 서브셋이 아니라 전체 transcriptome 규모의 module discovery가 필요해질 때
- Module별 자동 GO/pathway enrichment 연계(예: CEMiTool의 내장 ORA)가 요구될 때

예상 난이도: 높음(soft-threshold power 선택, module merging 파라미터 튜닝 등), 예상 비용: 10시간+ — 별도 논의 필요. 도입 시에도 이번에 만든 `write_parquet_dataset()`/staging JSON 패턴을 그대로 재사용해 CMG-SeqViewer 연동 가능.

---

## 변경/신규 파일 요약

| 파일 | 변경 내용 |
|---|---|
| `src/analysis/01c_run_masigpro_timeseries.R` | 신규 (분석 + seqviewer export 포함) |
| `src/analysis/10_run_coexpression_modules.R` | 신규 (분석 + seqviewer export 포함) |
| `Snakefile` | rule 2개 추가 + `rule all` 조건부 output 추가 + `aggregate_seqviewer` input 확장 |
| `environment.yml` | `bioconductor-masigpro`, `bioconductor-degreport`, (누락 보정) `r-arrow`, `r-jsonlite` 추가 |
| `configs/template/config.yml` | `time_series`, `coexpression_modules` 섹션 추가 (`export_seqviewer` 옵션 포함) |
| `src/analysis/07_generate_summary_report.R` | 조건부 섹션 추가 |
| `src/analysis/08_generate_methods_section.R` | 조건부 섹션/인용 추가 |
| (테스트용) 메타데이터 CSV 사본 | Phase 8에서 `time` 컬럼 추가한 복사본 생성 (원본 미변경) |

## 검증 방법

1. Phase 1 후: `conda run -n rna-seq-de-go-analysis Rscript -e 'library(maSigPro); library(DEGreport); library(MASS); library(arrow); library(jsonlite)'`로 로딩 확인.
2. Phase 2 후: GSE205486(time point 5개: D1/D2/D3/D5/D7)의 테스트용 메타데이터로 `time_series.enabled: true` 설정 후 `snakemake`로 `01c` rule 단독 실행 → 유의 유전자 CSV/패턴 플롯 생성 및 로그에 에러 없는지 확인.
3. Phase 4 후: 같은 실행에서 `seqviewer/datasets/*.parquet`와 `seqviewer/staging/time_series_entries.json`이 생성되는지, `arrow::read_parquet()`로 다시 읽었을 때 컬럼/행 수가 기대와 일치하는지 확인.
4. Phase 5 후: `run_omnibus_test: true` + `coexpression_modules.enabled: true`로 전체 파이프라인 실행 → module CSV/heatmap 생성 확인, `min_genes` 미달 케이스(예: padj_cutoff를 매우 낮게)로 스킵 로직도 별도 확인.
5. Phase 7 후: `seqviewer/staging/coexpression_modules_entries.json` 생성 및 parquet 재로드 검증.
6. Phase 8 후: 두 모듈 + 두 seqviewer export 모두 활성화한 상태로 전체 `snakemake` 실행이 끝까지 완료되는지, `06b_aggregate_seqviewer.R`가 만든 `seqviewer/metadata.json`에 `time_series`/`coexpression_module` 타입 항목이 기존 항목과 함께 정상 병합되는지, `summary_report.html`/`methods_section.md`에 새 섹션이 반영되는지 확인.

## 실제 테스트 결과 (Phase 8, 2026-07-22)

- 데이터셋: `GSE205486-spinalcord-mouse` (29 샘플, Lesion/Naive × D1/D2/D3/D5/D7), 메타데이터 사본에 `time`(1/2/3/5/7)·`series`(Lesion/Naive) 컬럼 추가해 테스트.
- maSigPro: `filterByExpr`(Time×Series 셀 기준)로 57180→39529 유전자 필터링 후 정상 완료, 879개 유의 유전자(R²≥0.6) 검출. CSV/xlsx/패턴 플롯/parquet/staging JSON 모두 생성 및 재검증 완료.
- Coexpression modules: omnibus 유의 유전자가 12211개로 매우 많아(강한 처치 효과 데이터셋) `degPatterns` 클러스터링이 상당히 느려짐 — `padj_cutoff`를 더 엄격하게(1e-6, 2737개) 잡아 검증했고, 2730개 유전자가 39개 모듈로 클러스터링됨. **참고**: 유의 유전자가 수천~수만 개로 매우 많은 프로젝트에서는 `coexpression_modules.padj_cutoff`를 기본값(전역 padj_cutoff)보다 엄격하게 오버라이드하는 것을 권장.
- `summary_report.html`, `methods_section.md` 모두 신규 섹션(Time-Series Analysis, Coexpression Modules / 2a, 2b) 정상 렌더링 확인.
