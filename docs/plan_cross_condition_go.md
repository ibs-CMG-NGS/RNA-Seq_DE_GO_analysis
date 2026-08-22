# GO Term Cross-Condition 비교 분석 계획서

**적용 대상:** mRNA-seq 파이프라인에서 다수 조건(≥ 2)의 GO enrichment 결과를 보유한 모든 프로젝트  
**온톨로지:** GO BP (primary) · KEGG (primary) · MF (optional) · CC (선택적)

---

## 목차

1. [배경 및 목적](#1-배경-및-목적)
2. [기존 파이프라인 접점](#2-기존-파이프라인-접점)
   - [기존 분석과 중첩 여부](#기존-분석과-중첩-여부)
   - [입력 소스 결정](#입력-소스-결정-raw-go-enrichment-vs-rrvgo-결과)
3. [분석 범위 — 온톨로지 선택](#3-분석-범위--온톨로지-선택)
4. [핵심 알고리즘](#4-핵심-알고리즘)
   - [A. 공통 Term](#a-공통-term--find_common_direction)
   - [B. 방향 역전](#b-방향-역전--find_direction_flip)
   - [C. 조건 특이적 Term](#c-조건-특이적-term--find_exclusive--find_mixed)
   - [D. Semantic Filtering](#d-semantic-filtering-rrvgo-재적용)
5. [주의 사항 (Known Issues)](#5-주의-사항-known-issues)
6. [추천 시각화](#6-추천-시각화)
7. [Phase별 구축 계획](#7-phase별-구축-계획)
8. [프로젝트 적용 체크리스트](#8-프로젝트-적용-체크리스트)

---

## 1. 배경 및 목적

표준 mRNA-seq 파이프라인(DEG → GO enrichment → GO clustering)은 각 조건을 **독립적으로** 기술한다. 조건별 dot plot과 treemap은 "이 조건에서 무슨 경로가 활성화됐는가"에 답하지만, 연구자가 실제로 알고 싶은 질문에는 직접 답하지 않는다.

**이 분석이 필요한 이유:**
- 조건별 GO dot plot을 나란히 놓고 눈으로 비교하는 수작업을 자동화한다.
- 입력은 이미 생성된 GO enrichment CSV이므로 추가 연산 비용이 없다.
- 조건 수가 많을수록(≥ 3) 수작업 비교의 한계가 커지고 이 분석의 가치가 높아진다.

**이 분석이 답하는 질문:**
- 어떤 GO term이 여러 조건에서 일관되게 상방/하방 조절되는가?
- 특정 조건 그룹(예: 초기 vs 후기, 처리 vs 대조)에서 방향이 역전되는 경로는 무엇인가?
- 특정 조건에만 나타나는 term은 무엇인가?
- 조건 간 GO term 중첩 구조(어떤 조건 조합에서 공유되는가)는 어떠한가?

---

## 2. 기존 파이프라인 접점

파이프라인 구조에서 Cross-Condition 비교 분석의 위치:

```
1단계  DEG 탐색           (per condition)
2단계  GO Enrichment      (per condition × direction)  ← 입력 소스
         예: go_bp_up_cluster_dot_bundle/inputs/data.csv
3단계  GO Clustering      (per condition × direction, rrvgo)
4단계  Cross-Condition 비교  ← 신규 (모든 조건 완료 후 synthesis)
         ├─ 공통 UP/DOWN (여러 조건에서 일관된 term)
         ├─ 방향 역전 (그룹 A에서 UP → 그룹 B에서 DOWN, 또는 반대)
         └─ 조건 특이적 term (특정 조건에만 존재하거나 특정 조건에서만 사라지는 term)
4b단계 Semantic Filtering + 시각화
         비교 결과에 rrvgo 재적용 → multi-condition dot plot
```

### 기존 분석과 중첩 여부

아래 표는 일반적인 time-course 또는 multi-condition mRNA-seq 파이프라인을 기준으로 한다. 프로젝트마다 포함된 분석 모듈이 다를 수 있으므로 실제 파이프라인 구성에 맞춰 확인한다.

| 질문 | Time Series 분석 | Co-expression Module 분석 | Cross-Condition GO |
|---|---|---|---|
| 조건 내 시간/그룹별 GO | ✓ 커버 가능 | — | 일부 겹침 |
| 두 조건 그룹 공통 GO | — | ✓ 커버 가능 | 일부 겹침 |
| 전체 조건 공통 GO | — | 제한적 | **새로 커버** |
| GO term 수준 방향 역전 | — | — | **새로 커버** |
| 조건 조합별 중첩 구조 | — | — | **새로 커버** |
| 전체 합성 시각화 | — | — | **새로 커버** |

> **고려사항:** co-expression module GO 분석이 이미 있는 경우 결과를 나란히 비교해 중복 정도를 확인한다. 완전히 독립적인 분석이 아니므로 해석 시 중복을 명시한다.

### 입력 소스 결정: Raw GO Enrichment vs rrvgo 결과

> **결정: 2단계 raw GO enrichment CSV를 입력으로 사용하고, semantic filtering은 출력 단에서 별도 적용한다.**

rrvgo 결과를 입력으로 쓰면 조건마다 다른 representative term이 선택될 수 있어 term_id 기준 집합 연산이 깨진다. 예컨대 조건 A에서는 `GO:0030198`이, 조건 B에서는 `GO:0043062`가 대표로 선택되면 같은 의미 클러스터가 "비공통"으로 분류된다.

Raw GO enrichment에서 term이 빠지지 않게 비교한 뒤, 결과 리스트에 rrvgo를 재적용해 중복을 정리하는 순서가 안전하다.

---

## 3. 분석 범위 — 온톨로지 선택

| 온톨로지 | 우선순위 | Cross-condition에서의 가치 | 중복 수준 |
|---|---|---|---|
| GO BP | **Primary** | 생물학적 과정 전체 커버. 해석 가치 최대. | 높음 → rrvgo 필수 |
| KEGG | **Primary** | pathway 수가 적고 큐레이션이 잘 되어 결과가 깔끔. BP를 actionable pathway로 보완. | 낮음 → filtering 최소 |
| GO MF | Optional | 분자 활성 수준 메커니즘 보완. 특정 연구 질문(효소 활성, 결합 등)이 있을 때 포함. | 중간 |
| GO CC | 낮음 | 세포 내 위치 정보. 특정 compartment 가설 검증 목적일 때만 포함. | 낮음 |

**온톨로지 선택 가이드:**
- 연구 질문이 "어떤 생물학적 과정이 변하는가"이면 BP + KEGG로 충분하다.
- 효소 활성, 분자 기능 수준의 메커니즘이 중요하면 MF를 추가한다.
- 세포 소기관/구획 변화가 핵심 가설이면 CC를 포함한다.

**Term 수 규모 기준 시각화 선택:**  
비교 결과 term 수가 50개 이하이면 기존 dot plot 포맷(조건을 열로 확장)이 gene ratio 정보까지 담아 더 정보량이 많다. 50개 초과이면 heatmap이 패턴 파악에 유리하다. Semantic filtering 후 대부분의 경우 50개 이하로 줄어든다.

---

## 4. 핵심 알고리즘

모든 함수는 순수 함수(side-effect 없음)로 설계한다. 입력은 term database(`term_id → {조건: {dir, fdr, fe, gene_count, genes}}`), 출력은 GoTerm 리스트.

**공통 파라미터:**
- `fdr_cutoff`: term 로딩 시 FDR 필터 임계값 (기본값 0.05)
- `min_conditions`: 공통 term 인정 최소 조건 수 (기본값: 전체 조건 수 × 0.75 반올림)
- `conditions`: 분석에 포함할 조건 목록 (config.yaml에서 지정)
- `groups`: 조건을 묶는 상위 그룹 정의 (방향 역전 분석에 사용)

### A. 공통 Term — `find_common_direction()`

N개 이상의 조건에서 동일 방향으로 유의미한 term을 추출한다.

```
common_UP   = { t | count(조건 c: t ∈ UP(c)) ≥ min_conditions }
common_DOWN = { t | count(조건 c: t ∈ DOWN(c)) ≥ min_conditions }

strict mode: 모든 entries가 동일 방향이어야 함 (혼재 term 제외)
             → 보수적 기준. 논문 main figure 용도.
loose mode:  방향 카운트만 기준 (다른 방향 혼재 허용)
             → 탐색적 분석 또는 조건 수가 많을 때.
```

- **정렬 기준:** 등장 조건 수 내림차순 → 최소 FDR 오름차순
- **활용:** 처리 조건 전체에서 일관되게 조절되는 핵심 경로 파악

### B. 방향 역전 — `find_direction_flip()`

그룹 A에서 X 방향, 그룹 B에서 반대 방향인 term. 두 그룹 내에서 혼재 없는 **깨끗한 역전**만 통과시킨다.

```
flip_UP→DOWN (group_A → group_B) =
  ( ⋃ UP(c) for c ∈ group_A )
  ∩ ( ⋃ DOWN(c) for c ∈ group_B )
  − ( ⋃ DOWN(c) for c ∈ group_A )   ← group_A 내 DOWN 혼재 제거
  − ( ⋃ UP(c) for c ∈ group_B )     ← group_B 내 UP 혼재 제거

flip_DOWN→UP (group_A → group_B) = 위의 방향 교환
```

- **groups 설정 예시:**
  - Time-course: `{early: [t1, t2], late: [t3, t4]}`
  - Dose-response: `{low_dose: [0.1uM, 1uM], high_dose: [10uM, 100uM]}`
  - Treatment vs control: `{control: [ctrl], treated: [treat_1h, treat_6h]}`

> **⚠ 해석 주의:** 같은 GO term이더라도 두 그룹에서 **다른 유전자 집합**이 구동할 수 있다. "방향 역전"이 아니라 "다른 유전자가 우연히 같은 GO term에 걸린 것"일 수 있으므로, 각 조건의 `gene_symbols`를 함께 출력하고 Jaccard overlap을 확인해야 한다.

### C. 조건 특이적 Term — `find_exclusive()` / `find_mixed()`

**`find_exclusive()`:** target 조건 전부 포함 & absent 조건 전부 미포함인 term.

```
exclusive(target, absent, direction=None) =
  { t | ∀c ∈ target: t ∈ cond(c) [and t.direction == direction if specified] }
  ∩ { t | ∀c ∈ absent: t ∉ cond(c) }
```

설정 예시 (config.yaml):

```yaml
exclusives:
  - label: early_only_up          # 초기 조건에만 UP
    target: [cond_A, cond_B]
    absent: [cond_C, cond_D]
    direction: UP

  - label: late_shared_down       # 후기 조건들에서 공통 DOWN
    target: [cond_C, cond_D]
    absent: [cond_A, cond_B]
    direction: DOWN
```

**`find_mixed()`:** 같은 그룹 내에서 UP과 DOWN이 혼재하는 term. 진동 패턴(oscillation) 또는 비선형적 반응 포착.

### D. Semantic Filtering (rrvgo 재적용)

비교 결과 term 리스트에 `rrvgo::reduceSimMatrix()`를 적용해 의미적 중복을 제거한다. 이 단계의 rrvgo는 3단계(조건별 treemap)와 **별도**로 실행되며, 비교 결과만을 대상으로 한다.

```
입력: cross-condition 비교 결과 term_id 목록 + 해당 term들의 best_fdr
설정: threshold=0.7 (조정 가능), orgdb: 종에 따라 선택
      - 마우스: org.Mm.eg.db
      - 인간:   org.Hs.eg.db
      - 제브라피쉬: org.Dr.eg.db
      - 초파리: org.Dm.eg.db
출력: representative term 목록 (중복 제거된 GO term)
```

R 환경 의존성이 생기므로 Python 파이프라인에서는 rpy2 또는 별도 R script 호출로 처리한다. R script 분리 방식이 환경 의존성 측면에서 더 안정적이다.

---

## 5. 주의 사항 (Known Issues)

### ⚠ 특정 조건의 GO term 수 희소성

일부 조건에서 enriched term 수가 매우 적을 경우(예: UP term이 10개 미만) 공통 term 분석 결과가 왜곡된다. 이것이 생물학적 사실인지, 실험 설계/배치 효과인지, 필터링 임계값 문제인지를 먼저 확인해야 한다.

**대응 방법:**
- 파이프라인에 조건별 term count를 자동 로깅하고, 설정 임계값(예: 20개) 이하이면 경고 플래그를 출력한다.
- `min_conditions_common`을 낮추거나 `fdr_cutoff`를 완화(0.1)해 탐색적 분석을 먼저 수행한다.
- 희소성이 있는 조건을 제외한 subset 분석을 병행한다.

### ⚠ 방향 역전 term의 유전자 집합 확인 필요

역전 term 리스트만으로는 진짜 생물학적 역전인지 판단할 수 없다. 출력에 각 조건의 `gene_symbols`를 항상 포함시키고, 두 그룹 간 gene Jaccard overlap을 계산한다.

```
Jaccard(A, B) = |genes_A ∩ genes_B| / |genes_A ∪ genes_B|

Jaccard < 0.3  → "유전자 집합 불일치" 플래그: 생물학적 역전보다 다른 유전자 구동 가능성
Jaccard ≥ 0.3  → 같은 유전자들이 방향을 바꾼 진짜 역전 가능성 높음
```

### ⚠ rrvgo의 조건별 representative term 불일치

3단계 rrvgo 결과에서 각 조건의 representative term이 다를 수 있다. 이 문서의 입력 소스 결정(§2)이 이 문제의 근본적 해결책이다. 이미 rrvgo 결과만 보유한 경우에는 parent term(상위 의미 클러스터)으로 매핑 후 비교하는 대안을 사용한다.

### ⚠ 조건 수와 엄격도 균형

조건이 2개이면 모든 term이 "공통" 또는 "exclusive"가 되어 구분 의미가 없다. 조건이 5개 이상이면 `min_conditions`를 전체 수로 설정하면 공통 term이 거의 없어진다. 프로젝트별로 `min_conditions_common` 파라미터를 조정하고, strict/loose 두 모드를 모두 보고한다.

---

## 6. 추천 시각화

기본 원칙: 기존 파이프라인의 dot plot 포맷과 시각 언어를 통일한다. 새 포맷을 도입하는 것은 term 수나 조건 수가 많아 dot plot이 한계에 달할 때로 한정한다.

### Primary: Multi-condition Dot Plot (GO BP)

행=GO term, 열=조건. 점 크기=gene ratio, 색=방향×−log₁₀(FDR) 발산 스케일(UP 녹색, DOWN 적색). 기존 포맷을 열 방향으로 확장한 것이므로 파이프라인 일관성이 유지된다.

- term ≤ 50일 때 권장
- gene ratio 정보 포함이 장점
- 방향 역전 term에 별도 마커(★ 등) 추가로 구분 가능

### Primary: Multi-condition Bar Chart (KEGG)

기존 KEGG bar chart를 조건별 grouped bar로 확장. 각 pathway의 FDR 또는 gene ratio를 조건별로 나란히 표시. 방향(UP/DOWN)은 색으로 구분.

- KEGG는 term 수가 적어 bar가 적합
- pathway 간 상대적 강도 비교에 유리

### Secondary: UpSet Plot

조건 간 GO term 중첩 구조를 표시. Venn은 조건이 4개 이상이면 가독성이 떨어지지만 UpSet은 깔끔하다. 어떤 조건 조합에서 term이 공유되는지 한눈에 파악.

- Python: `upsetplot` 패키지
- R: ComplexHeatmap의 UpSet 함수

### 대규모 term 대응: Direction Heatmap

semantic filtering 후에도 term이 50개 초과이거나 조건이 6개 이상으로 늘어날 때. dot plot이 빽빽해지면 heatmap으로 전환. gene ratio 정보는 잃는다.

---

## 7. Phase별 구축 계획

### Phase 1 — 핵심 로직 구현 ✦ 난이도: 낮음

**목표:** "어떤 BP term이 여러 조건에서 일관되게 UP/DOWN 조절되는가"를 Python 코드로 안정적으로 산출한다.

**구현 항목:**
- `models.py`: `ConditionEntry`, `GoTerm`, `AnalysisConfig` 데이터 모델
- `loader.py`: raw GO enrichment CSV 로딩, 조건별 term_id set 구성, fdr_cutoff 필터, 컬럼명 매핑(col_map)
- `analyzer.find_common_direction()`: strict/loose 모드, min_conditions 파라미터화
- `reporter.py`: 조건별 FDR/FE/gene_count 컬럼 포함 CSV 및 Excel 출력
- 조건별 term count 자동 로깅, 희소성 경고 플래그
- `config.yaml` 기반 경로/조건명 설정 (하드코딩 없음)

**config.yaml 최소 구성:**
```yaml
project:
  name: my-project
  root: /path/to/data
  output_dir: ./output

conditions: [cond_A, cond_B, cond_C, cond_D]

directions: [UP, DOWN]

file_template:
  UP:   "{condition}/go_bp_up_cluster_dot_bundle/inputs/data.csv"
  DOWN: "{condition}/go_bp_down_cluster_dot_bundle/inputs/data.csv"

thresholds:
  fdr_cutoff: 0.05
  min_conditions_common: 3
```

**검증 방법:**
- 수작업 교집합 결과와 일치 여부 확인 (소규모 term으로 cross-check)
- 조건별 term count 로그가 기대값과 일치하는지 확인

**산출물:** `common_up_strict.csv`, `common_down_strict.csv`, `condition_count_log.txt`

---

### Phase 2 — 방향 역전 분석 + Gene Overlap 검증 ✦ 난이도: 중간

**목표:** 두 조건 그룹 사이에 방향이 역전되는 GO term을 추출하고, 유전자 집합 overlap을 함께 제공해 해석 신뢰성을 높인다.

**구현 항목:**
- `find_direction_flip()`: 집합 연산 기반 깨끗한 역전 탐지
- `find_exclusive()`: 특정 조건에만 존재하는 term (config의 `exclusives` 섹션으로 유연하게 정의)
- `find_mixed()`: 그룹 내 진동 패턴 term
- 각 역전 term에 대해 두 그룹 간 **gene Jaccard overlap** 계산 및 출력
- Jaccard < 0.3인 경우 `low_overlap` 플래그 추가

**config.yaml 추가 구성:**
```yaml
groups:
  group_A: [cond_A, cond_B]
  group_B: [cond_C, cond_D]

analysis:
  flips:
    - [group_A, group_B, UP]    # group_A UP → group_B DOWN
    - [group_A, group_B, DOWN]  # group_A DOWN → group_B UP

  exclusives:
    - label: group_A_only_up
      target: [cond_A, cond_B]
      absent: [cond_C, cond_D]
      direction: UP
```

**기술적 고려사항:** 유전자 목록이 구분자(예: `/`, `,`)로 연결된 문자열로 저장되어 있으므로, split 후 set 변환 로직이 필요하다. 종(species)에 따른 유전자 심볼 대소문자 일관성 확인 필요.

**산출물:** `flip_UP_groupA_to_DN_groupB.csv`, `flip_DN_groupA_to_UP_groupB.csv`, `exclusive_*.csv`, `gene_overlap_summary.csv`

---

### Phase 3 — 시각화 구현 ✦ 난이도: 중간

**목표:** Phase 1–2 결과를 논문 figure 수준의 시각화로 변환한다. 기존 파이프라인의 시각 언어와 통일한다.

**구현 항목:**
- **Multi-condition dot plot**: matplotlib/seaborn 기반. 조건=열, term=행, 점 크기=gene ratio, 색=방향(UP/DOWN) × -log₁₀(FDR) 발산 스케일. 기존 dot bundle과 같은 color convention 사용.
- **UpSet plot**: `upsetplot` 패키지. 조건별 UP/DOWN term 중첩 구조.
- 방향 역전 term 하이라이트: dot plot에 역전 term에 별도 마커 추가.
- 출력 포맷: SVG + PNG (300 dpi). 논문 figure용 크기 설정 파라미터화.

**기술적 고려사항:**
- 발산 색상 스케일 구현 시 중앙(유의하지 않음/absent)을 명확히 구분해야 한다.
- 조건에 term이 없는 경우(absent)와 FDR > threshold인 경우를 시각적으로 다르게 처리한다.
- 조건 수가 많아지면(≥ 6) figure 너비를 자동 조정하는 로직 필요.

**산출물:** `common_up_dotplot.svg`, `flip_dotplot.svg`, `upset_bp.svg`

---

### Phase 4 — KEGG 통합 + Semantic Filtering ✦ 난이도: 중간–높음

**목표:** KEGG를 BP와 병행 분석하고, BP 결과에 rrvgo semantic filtering을 적용해 최종 결과물을 압축한다.

**구현 항목:**
- **KEGG 통합**: config.yaml에 KEGG file template 추가. 동일 로직으로 분석. KEGG pathway ID 체계(hsa/mmu 등 species prefix) 처리.
- **MF 선택적 포함**: config에 on/off 플래그.
- **Semantic filtering**: Python에서 R `rrvgo` 호출 (rpy2 또는 subprocess). orgdb는 config에서 지정. threshold 파라미터화(기본 0.7).
- KEGG는 term 수가 적어 semantic filtering 선택적 적용.

**config.yaml KEGG 추가 예시:**
```yaml
file_template:
  UP:         "{condition}/go_bp_up_cluster_dot_bundle/inputs/data.csv"
  DOWN:       "{condition}/go_bp_down_cluster_dot_bundle/inputs/data.csv"
  KEGG_UP:    "{condition}/kegg_up_cluster_dot_bundle/inputs/data.csv"
  KEGG_DOWN:  "{condition}/kegg_down_cluster_dot_bundle/inputs/data.csv"

semantic_filtering:
  enabled: true
  threshold: 0.7
  orgdb: org.Mm.eg.db    # 마우스: org.Mm.eg.db / 인간: org.Hs.eg.db
```

**기술적 고려사항 (난이도 상승 요인):** Python–R 인터페이스(rpy2)가 환경에 따라 설치 까다로움. 비교 결과 CSV를 R script에 넘겨 filtering 후 돌려받는 파이프라인 설계가 더 안정적이다.

**산출물:** `kegg_flip_*.csv`, `bp_common_up_rrvgo.csv`, `final_dotplot_bp.svg`, `final_dotplot_kegg.svg`

---

### Phase 5 — 파이프라인 통합 ✦ 난이도: 높음

**목표:** Cross-condition 분석을 3단계 완료 후 자동으로 트리거되는 4단계로 파이프라인에 편입한다. 새 데이터셋에 `config.yaml`만 교체하면 재실행 가능하도록 한다.

**구현 항목:**
- 파이프라인 오케스트레이터(Snakemake 또는 Nextflow)에 `cross_condition` rule 추가. 모든 조건의 GO enrichment 완료를 dependency로 설정.
- `config.yaml` 검증 로직: conditions, groups, file_template, fig_prefix_map 일관성 자동 체크.
- species 파라미터 추가: orgdb 자동 선택(마우스/인간/기타).
- Phase 1–4 출력물을 종합한 **Excel 멀티시트 리포트** 자동 생성.
- 문서화: 새 데이터셋 적용 방법, config 수정 가이드.

**기술적 고려사항:** 파이프라인 프레임워크 의존성이 생긴다. 없는 환경이면 `run_analysis.py --config config.yaml` 단독 실행으로도 작동하도록 독립성을 유지한다.

**산출물:** `{project_name}_go_comparison.xlsx`, `Snakefile (또는 .nf)`, `config_template.yaml`

---

### 난이도 요약 및 권장 순서

Phase 1·2는 순수 Python 집합 연산으로 의존성이 없어 빠르게 구축 가능하다. Phase 3는 matplotlib 숙련도에 따라 1–3일 소요. Phase 4의 rpy2 인터페이스가 환경 설정 시간을 가장 많이 잡아먹는 구간이다. R script 분리 방식이 더 안정적. Phase 5는 파이프라인 전체 구조를 이해해야 하므로 마지막에 진행하는 것이 맞다.

**권장 순서:** Phase 1·2 → Phase 3 (dot plot 먼저) → Phase 4 (KEGG 먼저, rrvgo 나중) → Phase 5

---

## 8. 프로젝트 적용 체크리스트

새 프로젝트에 이 분석을 적용할 때 확인할 항목:

**데이터 준비:**
- [ ] 모든 조건의 GO enrichment CSV가 완성되었는가?
- [ ] CSV 컬럼명이 표준 형식과 일치하는가? (다르면 `col_map`에서 매핑)
- [ ] 조건별 enriched term 수가 충분한가? (방향별 ≥ 20개 권장)

**config.yaml 설정:**
- [ ] `conditions` 목록이 분석할 모든 조건을 포함하는가?
- [ ] `groups`가 방향 역전 분석에 맞게 논리적으로 정의되었는가?
- [ ] `file_template`이 실제 파일 경로 패턴과 일치하는가?
- [ ] `fig_prefix_map`이 필요한 경우(폴더명 ≠ 조건명) 설정되었는가?
- [ ] `min_conditions_common`이 조건 수에 맞게 조정되었는가?
- [ ] `orgdb`가 해당 종(species)에 맞게 설정되었는가?

**분석 실행:**
- [ ] Phase 1 실행 후 조건별 term count 로그 확인: 희소 조건이 있는가?
- [ ] 희소 조건이 있다면 `min_conditions_common` 조정 또는 해당 조건 제외 검토
- [ ] Phase 2 실행 후 역전 term의 gene Jaccard overlap 확인: 낮은 overlap은 별도 해석
- [ ] Semantic filtering 결과가 논문에 쓰기 적절한 수(10–30개)로 줄었는가?

**결과 해석:**
- [ ] 기존 time series / co-expression module GO 결과와 중복 여부 비교했는가?
- [ ] 역전 term 중 Jaccard < 0.3인 항목은 해석 시 주의 표기했는가?
- [ ] 온톨로지별(BP, KEGG) 결과가 서로 일관된 생물학적 내러티브를 지지하는가?
