# GSE 공공데이터 메타분석 — Cross-Dataset 매핑안

## 목적 및 배경

CNS(중추신경계) 손상 관련 공공 GEO RNA-seq 데이터셋 8개를 확보해 각 프로젝트별
pairwise DE / time-series / coexpression module 분석을 파이프라인화했다. 이 문서는
그 다음 단계 — 프로젝트를 넘나드는 **cross-dataset 메타분석**을 어떤 축으로,
어떤 비교 쌍으로 진행할지에 대한 매핑안이다. 실제 실행 스크립트
(`run_cross_dataset_go_comparison.R`, `run_cluster_cross_dataset_comparison.R`)의
사용법은 `docs/CROSS_DATASET_GUIDE.md`를 참고.

**중요**: 이 문서의 시점(timepoint)·조직·축 배정은 각 데이터셋의 메타데이터
파일과 config에서 기계적으로 추출한 것이며, 원 논문의 실험적 맥락(예: 손상
중증도, 마취 방법, 조직 채취 정밀 시점)까지 반영한 것은 아니다. **다른 연구자의
검토·수정을 전제로 한 초안**이다.

## 2. 사용 가능 티어 요약

아래 3~4절의 상세 근거를 먼저 한눈에 보기 위한 요약. 각 데이터셋이 이번 매핑안
(6절의 축 1~4 + 보조 비교)에서 어느 수준으로 쓰일 수 있는지 기준.

| 티어 | 데이터셋 | 근거 |
|---|---|---|
| **Tier 1** — 즉시 사용, 신뢰도 높음 | GSE104036, GSE205486 | replicate 정상(조건당 n=3 이상), 설계 문제 없음, time_series/coexpr 트랙 완료(8절) — 축 1~4의 중심축 |
| **Tier 2** — 사용 가능, 캐비어트 동반 | GSE142445 | 트랙 보유·축 1~3 참여 가능하나 조건당 n=1이라(4-2절) pairwise(run_cross_dataset_go_comparison.R)는 보조적으로만, time_series(run_cluster_cross_dataset_comparison.R) 우선 신뢰 |
| | GSE326470 | replicate 정상이나 원발 병소가 아닌 원격 조직 + rat — 축 4(탐색적)에만 |
| | GSE286075 | replicate 정상이나 bulk가 아닌 세포타입(성상세포) 특이 신호 — 보조 비교에만 |
| **Tier 3** — 이번 라운드 제외, 재작업 시 재검토 | GSE173544 | 실제 시점(7주)이 현재 축의 acute/chronic 기준과 맞지 않음(4-1절) — "매우 만성" 축이 생기면 재검토 |
| | GSE155610 | 조직 confounding + 애초 다른 질문(치료 개입 효과)의 데이터(4-3절) — 재설계·별도 축 신설 필요 |
| **Tier 4** — 완전 배제 | GSE234052 | replicate 없음, DESeq2 구동 자체가 불가 — 재시퀀싱 전까지 불가 |

## 3. 데이터셋 인벤토리

| 데이터셋 | 종 | 조직 | 손상 모델 | 시점(시간 환산) | time_series/coexpr 트랙 | 비고 |
|---|---|---|---|---|---|---|
| GSE104036 | mouse | 뇌(MCAO) | 허혈성 뇌졸중 | 3, 6, 12, 24h (ipsi/contra vs sham) | 있음 (ipsi+contra) | 급성기만, 만성 시점 없음 |
| GSE142445 | mouse | 뇌(MCAO) | 허혈성 뇌졸중 | 4, 24, 72, 168h (ipsi/contra vs uninjured) | 있음 (ipsi+contra) | 급성+만성 모두 포함, 각 시점×반구 **n=1**(uninjured만 n=3) — 4절 참고 |
| GSE205486 | mouse | 척수 | 척수손상(contusion) | 24, 48, 72, 120, 168h (Lesion vs Naive) | 있음 (Lesion 주/Naive 보조) | D1만 급성 경계, 나머지 만성 |
| GSE326470 | **rat** | 대뇌피질 (SCI의 **원격** 반응 — 원발 병소 아님) | 척수손상(SCI) | 168h(D7, Male만), 720h(D30, Male+Female) | 없음(단일 비교만) | 원발 부위 아닌 원격 조직 반응이라는 점이 특징 |
| GSE173544 | mouse | 뇌(distal MCAO+저산소) | 만성 허혈성 뇌졸중 | **7주(약 1176h)** — 확인됨(4-1절) | 없음 | 다른 데이터셋의 "chronic"(7d)보다 훨씬 늦은 시점 — 별도 취급 |
| GSE286075 | mouse | 뇌(MCAO), **성상세포 endfeet 특이 translatome(RiboTag)** | 허혈성 뇌졸중 | **6h 재관류(2h 허혈+6h reperfusion)** — 확인됨(4-1절) | 없음 | bulk가 아니라 세포타입 특이 신호 — 급성 초기(3~6h대) 매칭 가능 |
| GSE155610 | **rat** | 척수 + 피질 백질(5개 하위 조직) | 경추 편측 SCI **전 개체 이미 손상**, 자극 유무만 비교 | 1주(손상 후, 자극군은 그 1주간 자극) | 없음 | **비교 축 자체가 다름** — 4-1절 참고, "손상 vs 정상"이 아니라 "손상+자극 vs 손상만" |
| GSE234052 | mouse | 뇌(MCAO) | 허혈성 뇌졸중 | 1, 12, 24h | — | **제외 확정** — 조건당 replicate 1개뿐이라 DESeq2 분산 추정 자체가 불가능(`checkForExperimentalReplicates` 에러). 재시퀀싱 데이터가 생기기 전까지 이번 라운드에서 계속 배제 |

## 4. 표본수/설계 감사

각 데이터셋 metadata의 조건별 표본수를 전수 확인한 결과, 8개 중 6개는 문제없고
(GSE104036/GSE173544/GSE286075/GSE326470은 조건당 n=3, GSE205486은 Naive_D7만
n=2로 경미), 나머지 2개(GSE142445, GSE155610)는 서로 다른 성격의 실질적 문제가
있다. GSE234052는 이미 제외 확정(7절).

### 4-1. 시점 확인 결과 (원 논문/GEO 레코드 조회)

3절에서 "불명"으로 표시했던 3개 데이터셋의 실제 시점을 GEO 레코드(SOFT 텍스트) +
연관 논문에서 확인했다.

- **GSE173544** — *RNA Sequencing of WT Chronic Stroke Infarcts vs. Contralateral
  Cortices*. Distal MCAO + 저산소(hypoxia) 모델, 3개월령 마우스, **뇌졸중 유발 후
  7주(약 49일, 1176h) 시점**에 조직 채취. 다른 데이터셋들의 "chronic" 기준점(7일)보다
  7배 가까이 늦은 시점이라, 기존 축 2/3의 `Chronic_vs_Control`(7일 기준)에 그대로
  끼워 넣으면 시점 차이가 너무 커서 왜곡된 해석을 유발할 수 있다. → **축 2/3의
  Acute/Chronic 이분법에서 제외**하고, 추후 "매우 만성(수주~수개월)" 축이 필요할 때
  단독 참고점으로 남겨둔다. 연관 논문: [Shared transcriptomic signatures in
  perilesional and contralesional cortex after ischemic stroke](https://pubmed.ncbi.nlm.nih.gov/42298604/).
- **GSE286075** — *RiboTag RNA Sequencing Identifies Local Translation of HSP70 in
  Astrocyte Endfeet After Cerebral Ischemia*. FVB 배경 astrocyte-specific RiboTag
  형질전환 마우스, 10–12주령 수컷, **2h 허혈(MCAO) + 6h 재관류** 후 뇌 미세혈관
  (성상세포 endfeet) 채취, 대조/허혈 각 n=3(메타데이터와 일치). GSE104036의 3–6h대,
  GSE142445의 4h와 시기상 겹치는 **급성 초기** 시점 — 다만 bulk 조직이 아니라
  성상세포 endfeet 특이 translatome이라는 점은 여전히 유효한 차이. 연관 논문:
  [PMC11720067](https://pmc.ncbi.nlm.nih.gov/articles/PMC11720067/),
  [PubMed 39416227](https://pubmed.ncbi.nlm.nih.gov/39416227/).
- **GSE155610** — *Transcriptome of Subcortical White Matter and Spinal Cord After
  Spinal Injury and Cortical Stimulation*. Long-Evans 암컷 랫 4마리 **전부**
  경추(C4) 편측 contusion 손상을 받고, 그중 2마리만 손상 후 1주간 운동피질
  전기자극(역치 80%)을 추가로 받음 — 나머지 2마리(Control)는 **자극 없이 손상만**.
  즉 `Control` 조건도 이미 손상된 개체다. 시점은 손상 후 1주(약 168h). 조직은
  양측 피질 백질(subcortical white matter ipsi/contra) + 척수 3구간(rostral/
  lesion epicenter/caudal), 총 5개. 연관 논문(Scientific Data,
  [s41597-021-00953-4](https://www.nature.com/articles/s41597-021-00953-4)).
  **이 데이터셋을 축 2/3에 넣으면 안 되는 더 근본적인 이유가 드러났다** — 4-3절 참고.

### 4-2. GSE142445 — 조건당 n=1 (통계적 검정력 문제)

`design_formula: ~ condition`, 조건 9개(uninjured n=3 + 8개 timepoint×side 각
n=1) = 샘플 11개, 모델 매트릭스 계수 9개 → **잔차 자유도(residual df) = 2**.
GSE234052(잔차 df=0)처럼 크래시는 안 나지만, 유전자별 분산 추정이 사실상 전체
유전자의 dispersion trend curve(shrinkage)에 거의 전적으로 의존한다는 뜻 —
개별 유전자의 padj를 액면 그대로 신뢰하기 어렵다.

- **pairwise DE**(`run_cross_dataset_go_comparison.R`이 쓰는 단위)가 가장
  취약함 — n=1 vs n=3 비교는 그 개체 하나의 우연한 변동에 좌우되기 쉬움.
- **time_series**(`run_cluster_cross_dataset_comparison.R`이 쓰는 단위)는
  상대적으로 낫다 — maSigPro 다항회귀가 9개 조건 전체의 시간축 정보를 함께
  써서 추세를 추정하므로 개별 시점 n=1이 부분적으로 보완됨.
- **권장**: GSE142445가 관여하는 모든 pairwise(run_cross_dataset_go_comparison.R) 기반 비교(축 1/2/3)는
  GSE104036·GSE205486처럼 진짜 replicate가 있는 데이터셋과 같은 방향으로 나오는
  term만 우선 신뢰하고, 가능하면 run_cluster_cross_dataset_comparison.R(time_series 클러스터) 결과로 교차검증할 것.
  padj 단독보다 log2FC 크기 + pairwise/time_series 양쪽 일관성을 함께 볼 것.

### 4-3. GSE155610 — 조직 confounding + 비교 축 자체가 다름 (2중 문제)

**문제 A (기존 발견, 여전히 유효) — 조직 하위구분을 무시한 설계**: metadata에
`tissue` 컬럼이 있고 실제로는 5개 하위 조직(`subcortical_white_matter_
ipsi/contra`, `spinal_cord_rostral/caudal/lesion_epicenter`) × Control/Treatment ×
n=2로 균형 잡힌 설계인데, 현재 config는 `design_formula: ~ condition`,
`group_variable: condition`으로 **tissue를 완전히 무시**하고 5개 조직을 그대로
합쳐 Control n=10 / Treatment n=10으로 취급하고 있다. 조직 간 발현 baseline
차이가 훨씬 클 가능성이 높은 조직들(백질 vs 척수 여러 분절)을 섞은 것이라,
현재 DE 결과는 자극 효과가 아니라 조직 구성비 차이에 의해 지배되고 있을 위험이
크다.

- **원인**: 표본수 부족이 아니라 **설계 자체의 confounding** — n=2×5=10은
  조직별로 보면 오히려 충분한 편.
  - **가능한 수정 방향 두 가지(택 1)**:
    1. `design_formula: ~ tissue + condition`으로 바꿔 조직을 공변량으로
       통제(균형 설계라 바로 적용 가능).
    2. 조직별로 config를 쪼개 5개의 독립 프로젝트(또는 파생 config)로 재실행
       — 조직마다 자극 반응이 다를 수 있다는 가설 자체를 검증할 수 있다는 장점.

**문제 B (원 논문 확인 후 새로 발견, 더 근본적) — `Control`도 이미 손상된
개체다**: 4-1절에서 확인했듯 이 연구는 랫 4마리 **전부** 경추 편측 손상을
받았고, `Treatment`는 "손상 + 1주간 운동피질 자극", `Control`은 "손상만,
자극 없음"이다. 즉 우리 축 2/3가 전제하는 "손상 vs 정상(sham/naive)" 비교가
애초에 아니라 **"치료적 개입(자극) 효과"** 비교다. 문제 A를 고쳐 조직
confounding을 없애더라도, 이 데이터셋을 다른 데이터셋들의 `Acute_vs_Control`/
`Chronic_vs_Control`과 같은 칸에 넣으면 안 된다 — 비교하는 대상 자체가 다르기
때문이다(손상 반응 vs 개입 효과는 다른 질문).

- **결론**: GSE155610은 문제 A를 고친 뒤에도 축 2/3에는 참여시키지 않는다.
  대신 "손상 후 치료적 개입이 손상 시그니처를 되돌리는가"라는 **별도의 질문**을
  던질 수 있는 유일한 데이터셋이므로, 향후 그 질문 전용의 새 축(예: 축 5 —
  개입/치료 효과, 이번 매핑안 범위 밖)으로 남겨둔다.

## 5. 매핑 메커니즘 (참고)

`cross_dataset_*.yaml`의 최상위 `pairs:`가 "정식(canonical) 비교 이름"이고,
각 `datasets[].pair_map`이 그 정식 이름을 그 프로젝트의 실제 pairwise 폴더명에
매핑한다. 즉 `pairs: ["Acute_vs_Control", "Chronic_vs_Control"]`처럼 데이터셋마다
실제 시점이 달라도 추상화된 라벨로 묶을 수 있다 — 아래 축들은 전부 이 메커니즘
그대로 표현 가능(코드 변경 불필요).

- `run_cross_dataset_go_comparison.R`: pairwise GO/KEGG term 비교(common/flip/
  exclusive/mixed). `pair_map` 기반 — time_series/coexpression 트랙이 없는
  데이터셋도 참여 가능.
- `run_cluster_cross_dataset_comparison.R`: time_series/coexpression 클러스터의
  GO term-set Jaccard 매칭. 데이터 기반 자동 매칭이라 `pair_map` 불필요하지만,
  **트랙이 있는 데이터셋끼리만** 비교 가능(GSE104036/142445/205486만 해당).

## 6. 제안 축

### 축 1 — 뇌졸중 급성기, 뇌 조직, 마우스 한정 (1차 검증용, 신뢰도 최고)

같은 손상모델·조직·종에 시점도 겹쳐 매칭 근거가 가장 탄탄하다.
`run_cross_dataset_go_comparison.R`, `run_cluster_cross_dataset_comparison.R`
둘 다 가능 — 파이프라인/매핑 로직 자체의 정합성을 먼저 이걸로 검증하는 것을
권장. **단, GSE142445는 조건당 n=1이라(4절) pairwise(run_cross_dataset_go_comparison.R)
결과는 보조적으로만 보고 time_series 클러스터(run_cluster_cross_dataset_comparison.R)
결과를 우선 신뢰할 것.**

| 정식 이름 | GSE104036 (ipsi) | GSE142445 (ipsi) |
|---|---|---|
| `Acute_Early_vs_Baseline` | `3hr_MCAO_ipsilateral_vs_sham` | `4hr_ipsilateral_vs_uninjured` |
| `Acute_Late_vs_Baseline` | `24hr_MCAO_ipsilateral_vs_sham` | `1d_ipsilateral_vs_uninjured` |

GSE104036의 6hr/12hr 시점은 GSE142445에 대응 시점이 없어 이 cross-dataset
축에서는 제외(각 프로젝트 자체 분석에서는 그대로 유지). contra 트랙도 동일한
패턴으로 별도 config(`Acute_Early_Contra_vs_Baseline` 등)로 추가 가능.

### 축 2 — Acute(≤24h) vs Chronic(>24h), 조직·종 초월

24h를 경계로 잡고, 각 데이터셋에서 그 경계에 가장 가까운 시점을 선택(급성은
"경계에 가장 가까운 늦은 시점", 만성은 "가장 늦은 시점"을 기본 규칙으로 채택 —
반응 크기가 시점 끝에서 가장 뚜렷하다는 가정). run_cross_dataset_go_comparison.R만 가능(run_cluster_cross_dataset_comparison.R은 데이터셋별
자체 시간축을 이미 쓰므로 이 축과는 별개 실행).

| 정식 이름 | GSE104036 | GSE142445 | GSE205486 | GSE326470 |
|---|---|---|---|---|
| `Acute_vs_Control` | `24hr_MCAO_ipsilateral_vs_sham` | `1d_ipsilateral_vs_uninjured` | `Lesion_D1_vs_Naive_D1` | (없음 — acute 데이터 없어 제외) |
| `Chronic_vs_Control` | (없음 — chronic 데이터 없어 제외) | `7d_ipsilateral_vs_uninjured` | `Lesion_D7_vs_Naive_D7` | `SCI_Day30_Male_vs_Sham_Day30_Male` |

**주의**: `Chronic_vs_Control`의 GSE326470은 30일 시점으로 다른 데이터셋의 7일과
꽤 차이가 있고, 종도 rat이라 유전자 수준 비교는 symbol 대소문자 일치 수준(오소로그
매핑 없음)이라는 한계가 있음 — GO term 수준 비교로 해석을 제한하는 게 안전.
`Acute_vs_Control`의 GSE142445(n=1, 4절)도 동일하게 보조적으로만 취급.

### 축 3 — 조직 공통/특이 (뇌 vs 척수), 급성기 고정

축 2의 `Acute_vs_Control` 매핑을 그대로 재사용하되, 데이터셋을 뇌(GSE104036,
GSE142445)와 척수(GSE205486)로 나눠서 비교. 뇌 2개에서 공통으로 바뀌면서
척수에서는 안 바뀌는(또는 반대 방향인) term = 조직 특이 신호, 뇌·척수 모두에서
같은 방향 = 손상 자체에 공통적인 신호로 해석.

GSE173544는 4-1절에서 확인된 실제 시점이 7주(급성이 아니라 훨씬 늦은 만성)라
이 축에서 **제외**한다 — 애초 "불명"일 때 잠정 배치했던 것을 되돌림.

### 축 4 (탐색적) — 원발 병소 vs 원격 반응

같은 "Day7 척수손상"을 원발 부위(GSE205486, 척수 자체)와 원격 부위(GSE326470,
대뇌피질)에서 각각 본다는 점에서 지난 디커세이션 논의와 직접 연결되는 축.
데이터셋이 각 방향에 1개씩뿐이라 통계적 근거는 약하고 **가설 생성용**으로만
사용 권장.

| 정식 이름 | GSE205486 (원발, 척수) | GSE326470 (원격, 대뇌피질) |
|---|---|---|
| `Day7_Injury_vs_Control` | `Lesion_D7_vs_Naive_D7` | `SCI_Day7_Male_vs_Sham_Day7_Male` |

### 보조 비교 — Bulk vs 세포타입 특이(성상세포)

같은 MCAO 뇌졸중이라는 전제로 GSE104036/GSE142445(bulk) vs GSE286075(성상세포
endfeet 특이 RiboTag translatome)를 비교해 bulk 시그널이 성상세포 기원인지
확인하는 용도. 4-1절에서 확인된 GSE286075의 시점(2h 허혈+6h 재관류)은
GSE104036의 3h·6h, GSE142445의 4h와 겹치는 급성 초기라 — 축 1의
`Acute_Early_vs_Baseline`과 **같은 시점 대역**으로 짝지어 비교 가능(단, bulk vs
성상세포 endfeet이라는 조직 수준 차이는 해석에 계속 반영해야 함).

| 정식 이름 | GSE104036 (ipsi, bulk) | GSE142445 (ipsi, bulk) | GSE286075 (성상세포 endfeet) |
|---|---|---|---|
| `Acute_Early_vs_Baseline` | `3hr_MCAO_ipsilateral_vs_sham` | `4hr_ipsilateral_vs_uninjured` | `ipsi_vs_contra` |

## 7. 제외/보류 항목

- **GSE234052**: replicate 없음 — 이번 라운드 전체 제외(위 표 참고).
- **GSE155610**: 조직 confounding(4-3절 문제 A) + 비교 축 자체가 다름(4-3절
  문제 B, "손상 vs 정상"이 아니라 "손상+개입 vs 손상만") — 재설계 및 별도 축
  신설 전까지 이번 매핑안에서 계속 제외.
- **GSE173544**: 시점 확인 결과(4-1절) 7주(급성 아님, 다른 데이터셋의 chronic
  기준 7일보다도 훨씬 늦음) — 현재 축 2/3의 acute/chronic 이분법과 맞지 않아
  이번 라운드는 제외, 향후 "매우 만성" 축이 생기면 단독 참고점으로 재검토.
- **GSE286075**: 시점 확인 완료(4-1절, 2h 허혈+6h 재관류) — "보조 비교"에
  `Acute_Early_vs_Baseline`으로 확정 배정(더 이상 잠정 아님).
- **GSE142445**: 제외는 아니지만 4-2절의 n=1 캐비어트를 전제로만 사용(축 1/2/3에
  전부 해당).

## 8. 현재 진행 상태

- GSE142445, GSE104036: time_series(ipsi+contra 2트랙 동등 실행)/coexpression
  모듈 분석 포함 전체 파이프라인 완료(GDrive 업로드까지).
- GSE326470: time_series는 maSigPro의 시점(time point) 최소 요구(3개 이상)를
  Male 스트라텀도 충족 못 해(Day7/Day30, 2개뿐) 비활성화 확정 — coexpression만
  완료.
- **GSE205486 — 정정 + 완료**: 앞서 이 절에서 "완료"라고 기록했었는데 실제로는
  config에 `time_series`/`coexpression_modules`가 설정만 되어 있고 한 번도
  실행된 적이 없었다(축 2/축 4 cross-dataset cluster 비교가 "0개 유의 클러스터"로
  조용히 건너뛰어지는 걸 보고 뒤늦게 발견). 실행 중 두 번째 문제도 발견 —
  `coexpression_modules`가 dict(단일 트랙)인데 `variant_label: lesion`이 남아있어
  R 스크립트가 쓰는 staging 파일명과 Snakefile이 기대하는 파일명이 어긋나
  Snakemake가 "output 미생성"으로 판단해 이미 완성된 CSV 등을 통째로 삭제하던
  버그였다 — `variant_label` 제거 후 재실행해 **완전히 완료**(GDrive 업로드까지).
- 축 1~4 + 보조 비교의 cross-dataset config 5개(`configs/cross_dataset_gse-*.yaml`)
  전부 **GO 비교(pairwise)와 cluster 비교(time_series/coexpr Jaccard) 둘 다
  실행 완료**. 주목할 만한 신규 발견: 축 3(뇌 vs 척수)에서 GSE142445(뇌졸중,
  뇌)와 GSE205486(척수손상) 사이 coexpression 모듈 한 쌍이 Jaccard 0.696(공유
  39/합집합 56 term)으로 겹침 — 서로 다른 조직의 손상 반응인데도 거의 동일한
  발현 패턴 모듈이 존재한다는 뜻. `run_cross_dataset_go_comparison.R`/
  `run_cluster_cross_dataset_comparison.R` 산출물은 각각
  `output/cross_dataset_gse-{acute-brain,acute-vs-chronic,brain-vs-spinalcord,
  primary-vs-remote,bulk-vs-astrocyte}/`에 있다.
- **1단계(term-level p-value 결합, Fisher's method/Stouffer's Z) 구현 완료** —
  `run_cross_dataset_go_comparison.R`에 `meta_analysis` 옵션으로 추가됨(순수
  additive, 기존 common/flip/exclusive/mixed 결과는 그대로 유지). GSE142445가
  관여하는 축 1~3은 `method: stouffer` + `weights: {GSE142445: 0.5}`로 n=1
  캐비어트(4-2절)를 실제 통계 결합에 반영해뒀다. 사용법은
  `docs/CROSS_DATASET_GUIDE.md`의 `meta_analysis` 절 참고.
- 2단계(유전자 log2FC 레벨 메타분석, 오소로그 매핑 선행 필요)는 아직 범위 밖.
