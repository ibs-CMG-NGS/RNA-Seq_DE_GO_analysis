# atac-seq-da-analysis 구조 정합 요청서

**작성 배경**: RNA-Seq_DE_GO_analysis에 cross-dataset 비교 도구
(`run_cross_dataset_go_comparison.R`, `run_cluster_cross_dataset_comparison.R`)
를 만들면서 atac-seq-da-analysis의 산출물을 직접 대조 확인했다. GO enrichment
결과 자체(clusterProfiler `enrichGO`/`enrichKEGG` 표)는 두 파이프라인이 이미
동일한 컬럼 스키마(`ID`/`Description`/`GeneRatio`/`BgRatio`/`p.adjust`/`geneID`/
`Count`)로 수렴하고, GO ID가 species/assay-agnostic 온톨로지라 비교 로직도
완전히 같다. 다만 **그 표를 어디서 찾을지, geneID가 어떤 형식인지**의 3가지
지점에서 두 레포가 갈려서, 지금은 `assay` 프리셋이라는 어댑터 코드로 차이를
흡수하고 있다(RNA-Seq_DE_GO_analysis, `src/analysis/18_....R`/`19_....R`의
`ASSAY_PRESETS`).

**요청 취지**: 이 어댑터는 실용적인 단기 해법이지 최종 목표가 아니다. 아래 변경이
atac-seq-da-analysis 쪽에 반영되면(전부 또는 일부라도) 어댑터의 `atac` 프리셋을
지워도 되고, 장기적으로 두 파이프라인을 정말 "하나의 cross-dataset 관리 체계"로
운용할 수 있다. 실제 코드 변경은 이 문서의 범위가 아니며, 반영 여부/시점은
atac-seq-da-analysis를 다루는 쪽(사용자 본인 또는 별도 세션)이 판단한다.

---

## 1. 확인된 차이 (실측)

| 항목 | RNA-Seq_DE_GO_analysis | atac-seq-da-analysis |
|---|---|---|
| pairwise GO/KEGG 경로 | `pairwise/{pair}/enrichment/go_enrichment_{dir}_{ont}.csv`, `.../enrichment/kegg_enrichment_{dir}.csv` | `pairwise/{pair}/go_enrichment_{dir}_{ont}.csv`, `.../kegg_enrichment_{dir}.csv` (서브폴더 없음) |
| geneID 컬럼 형식 | Entrez ID (`enrichGO(..., readable=FALSE)`, `03_enrichment_analysis.R`) | Gene SYMBOL (`enrichGO(..., readable=TRUE)`, 여러 스크립트에 명시) |
| time-series 클러스터 GO | `time_series[_{variant}]/go_termcluster_cluster{N}_BP.csv`(variant=폴더명 접미사, 서브폴더 없음, variant 없는 base 케이스도 있음) | `time_series/{variant}/go_enrichment/masigpro_cluster_{k}_BP.csv`(variant=항상 중첩 서브폴더, base 케이스 없음) |
| coexpression 모듈 GO | `coexpression_modules[_{variant}]/go_termcluster_module{N}_BP.csv` | `coexpression_modules/[{variant}/]go_enrichment/module_{id}_BP.csv`(base도 자체 go_enrichment/ 보유) |
| 나머지 컬럼 | 동일 | 동일 |

## 2. 요청 사항 (비용순)

### 2a. [저비용] pairwise GO/KEGG를 `enrichment/` 서브폴더로

RNA-Seq_DE_GO_analysis는 이미 이 구조 전환을 6개 프로젝트에 실제로 적용해봤다
(단순 `mv`, 파일명 접두사 기준 분류). 참고용 패턴:

```bash
for pair in <pairwise 폴더 목록>; do
  d="output/<project>/pairwise/$pair"
  mkdir -p "$d/enrichment" "$d/plots"
  mv $d/go_enrichment_*.csv $d/go_rrvgo_*.csv $d/go_slim_*.csv $d/go_termcluster_*.csv $d/kegg_enrichment_*.csv "$d/enrichment/"
  mv $d/go_barplot_*.png $d/go_dotplot_*.png $d/go_rrvgo_scatter_*.png $d/go_rrvgo_treemap_*.png $d/go_slim_overview_*.png $d/go_termcluster_*.png $d/kegg_dotplot_*.png "$d/plots/"
done
```

기존 프로젝트는 그대로 두고(§3 참고), 이 구조를 생성하는 스크립트(pairwise GO/KEGG
enrichment를 도는 단계)의 출력 경로만 신규 프로젝트부터 `enrichment/`를 거치도록
바꾸면 된다.

### 2b. [저비용] time_series/coexpression_modules variant 구조를 폴더명 접미사로

`time_series/{variant}/go_enrichment/...` → `time_series_{variant}/...`,
`coexpression_modules/{variant}/go_enrichment/...` → `coexpression_modules_{variant}/...`
로 평탄화하고 `go_enrichment/` 중첩 래퍼를 제거. ATAC 쪽 masigpro/degPatterns
클러스터링 스크립트가 출력 경로를 만드는 지점만 수정하면 되는, 로직 변경이 아닌
경로 변경.

### 2c. [고비용, 즉시 강제하지 않음] geneID 포맷(`readable=`) 통일

RNA는 Entrez(`readable=FALSE`), ATAC은 SYMBOL(`readable=TRUE`)이다. 어느 쪽으로
통일하든 **이미 완료된 모든 프로젝트의 enrichGO/enrichKEGG를 다시 돌려야
소급 적용된다** — 저비용 항목들과 달리 새 프로젝트만 바꿔서 될 일이 아니다(과거
결과와 신규 결과의 geneID 포맷이 프로젝트마다 달라지면 오히려 더 혼란스러움).
지금은 어댑터가 이 차이를 흡수하고 있으니 당장 강제하지 않는 것을 권장한다.
파이프라인을 어떤 이유로든 대규모 재실행할 계기(예: reference genome 업데이트)가
생기면 그때 같이 논의.

참고로 통일한다면 SYMBOL(`readable=TRUE`, ATAC 쪽 현재 방식) 쪽을 표준으로
삼는 것을 제안한다 — CSV를 사람이 직접 열어봤을 때 바로 읽히고, cross-dataset
비교에서 종간 대소문자 정규화(`toupper()`) 한 단계만 거치면 되어 Entrez처럼 별도
매핑 테이블 구축이 필요 없다.

## 3. 적용 범위 권고

**기존에 완료된 프로젝트는 레거시로 그대로 두고, 신규 프로젝트부터만 새 구조를
적용한다.** RNA-Seq_DE_GO_analysis 자신도 이 원칙으로 6개 프로젝트만 선별
백필하고 나머지는 손대지 않았다 — 이미 나간 산출물/다운스트림 분석을 망가뜨릴
이유가 없고, cross-dataset 어댑터가 있는 한 구버전 구조도 계속 지원 가능하다
(RNA-Seq_DE_GO_analysis의 `2026-06-hiy-human-rna`가 구버전 flat 구조로 남아있다가
이번에 뒤늦게 발견되어 수동 마이그레이션한 사례 참고 — 정합이 완전히 끝나기
전까지는 이런 누락이 또 생길 수 있음을 감안할 것).

## 4. 기대 효과

2a/2b가 반영된 신규 프로젝트에 한해서는 `run_cross_dataset_go_comparison.R`/
`run_cluster_cross_dataset_comparison.R`의 `ASSAY_PRESETS`에서 `atac` 프리셋이
`rna` 프리셋과 동일해지므로, 어댑터 분기 자체를 지울 수 있다. 2c까지 반영되면
Entrez/SYMBOL 변환 로직도 제거 가능 — 결과적으로 "cross-dataset 비교"가 정말
`assay` 구분 없이 config(어떤 프로젝트 config를 참조하느냐)만으로 관리되는 하나의
도구가 된다.
