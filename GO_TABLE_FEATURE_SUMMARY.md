# GO Enrichment Summary Table 기능 추가 완료 ✨

## 📋 개요

논문 supplementary material 형식의 GO enrichment 결과 통합 Excel 파일을 자동 생성하는 기능이 추가되었습니다.

## 🎯 추가된 기능

### 1. 새로운 스크립트: `05_generate_go_table.R`

**위치**: `src/analysis/05_generate_go_table.R`

**기능**:
- 모든 GO enrichment CSV 결과를 하나의 Excel 파일로 통합
- 다중 워크시트 구성으로 체계적인 정리
- 논문 제출용 전문 서식 자동 적용

**출력 파일**: `GO_enrichment_summary.xlsx`

### 2. Excel 파일 구조

파일에는 다음 9개의 워크시트가 포함됩니다:

| 시트 이름 | 내용 |
|-----------|------|
| **All_Results** | 모든 GO enrichment 결과 통합 (메인 시트) |
| **UP_regulated** | 상향 조절 유전자의 GO terms |
| **DOWN_regulated** | 하향 조절 유전자의 GO terms |
| **TOTAL_regulated** | 전체 유의 유전자의 GO terms |
| **Biological_Process** | BP ontology의 모든 결과 |
| **Cellular_Component** | CC ontology의 모든 결과 |
| **Molecular_Function** | MF ontology의 모든 결과 |
| **Top_Terms** | 각 카테고리의 상위 20개 GO term 요약 |
| **Analysis_Info** | 분석 파라미터 및 메타데이터 |

### 3. 컬럼 구성

각 시트는 다음 컬럼을 포함합니다:

- **Gene Set**: UP, DOWN, TOTAL
- **Ontology**: BP, CC, MF
- **GO ID**: GO term 식별자 (예: GO:0030072)
- **GO Term**: GO term 설명
- **Gene Ratio**: 내 유전자 중 해당 GO에 포함된 비율
- **Background Ratio**: 전체 유전자 중 해당 GO에 포함된 비율
- **P-value**: 통계적 유의성 (scientific notation)
- **Adjusted P-value**: FDR 보정된 p-value
- **Q-value**: 보정된 q-value
- **Gene Count**: 해당 GO에 포함된 유전자 개수
- **Gene IDs**: 해당하는 유전자 ID 목록 (세미콜론으로 구분)

### 4. 서식 및 스타일

- ✅ **헤더**: 파란색 배경, 흰색 글자, 굵게
- ✅ **테이블**: 테두리 적용, 교대 색상 없음 (깔끔한 디자인)
- ✅ **P-value**: Scientific notation (예: 1.80E-06)
- ✅ **자동 열 너비**: 내용에 맞게 자동 조정
- ✅ **Freeze pane**: 헤더 행 고정 (스크롤 시 헤더 유지)
- ✅ **GO Term**: 50자 열 너비 (긴 설명도 잘 보임)
- ✅ **Gene IDs**: 60자 열 너비 (여러 유전자 ID 표시)

## 📦 Pipeline 통합

### Snakefile 수정사항

**새로운 규칙 추가**: `generate_go_summary_table`

```python
# Rule 5: Generate GO Summary Table for Publication
rule generate_go_summary_table:
    input:
        script = "src/analysis/05_generate_go_table.R",
        config_file = "config_H2O2_Neuron.yml",
        enrichment_flag = OUTPUT_DIR / "pairwise/{pair}/.enrichment_done.flag",
        go_csvs = lambda wildcards: expand(...)
    output:
        excel = OUTPUT_DIR / "pairwise/{pair}/GO_enrichment_summary.xlsx"
    ...
```

**target rule 업데이트**: `rule all`에 새 출력 파일 추가

```python
expand(OUTPUT_DIR / "pairwise/{pair}/GO_enrichment_summary.xlsx", pair=PAIRS)
```

### 실행 방법

#### 기존 프로젝트에서 GO table만 다시 생성:

```bash
# Snakemake 환경 활성화
conda activate snakemake_env

# GO table만 강제로 재생성
snakemake --cores 4 --use-conda --forcerun generate_go_summary_table

# 또는 특정 비교만
snakemake --cores 4 --use-conda --force output/H2O2_Neuron/pairwise/H2O2_vs_Control/GO_enrichment_summary.xlsx
```

#### 전체 파이프라인 실행:

```bash
# 기존과 동일하게 실행하면 자동으로 포함됨
snakemake --cores 4 --use-conda
```

## 🔧 Snakefile 업데이트 방법

`snakemake_rule_go_table.txt` 파일에 새 규칙이 작성되어 있습니다. 

**수동 추가 방법**:
1. `Snakefile` 파일을 열기
2. 파일 끝에 `snakemake_rule_go_table.txt`의 내용을 복사하여 붙여넣기
3. 저장

**또는 명령어로 추가**:
```bash
cd \\wsl.localhost\Ubuntu\home\ygkim\ngs_pipeline\RNA-Seq_DE_GO_analysis
cat snakemake_rule_go_table.txt >> Snakefile
```

## 📁 결과 파일 위치

```
output/
└── H2O2_Neuron/
    └── pairwise/
        ├── H2O2_vs_Control/
        │   ├── final_de_results.csv
        │   ├── go_enrichment_*.csv      # 개별 CSV 파일들
        │   ├── GO_enrichment_summary.xlsx  # ✨ 새로 생성되는 통합 Excel
        │   └── ...
        └── GABA_vs_Control/
            ├── ...
            └── GO_enrichment_summary.xlsx  # ✨ 각 비교마다 생성
```

## 💡 사용 예시

### 연구자에게 결과 전달

```
안녕하세요,

H2O2 vs Control 비교 분석 결과를 첨부합니다.

- final_de_results.xlsx: 차등 발현 유전자 목록
- GO_enrichment_summary.xlsx: GO enrichment 분석 결과 (모든 ontology 포함)

GO 파일의 'Top_Terms' 시트를 먼저 확인하시면 주요 결과를 빠르게 파악하실 수 있습니다.
각 ontology별로 상세한 결과는 별도 시트에 정리되어 있습니다.
```

### 논문 Supplementary Material

```
Supplementary Table S3. 
Gene Ontology (GO) enrichment analysis of differentially expressed genes 
in H2O2-treated neurons compared to control.

파일: GO_enrichment_summary.xlsx
- Sheet 1: All results
- Sheet 2-4: Gene set specific results (UP/DOWN/TOTAL)
- Sheet 5-7: Ontology specific results (BP/CC/MF)
- Sheet 8: Top 20 most significant terms
- Sheet 9: Analysis parameters
```

## ✅ 장점

1. **연구자 친화적**: Excel 형식으로 누구나 쉽게 열람 가능
2. **체계적 구성**: 워크시트별로 정리되어 원하는 정보 빠른 검색
3. **논문 준비 완료**: 바로 supplementary material로 제출 가능
4. **재현성**: 분석 파라미터가 파일 내 자동 기록됨
5. **전문성**: 깔끔한 서식으로 전문적인 인상
6. **시간 절약**: 수동으로 정리할 필요 없음

## 📚 문서 업데이트

README.md 파일이 다음 내용으로 업데이트되었습니다:

1. **스크립트 설명 섹션**: `05_generate_go_table.R` 추가
2. **폴더 구조**: 새 스크립트 파일 포함
3. **Snakemake 사용 가이드**: GO summary table 자동 생성 설명 추가

## 🎓 기술 상세

### 사용 패키지
- `openxlsx`: Excel 파일 생성 및 고급 서식 적용
- `dplyr`: 데이터 조작 및 필터링
- `yaml`: Config 파일 읽기

### 서식 적용 코드 예시
```r
header_style <- createStyle(
  fontSize = 11,
  fontName = "Arial",
  textDecoration = "bold",
  halign = "center",
  fgFill = "#4472C4",
  fontColour = "#FFFFFF",
  border = "TopBottomLeftRight"
)
```

### 데이터 변환
- GeneRatio, BgRatio: 문자열로 유지 (예: "4/9")
- P-values: Scientific notation (0.00E+00 형식)
- Count: 정수형
- Gene IDs: 세미콜론으로 연결된 문자열

## 🚀 향후 개선 가능사항

1. **KEGG 결과 통합**: GO와 마찬가지로 KEGG enrichment 결과도 별도 Excel로 생성
2. **시각화 포함**: Excel 내에 차트 자동 삽입
3. **필터링 옵션**: 사용자가 top N개 term 개수 조정 가능
4. **하이퍼링크**: GO ID를 클릭하면 GO 웹사이트로 이동
5. **조건부 서식**: P-value 크기에 따라 색상 변화

## ❓ 문제 해결

### Excel 파일이 생성되지 않는 경우

1. **GO 결과 파일 확인**:
   ```bash
   ls output/H2O2_Neuron/pairwise/H2O2_vs_Control/go_enrichment_*.csv
   ```

2. **로그 파일 확인**:
   ```bash
   cat output/H2O2_Neuron/pairwise/H2O2_vs_Control/logs/05_generate_go_table.log
   ```

3. **openxlsx 패키지 확인**:
   ```r
   # R에서
   library(openxlsx)
   ```

### Snakemake 규칙이 실행되지 않는 경우

```bash
# Dry-run으로 확인
snakemake -np --use-conda

# 강제 실행
snakemake --cores 4 --use-conda --forcerun generate_go_summary_table
```

---

**작성일**: 2025-12-01  
**작성자**: GitHub Copilot  
**관련 이슈**: GO analysis 결과 통합 및 배포 개선
