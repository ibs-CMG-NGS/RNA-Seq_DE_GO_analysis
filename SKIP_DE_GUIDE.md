# 외부 DE 분석 결과로 Downstream 분석 실행 가이드

이미 수행한 DE 분석 결과가 있는 경우, 파이프라인의 DE 분석 단계를 건너뛰고 enrichment analysis와 visualization만 실행할 수 있습니다.

## 📋 목차

1. [개요](#개요)
2. [필수 파일 형식](#필수-파일-형식)
3. [방법 1: 자동 스크립트 사용 (권장)](#방법-1-자동-스크립트-사용-권장)
4. [방법 2: 수동 배치](#방법-2-수동-배치)
5. [Downstream 분석 실행](#downstream-분석-실행)
6. [문제 해결](#문제-해결)

---

## 개요

### Downstream 분석에 포함되는 단계

DE 분석 이후 파이프라인에서 실행되는 단계들:

1. **Enrichment Analysis**
   - GO Enrichment (BP, CC, MF)
   - KEGG Pathway Enrichment
   - Dotplot 생성

2. **GO Visualization**
   - Barplot (Top 10 terms)
   - Combined plots

3. **Summary Table 생성**
   - `final_go_results.xlsx` (논문용 통합 테이블)

4. **QC Plots** (선택사항)
   - Volcano plot
   - MA plot
   - P-value histogram 등

---

## 필수 파일 형식

외부 DE 결과 파일(`final_de_results.csv`)에 포함되어야 할 컬럼:

### 필수 컬럼

| 컬럼명 | 설명 | 예시 |
|--------|------|------|
| `symbol` | 유전자 심볼 | "BRCA1", "TP53", "Cdhr1" |
| `log2FoldChange` | Log2 fold change | 2.408, -1.210 |
| `padj` | Adjusted p-value | 1.19e-24, 0.0038 |

### 권장 컬럼 (QC plots 생성시 필요)

| 컬럼명 | 설명 |
|--------|------|
| `baseMean` | 평균 발현량 |
| `pvalue` | Raw p-value |
| `lfcSE` | Log fold change standard error |
| `stat` | Test statistic |

### 파일 형식 예시

```csv
,symbol,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
ENSMUSG00000021803,Cdhr1,607.45,2.408,0.215,11.20,3.99e-29,1.19e-24
ENSMUSG00000001870,Ltbp1,1296.47,0.817,0.126,6.50,7.79e-11,1.16e-06
ENSMUSG00000049796,Crh,93.02,2.732,0.545,5.01,5.43e-07,0.0038
```

**중요**: 
- 첫 번째 컬럼은 유전자 ID (ENSEMBL ID 등)
- CSV 형식 (쉼표로 구분)
- 헤더 행 필수

---

## 방법 1: 자동 스크립트 사용 (권장)

### 1단계: 스크립트 실행

```bash
python prepare_external_de_results.py \
    --input your_de_results.csv \
    --comparison H2O2_vs_Control \
    --output-dir output/H2O2_Neuron \
    --config config_H2O2_Neuron.yml
```

### 파라미터 설명

- `--input`: 외부 DE 결과 CSV 파일 경로
- `--comparison`: 비교쌍 이름 (예: `H2O2_vs_Control`, `GABA_vs_Control`)
  - 형식: `{compare_group}_vs_{base_group}`
- `--output-dir`: 출력 디렉토리 (config의 `output_dir`와 동일)
- `--config`: Config YAML 파일 (검증용, 선택사항)

### 2단계: 검증

스크립트가 자동으로:
- ✅ 필수 컬럼 확인
- ✅ Config에 비교쌍이 정의되어 있는지 확인
- ✅ 올바른 경로에 파일 배치
- ✅ 다음 명령어 제안

### 예시 출력

```
============================================================
DE 결과 파일 준비 시작
============================================================

1. 입력 파일 읽기: my_de_results.csv
   - 행 수: 57,012
   - 컬럼 수: 8

2. 데이터 검증
✓ 필수 컬럼 확인 완료: symbol, log2FoldChange, padj
✓ 권장 컬럼 포함: baseMean, pvalue, lfcSE, stat

3. Config 파일 검증: config_H2O2_Neuron.yml
   ✓ 비교쌍 'H2O2_vs_Control'이 config에 정의되어 있습니다

4. 출력 디렉토리 생성: output/H2O2_Neuron/pairwise/H2O2_vs_Control

5. 파일 저장: output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv
   ✓ 파일 저장 완료

============================================================
✅ 준비 완료!
============================================================
```

---

## 방법 2: 수동 배치

스크립트 없이 수동으로 파일을 배치할 수도 있습니다.

### 1단계: 디렉토리 구조 생성

```bash
# Windows PowerShell
mkdir -p output/H2O2_Neuron/pairwise/H2O2_vs_Control

# 또는 WSL
mkdir -p output/H2O2_Neuron/pairwise/H2O2_vs_Control
```

### 2단계: 파일 복사

```bash
# Windows PowerShell
cp your_de_results.csv output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv

# 또는 WSL
cp your_de_results.csv output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv
```

### 3단계: Config 파일 복사 (선택사항)

```bash
cp config_H2O2_Neuron.yml output/H2O2_Neuron/pairwise/H2O2_vs_Control/config_used.yml
```

---

## Downstream 분석 실행

### 옵션 1: 특정 비교쌍의 모든 downstream 분석

```bash
snakemake --cores 4 \
    output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_go_results.xlsx
```

이 명령은 다음을 모두 실행합니다:
- GO/KEGG enrichment
- Dotplot 생성
- Barplot 생성
- Summary table 생성

### 옵션 2: Enrichment 분석만

```bash
snakemake --cores 4 \
    output/H2O2_Neuron/pairwise/H2O2_vs_Control/.enrichment_done.flag
```

### 옵션 3: 특정 단계만 강제 재실행

```bash
# GO enrichment만 다시 실행
snakemake --cores 4 --forcerun go_enrichment

# KEGG enrichment만 다시 실행
snakemake --cores 4 --forcerun kegg_enrichment

# Barplot만 다시 실행
snakemake --cores 4 --forcerun go_barplots
```

### 옵션 4: 여러 비교쌍 동시 실행

Config에 정의된 모든 비교쌍 실행:

```bash
snakemake --cores 4 --use-conda
```

특정 비교쌍들만:

```bash
snakemake --cores 4 \
    output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_go_results.xlsx \
    output/H2O2_Neuron/pairwise/GABA_vs_Control/final_go_results.xlsx
```

---

## 고급 사용법

### 1. Dry-run으로 실행 계획 확인

```bash
snakemake --dryrun --printshellcmds \
    output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_go_results.xlsx
```

### 2. DAG 시각화

```bash
snakemake --dag \
    output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_go_results.xlsx \
    | dot -Tpng > workflow_dag.png
```

### 3. 특정 규칙만 강제 재실행

```bash
# GO barplot과 summary table만 다시 생성
snakemake --cores 4 \
    --forcerun go_barplots generate_go_summary_table \
    output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_go_results.xlsx
```

### 4. Conda 환경 자동 생성

```bash
snakemake --cores 4 --use-conda \
    output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_go_results.xlsx
```

---

## 문제 해결

### ❌ "Input files not found" 오류

**원인**: DE 결과 파일이 올바른 경로에 없음

**해결**:
```bash
# 파일 존재 확인
ls output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv

# 없다면 다시 배치
python prepare_external_de_results.py --input your_de_results.csv ...
```

### ❌ "Column 'symbol' not found" 오류

**원인**: DE 결과 파일에 필수 컬럼이 없음

**해결**:
1. 파일 헤더 확인:
   ```bash
   head -1 your_de_results.csv
   ```

2. 컬럼명이 다른 경우 변경:
   ```python
   import pandas as pd
   df = pd.read_csv('your_de_results.csv')
   df.rename(columns={'gene_symbol': 'symbol', 'adj_pval': 'padj'}, inplace=True)
   df.to_csv('your_de_results_fixed.csv', index=False)
   ```

### ❌ Config에 비교쌍이 정의되지 않음

**원인**: `config.yml`의 `pairwise_comparisons`에 해당 비교쌍이 없음

**해결**:
```yaml
# config.yml에 추가
de_analysis:
  pairwise_comparisons:
    - ["H2O2", "Control"]
    - ["GABA", "Control"]
    # 새로운 비교쌍 추가
    - ["YourCompare", "YourBase"]
```

### ⚠️ Enrichment 결과가 비어있음

**원인**: 
- 유전자 심볼이 매칭되지 않음 (대소문자, species 불일치)
- padj 기준이 너무 엄격함

**해결**:
1. Species 확인:
   ```yaml
   # config.yml
   species: "mouse"  # 또는 "human"
   ```

2. padj cutoff 완화:
   ```yaml
   de_analysis:
     padj_cutoff: 0.1  # 기본값 0.05에서 완화
   ```

3. 유전자 심볼 확인:
   ```bash
   # DE 결과 파일의 심볼 몇 개 확인
   head your_de_results.csv
   # Mouse: Cdhr1, Ltbp1 (첫 글자만 대문자)
   # Human: BRCA1, TP53 (모두 대문자)
   ```

---

## 완전한 예시 워크플로우

### 시나리오: 3개의 비교쌍 분석

```bash
# 1. 준비: 각 비교쌍의 DE 결과 배치
python prepare_external_de_results.py \
    --input de_h2o2_vs_ctrl.csv \
    --comparison H2O2_vs_Control \
    --output-dir output/H2O2_Neuron \
    --config config_H2O2_Neuron.yml

python prepare_external_de_results.py \
    --input de_gaba_vs_ctrl.csv \
    --comparison GABA_vs_Control \
    --output-dir output/H2O2_Neuron \
    --config config_H2O2_Neuron.yml

python prepare_external_de_results.py \
    --input de_h2o2_vs_gaba.csv \
    --comparison H2O2_vs_GABA \
    --output-dir output/H2O2_Neuron \
    --config config_H2O2_Neuron.yml

# 2. Config 확인
cat config_H2O2_Neuron.yml

# 3. Dry-run으로 실행 계획 확인
snakemake --dryrun --printshellcmds

# 4. 전체 실행
snakemake --cores 4 --use-conda

# 5. 결과 확인
ls -lh output/H2O2_Neuron/pairwise/*/final_go_results.xlsx
```

---

## 생성되는 출력 파일

각 비교쌍(`{compare}_vs_{base}`)에 대해 다음 파일들이 생성됩니다:

### Enrichment 결과 (CSV)
```
output/{output_dir}/pairwise/{compare}_vs_{base}/
├── go_enrichment_up_BP.csv
├── go_enrichment_up_CC.csv
├── go_enrichment_up_MF.csv
├── go_enrichment_down_BP.csv
├── go_enrichment_down_CC.csv
├── go_enrichment_down_MF.csv
├── go_enrichment_total_BP.csv
├── go_enrichment_total_CC.csv
├── go_enrichment_total_MF.csv
├── kegg_enrichment_up.csv
├── kegg_enrichment_down.csv
└── kegg_enrichment_total.csv
```

### Visualization (PNG)
```
├── go_dotplot_up_BP.png
├── go_dotplot_up_CC.png
├── go_dotplot_up_MF.png
├── go_dotplot_down_BP.png
├── go_dotplot_down_CC.png
├── go_dotplot_down_MF.png
├── kegg_dotplot_up.png
├── kegg_dotplot_down.png
├── go_barplot_up_BP.png
├── go_barplot_up_CC.png
├── go_barplot_up_MF.png
├── go_barplot_down_BP.png
├── go_barplot_down_CC.png
└── go_barplot_down_MF.png
```

### Summary Table (Excel)
```
└── final_go_results.xlsx
    ├── Sheet: Up_BP
    ├── Sheet: Up_CC
    ├── Sheet: Up_MF
    ├── Sheet: Down_BP
    ├── Sheet: Down_CC
    └── Sheet: Down_MF
```

---

## 참고 자료

- [메인 README](README.md)
- [Snakemake 문서](https://snakemake.readthedocs.io/)
- [GO Enrichment 테이블 기능 설명](GO_TABLE_FEATURE_SUMMARY.md)
