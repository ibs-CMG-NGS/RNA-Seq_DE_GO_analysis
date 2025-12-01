# Snakemake Bridge: 사용 가이드

## 📋 개요

Bash 스크립트 대신 **Snakemake 기반 통합 워크플로우**를 사용하여 두 파이프라인을 연결합니다.

```
RNA-Seq_DE_GO_analysis  →  Bridge Snakefile  →  RNA-Seq_GO_GSEA_analysis
   (Snakemake)              (자동 변환)            (Snakemake)
```

### 장점

✅ **의존성 관리**: Snakemake가 파일 의존성 자동 추적  
✅ **병렬 처리**: 여러 비교군 동시 변환 가능  
✅ **재현성**: 워크플로우 정의가 코드로 명시  
✅ **부분 재실행**: 변경된 파일만 자동으로 재처리  
✅ **로그 관리**: 모든 단계별 로그 자동 생성  
✅ **통합 환경**: 이미 익숙한 Snakemake 환경 활용  

---

## 🚀 빠른 시작

### ⚠️ 중요 사항

1. **`--cores` 옵션 필수**: Snakemake는 항상 `--cores` 옵션 필요
2. **`--use-conda` 옵션 권장**: Python 패키지 자동 관리

```bash
# ✓ 올바른 사용 (conda 환경 자동 관리)
snakemake -s bridge/Snakefile --use-conda --cores 2

# ✓ 올바른 사용 (수동 패키지 관리)
snakemake -s bridge/Snakefile --cores 2

# ✗ 오류 발생
snakemake -s bridge/Snakefile
# Error: cores have to be specified for local execution
```

### 0. 첫 실행: 설정 확인 (권장)

```bash
# DE 분석 파이프라인 루트에서 실행
cd /path/to/RNA-Seq_DE_GO_analysis

# 현재 설정 확인
snakemake -s bridge/Snakefile show_config --cores 1
```

**출력 예시:**
```
======================================================================
Bridge Pipeline Configuration
======================================================================
Working Directory:   /home/ygkim/ngs_pipeline/RNA-Seq_DE_GO_analysis
DE Pipeline Root:    /home/ygkim/ngs_pipeline/RNA-Seq_DE_GO_analysis
Experiment:          H2O2_Neuron
DE Output Dir:       /home/.../output/H2O2_Neuron
----------------------------------------------------------------------
Pairwise Dir:        /home/.../output/H2O2_Neuron/pairwise
  Exists:            Yes ✓  (또는 No ✗)
----------------------------------------------------------------------
Detected Comparisons: 2
  - H2O2_vs_Control
  - GABA_vs_Control
======================================================================
```

### 1. 비교 그룹 목록 확인

```bash
snakemake -s bridge/Snakefile list_comparisons --cores 1
```

**출력 예시:**
```
Available comparisons (in output/H2O2_Neuron/pairwise/):
  1. GABA_vs_Control
  2. H2O2_vs_Control
```

---

### 2. 특정 비교 그룹 변환

#### 방법 A: Conda 환경 자동 관리 (권장) ✨

```bash
# 예: H2O2_vs_Control 변환
snakemake -s bridge/Snakefile \
  convert_H2O2_vs_Control \
  --use-conda \
  --cores 2
```

**첫 실행 시:**
- Conda 환경 자동 생성 (`bridge_converter`)
- pandas, openpyxl 자동 설치
- 이후 실행부터는 기존 환경 재사용

#### 방법 B: 수동 패키지 관리

```bash
# Python 패키지 직접 설치
pip install pandas openpyxl

# 변환 실행 (--use-conda 없이)
snakemake -s bridge/Snakefile \
  convert_H2O2_vs_Control \
  --cores 2
```

**생성되는 파일:**
```
../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/
  └── H2O2_vs_Control_DE_results.xlsx
```

---

### 3. 모든 비교 그룹 변환

```bash
# 모든 pairwise 비교 자동 변환 (conda 환경 자동 관리)
snakemake -s bridge/Snakefile batch_convert --use-conda --cores 2
```

**생성되는 파일:**
```
../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/
  ├── H2O2_vs_Control_DE_results.xlsx
  └── GABA_vs_Control_DE_results.xlsx
```

---

### 4. Dry-run (실제 실행 전 미리보기)

```bash
# 어떤 작업이 실행될지 확인
snakemake -s bridge/Snakefile --use-conda --cores 2 -n

# 상세 정보와 함께 확인
snakemake -s bridge/Snakefile --use-conda --cores 2 -n -p
```

---

## 📖 상세 사용법

### Conda 환경 관리

#### 자동 환경 생성 (권장)

```bash
# --use-conda 플래그 사용
snakemake -s bridge/Snakefile batch_convert --use-conda --cores 2
```

**장점:**
- ✅ 패키지 버전 자동 관리
- ✅ 재현 가능한 환경
- ✅ 시스템 Python과 격리
- ✅ 한 번 생성하면 재사용

**환경 파일 위치:**
```
bridge/environment.yml
```

**환경 내용:**
```yaml
name: bridge_converter
channels:
  - conda-forge
  - defaults
dependencies:
  - python>=3.8
  - pandas>=1.3
  - openpyxl>=3.0
```

#### 환경 초기화

```bash
# Conda 환경 삭제 (필요 시)
conda env remove -n bridge_converter

# 다음 실행 시 자동 재생성
snakemake -s bridge/Snakefile batch_convert --use-conda --cores 2
```

---

### Rule 목록

| Rule | 설명 | 사용 예시 |
|------|------|----------|
| `batch_convert` (기본) | 모든 비교군 변환 | `snakemake -s bridge/Snakefile --use-conda --cores 2` |
| `convert_{comparison}` | 단일 비교군 변환 | `snakemake -s bridge/Snakefile convert_H2O2_vs_Control --use-conda --cores 2` |
| `show_config` | 설정 표시 | `snakemake -s bridge/Snakefile show_config --cores 1` |
| `list_comparisons` | 비교군 목록 | `snakemake -s bridge/Snakefile list_comparisons --cores 1` |
| `validate_inputs` | 입력 파일 검증 | `snakemake -s bridge/Snakefile validate_inputs --cores 1` |
| `clean_converted_files` | 변환 파일 삭제 | `snakemake -s bridge/Snakefile clean_converted_files --cores 1` |

### 설정 옵션 (--config)

```bash
snakemake -s bridge/Snakefile \
  --config \
    experiment=H2O2_Neuron \              # 실험 디렉토리
    comparison=H2O2_vs_Control \          # 특정 비교군 (선택)
    source_pipeline=deseq2 \              # DE 도구 (deseq2/edger/limma)
    de_output_dir=/custom/path \          # 커스텀 DE 출력 경로
    gsea_input_dir=/custom/gsea/input \   # 커스텀 GSEA 입력 경로
  --cores 1
```

---

## 📁 입력/출력 구조

### 입력 (DE 파이프라인 결과)

```
output/
└── H2O2_Neuron/
    └── pairwise/
        ├── H2O2_vs_Control/
        │   └── final_de_results.csv  ← 입력
        └── GABA_vs_Control/
            └── final_de_results.csv  ← 입력
```

### 출력 (GSEA 파이프라인 입력)

```
../RNA-Seq_GO_GSEA_analysis/
└── data/
    └── from_de_pipeline/
        ├── H2O2_vs_Control_DE_results.xlsx  ← 출력
        └── GABA_vs_Control_DE_results.xlsx  ← 출력
```

### 로그 파일

```
output/
└── H2O2_Neuron/
    └── logs/
        ├── bridge_convert_H2O2_vs_Control.log
        └── bridge_convert_GABA_vs_Control.log
```

---

## 💡 실전 사용 예시

### 예시 1: 표준 워크플로우

```bash
# 1. DE 분석 실행
snakemake --configfile config_H2O2_Neuron.yml --cores 4

# 2. 결과 변환
snakemake -s bridge/Snakefile --cores 1

# 3. GSEA 파이프라인에서 분석
cd ../RNA-Seq_GO_GSEA_analysis
# GSEA Jupyter notebook 또는 Snakefile 실행
```

### 예시 2: 한 번에 실행 (통합 워크플로우)

```bash
# DE 분석부터 변환까지 한 번에
snakemake -s bridge/Snakefile run_full_workflow --cores 4
```

### 예시 3: 디버깅 및 테스트

```bash
# 1. 설정 확인
snakemake -s bridge/Snakefile show_config

# 2. 비교군 목록 확인
snakemake -s bridge/Snakefile list_comparisons

# 3. 입력 파일 검증
snakemake -s bridge/Snakefile validate_inputs

# 4. Dry-run으로 미리보기
snakemake -s bridge/Snakefile -n -p

# 5. 특정 비교군만 테스트
snakemake -s bridge/Snakefile \
  --config comparison=H2O2_vs_Control \
  --cores 1
```

### 예시 4: edgeR 또는 limma 결과 변환

```bash
# edgeR 결과
snakemake -s bridge/Snakefile \
  --config source_pipeline=edger \
  --cores 1

# limma 결과
snakemake -s bridge/Snakefile \
  --config source_pipeline=limma \
  --cores 1
```

### 예시 5: 병렬 처리 (여러 비교군)

```bash
# 4개 코어로 병렬 변환
snakemake -s bridge/Snakefile --cores 4

# 최대 코어 사용
snakemake -s bridge/Snakefile --cores all
```

---

## 🔧 고급 사용법

### 1. 커스텀 경로 설정

```bash
snakemake -s bridge/Snakefile \
  --config \
    de_pipeline_root=/custom/de/path \
    gsea_pipeline_root=/custom/gsea/path \
    de_output_dir=/custom/de/output \
    gsea_input_dir=/custom/gsea/input \
  --cores 1
```

### 2. 특정 파일만 강제 재실행

```bash
# H2O2_vs_Control만 재변환
snakemake -s bridge/Snakefile \
  --forcerun convert_de_to_gsea \
  --config comparison=H2O2_vs_Control \
  --cores 1
```

### 3. 그래프 생성 (워크플로우 시각화)

```bash
# DAG (Directed Acyclic Graph) 생성
snakemake -s bridge/Snakefile --dag | dot -Tpng > workflow_dag.png

# Rule graph 생성
snakemake -s bridge/Snakefile --rulegraph | dot -Tpng > workflow_rules.png
```

### 4. 클러스터 환경에서 실행

```bash
# SLURM 클러스터 예시
snakemake -s bridge/Snakefile \
  --cluster "sbatch -p normal -n 1 -c 1" \
  --jobs 10
```

---

## 🐛 문제 해결

### 문제 1: "cores have to be specified" 오류

**증상**:
```
Error: cores have to be specified for local execution (use --cores N with N being a number >= 1 or 'all')
```

**원인**: Snakemake는 항상 `--cores` 옵션이 필요함

**해결**:
```bash
# ✗ 잘못된 사용
snakemake -s bridge/Snakefile show_config

# ✓ 올바른 사용 (--cores 1 추가)
snakemake -s bridge/Snakefile show_config --cores 1
```

### 문제 2: "No comparisons detected" 또는 "Pairwise directory not found"

**증상**:
```
Warning: Pairwise directory not found: /home/ygkim/ngs_pipeline/output/H2O2_Neuron/pairwise
Detected comparisons: []
```

**원인**: 
- DE 분석이 아직 실행되지 않았음
- `experiment` 설정이 잘못됨
- 다른 경로에서 실행됨

**해결**:

```bash
# 1. 설정 확인 (--cores 1 필수!)
snakemake -s bridge/Snakefile show_config --cores 1

# 2. 출력에서 확인할 사항:
#    - DE Output Dir이 올바른지
#    - Pairwise Dir Exists가 "Yes ✓"인지

# 3. 경로가 틀렸다면 experiment 지정:
snakemake -s bridge/Snakefile \
  --config experiment=H2O2_Neuron \
  show_config --cores 1

# 4. DE 분석이 안 됐다면 먼저 실행:
snakemake --configfile config_H2O2_Neuron.yml --cores 4

# 5. 수동으로 경로 확인:
ls -lh output/H2O2_Neuron/pairwise/
# → H2O2_vs_Control, GABA_vs_Control 등이 보여야 함

ls -lh output/H2O2_Neuron/pairwise/*/final_de_results.csv
# → final_de_results.csv 파일들이 있어야 함
```

### 문제 3: "Input file not found"

**원인**: 특정 비교군의 final_de_results.csv가 없음

**해결**:
```bash
# 입력 파일 검증
snakemake -s bridge/Snakefile validate_inputs --cores 1

# 파일 확인
ls -lh output/H2O2_Neuron/pairwise/*/final_de_results.csv
```

### 문제 4: Python 패키지 오류 (pandas, openpyxl)

**원인**: 필요한 Python 패키지가 설치되지 않음

**해결 방법 1: Conda 자동 관리 (권장) ✨**
```bash
# --use-conda 플래그 추가
snakemake -s bridge/Snakefile batch_convert --use-conda --cores 2

# 첫 실행 시 conda 환경 자동 생성됨
# 이후 실행부터는 기존 환경 재사용
```

**해결 방법 2: 수동 설치**
```bash
# pip로 직접 설치
pip install pandas openpyxl

# 또는 conda로 설치
conda install pandas openpyxl

# 설치 확인
python -c "import pandas; import openpyxl; print('OK')"
```

### 문제 5: Python 스크립트 오류

**원인**: convert_de_to_gsea.py에서 오류 발생

**해결**:
```bash
# 로그 확인
cat output/H2O2_Neuron/logs/bridge_convert_*.log

# 스크립트 직접 테스트
python3 bridge/convert_de_to_gsea.py --help
```

### 문제 4: 권한 문제

**원인**: GSEA 파이프라인 디렉토리 쓰기 권한 없음

**해결**:
```bash
# 디렉토리 생성 및 권한 확인
mkdir -p ../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline
ls -ld ../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline
```

---

## 📊 성능 최적화

### 병렬 처리 가이드

```bash
# 비교군 수에 따른 권장 코어 수
# 2-5개 비교군: --cores 2
# 6-10개: --cores 4
# 11-20개: --cores 8
# 20개 이상: --cores all

# 예: 8개 비교군, 4코어 병렬
snakemake -s bridge/Snakefile --cores 4
```

### 대용량 데이터 처리

```bash
# batch_convert rule 사용 (매우 많은 비교군)
snakemake -s bridge/Snakefile batch_convert --cores 1
```

---

## 🔗 다른 파이프라인과 통합

### 기존 Snakefile에 통합

메인 `Snakefile`에 다음 추가:

```python
# 브릿지 레이어 포함
include: "bridge/Snakefile"

# 전체 워크플로우에 추가
rule all:
    input:
        # 기존 DE 분석 출력
        expand("output/{experiment}/pairwise/{comparison}/final_de_results.csv",
               experiment=EXPERIMENTS, comparison=COMPARISONS),
        # 브릿지 변환 출력 추가
        expand("../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/{comparison}_DE_results.xlsx",
               comparison=COMPARISONS)
```

---

## 📝 로그 및 디버깅

### 로그 확인

```bash
# 모든 브릿지 로그 확인
ls -lh output/H2O2_Neuron/logs/bridge_*.log

# 특정 비교군 로그 확인
cat output/H2O2_Neuron/logs/bridge_convert_H2O2_vs_Control.log

# 실시간 모니터링
tail -f output/H2O2_Neuron/logs/bridge_convert_H2O2_vs_Control.log
```

### 디버깅 모드

```bash
# 상세 출력 모드
snakemake -s bridge/Snakefile --cores 1 -p -r

# 특정 rule 디버깅
snakemake -s bridge/Snakefile \
  --debug-dag \
  --printshellcmds \
  --cores 1
```

---

## ✅ 체크리스트

변환 전:
- [ ] DE 분석 완료
- [ ] final_de_results.csv 파일 존재 확인
- [ ] GSEA 파이프라인 경로 확인
- [ ] Python 환경 확인

변환 실행:
- [ ] `show_config`로 설정 확인
- [ ] `validate_inputs`로 입력 검증
- [ ] Dry-run (`-n`) 먼저 실행
- [ ] 실제 변환 실행

변환 후:
- [ ] 출력 Excel 파일 존재 확인
- [ ] 로그 파일 확인
- [ ] Excel 파일 내용 검증
- [ ] GSEA 파이프라인에서 사용

---

## 📚 관련 문서

- **Python 스크립트 가이드**: [convert_de_to_gsea.py 문서](README.md)
- **사용 예시**: [EXAMPLES.sh](EXAMPLES.sh)
- **워크플로우 다이어그램**: [WORKFLOW_DIAGRAM.md](WORKFLOW_DIAGRAM.md)
- **Snakemake 공식 문서**: https://snakemake.readthedocs.io/

---

## 🎯 요약

### 가장 많이 사용하는 명령어

```bash
# 1. 기본 사용 (모든 비교군 변환)
snakemake -s bridge/Snakefile --cores 1

# 2. 설정 확인
snakemake -s bridge/Snakefile show_config

# 3. Dry-run
snakemake -s bridge/Snakefile -n

# 4. 특정 비교군
snakemake -s bridge/Snakefile --config comparison=H2O2_vs_Control --cores 1

# 5. 병렬 처리
snakemake -s bridge/Snakefile --cores 4
```

**이제 Bash 스크립트 없이 Snakemake만으로 파이프라인을 연결할 수 있습니다!** 🎉
