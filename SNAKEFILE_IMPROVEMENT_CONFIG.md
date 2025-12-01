# Snakefile 개선: Config 파일 중앙 관리

## 🎯 개선 사항

### 문제점
이전에는 `config_H2O2_Neuron.yml` 파일명이 Snakefile의 여러 rule에 하드코딩되어 있었습니다:
- `rule run_omnibus_test` → `config_file = "config_H2O2_Neuron.yml"`
- `rule run_pairwise_de` → `config_file = "config_H2O2_Neuron.yml"`
- `rule generate_global_pca` → `config_file = "config_H2O2_Neuron.yml"`
- ... (총 7개 rule)

**문제**:
- ❌ 다른 config 파일 사용 시 7군데를 모두 수정해야 함
- ❌ 실수로 일부만 수정하면 오류 발생
- ❌ 번거롭고 비효율적

### 해결 방법

**상단에 변수 선언**:
```python
# Snakefile (Line 6)
CONFIG_FILE = "config_H2O2_Neuron.yml"  # ← 여기만 수정!

configfile: CONFIG_FILE
```

**모든 rule에서 변수 사용**:
```python
rule run_omnibus_test:
    input:
        config_file = CONFIG_FILE  # ← 변수 참조

rule run_pairwise_de:
    input:
        config_file = CONFIG_FILE  # ← 변수 참조
        
# ... 모든 rule에 동일하게 적용
```

## ✅ 장점

1. **편의성**: 한 곳만 수정하면 됨 (6번째 줄)
2. **안전성**: 실수로 일부만 수정할 위험 제거
3. **유지보수**: 코드 관리가 훨씬 쉬워짐
4. **가독성**: Config 파일 위치가 명확함

## 📝 사용 방법

### 다른 프로젝트 분석하기

```python
# Snakefile 6번째 줄만 수정
CONFIG_FILE = "config_Shank2.yml"  # 변경!
```

그리고 실행:
```bash
snakemake --cores 4 --use-conda
```

### 파라미터 테스트

```python
# prefilter 값 비교 테스트
CONFIG_FILE = "config_test_prefilter0.yml"   # 테스트 1
# → 실행 후
CONFIG_FILE = "config_test_prefilter10.yml"  # 테스트 2
```

### 여러 데이터셋 배치 분석

```bash
# 방법 1: Snakefile 수정 + 실행 반복
# Snakefile → CONFIG_FILE = "config_A.yml"
snakemake --cores 4 --use-conda
# Snakefile → CONFIG_FILE = "config_B.yml"  
snakemake --cores 4 --use-conda

# 방법 2: 명령줄 옵션 (고급)
snakemake --cores 4 --use-conda --configfile config_A.yml
snakemake --cores 4 --use-conda --configfile config_B.yml
```

## 🔧 수정된 Rule 목록

총 **7개 rule**이 업데이트되었습니다:

1. ✅ `run_omnibus_test`
2. ✅ `run_pairwise_de`
3. ✅ `generate_global_pca`
4. ✅ `generate_pairwise_volcano`
5. ✅ `go_enrichment`
6. ✅ `kegg_enrichment`
7. ✅ `go_barplots`
8. ✅ `generate_go_summary_table`

모두 `config_file = CONFIG_FILE`로 변경되었습니다.

## 📚 문서 업데이트

**README.md** 섹션 추가:
- "🔧 Config 파일 변경하기" 섹션 신설
- 이전 방식과 현재 방식 비교 설명
- 사용 예시 추가

## 🎓 기술 세부사항

### Python 변수 활용
```python
# 상단에서 정의
CONFIG_FILE = "config_H2O2_Neuron.yml"

# configfile directive에서 사용
configfile: CONFIG_FILE

# 각 rule의 input에서 사용
input:
    config_file = CONFIG_FILE
```

### Snakemake의 변수 스코프
- `CONFIG_FILE`은 전역 변수로 모든 rule에서 접근 가능
- Python 문법을 그대로 사용 가능
- 다른 유용한 전역 변수들:
  - `OUTPUT_DIR` (이미 사용 중)
  - `R_ENV_NAME` (이미 사용 중)
  - `CONFIG_FILE` (새로 추가)

## 🚀 향후 개선 가능사항

1. **환경 변수로 지정**:
   ```bash
   export PIPELINE_CONFIG="config_A.yml"
   snakemake --cores 4 --use-conda
   ```

2. **Snakemake 프로파일**:
   ```yaml
   # profiles/default/config.yaml
   configfile: config_H2O2_Neuron.yml
   cores: 4
   use-conda: true
   ```

3. **자동 감지**:
   ```python
   import os
   CONFIG_FILE = os.getenv("PIPELINE_CONFIG", "config_H2O2_Neuron.yml")
   ```

## ✅ 검증

### Before (7번 수정 필요)
```bash
grep -n "config_H2O2_Neuron.yml" Snakefile
# 출력: 7개 라인 발견
```

### After (1번만 수정)
```bash
grep -n "config_H2O2_Neuron.yml" Snakefile
# 출력: Line 6: CONFIG_FILE = "config_H2O2_Neuron.yml"
```

모든 다른 출현은 `CONFIG_FILE` 변수로 대체됨!

---

**날짜**: 2025-12-01  
**개선 내용**: Config 파일 경로 중앙 집중화  
**수정 파일**: `Snakefile`, `README.md`  
**영향 받는 Rule**: 8개 (모든 분석 rule)
