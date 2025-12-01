# Bridge Layer: Snakemake Integration Complete! 🎉

## ✅ 업데이트 완료

Bash 스크립트 대신 **Snakemake 기반 통합 워크플로우**로 전환했습니다!

---

## 📁 새로 추가된 파일

```
bridge/
├── Snakefile ⭐ NEW!              # Snakemake 통합 워크플로우
├── SNAKEMAKE_GUIDE.md ⭐ NEW!     # Snakemake 사용 가이드 (상세)
├── convert_de_to_gsea.py          # Python 변환 스크립트 (기존)
├── run_downstream_analysis.sh     # Bash wrapper (대체 방법)
├── README.md                      # Python 스크립트 가이드
├── EXAMPLES.sh                    # 사용 예제 모음
├── WORKFLOW_DIAGRAM.md            # 워크플로우 다이어그램
├── IMPLEMENTATION_SUMMARY.md      # 기술 문서
└── INDEX.md                       # 네비게이션 (업데이트됨)
```

---

## 🚀 사용 방법

### ⭐ 방법 1: Snakemake (권장)

```bash
# 1. 모든 비교군 자동 변환
snakemake -s bridge/Snakefile --cores 1

# 2. 특정 비교군만 변환
snakemake -s bridge/Snakefile \
  --config comparison=H2O2_vs_Control experiment=H2O2_Neuron \
  --cores 1

# 3. 병렬 처리 (4개 코어)
snakemake -s bridge/Snakefile --cores 4

# 4. 설정 확인
snakemake -s bridge/Snakefile show_config

# 5. 비교군 목록 확인
snakemake -s bridge/Snakefile list_comparisons

# 6. Dry-run (실행 전 미리보기)
snakemake -s bridge/Snakefile -n
```

### 방법 2: Python 스크립트 (필요시)

```bash
# Python 스크립트로 직접 실행 (기존 방법 유지)
cd bridge
python3 convert_de_to_gsea.py --batch \
  --de-output-dir ../output/H2O2_Neuron \
  --gsea-input-dir ../../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline
```

---

## ✨ Snakemake 방식의 장점

| 특징 | Bash Script | Snakemake ⭐ |
|------|-------------|--------------|
| 의존성 관리 | 수동 | 자동 ✓ |
| 병렬 처리 | 불가능 | 가능 ✓ |
| 부분 재실행 | 전체 재실행 필요 | 변경된 파일만 ✓ |
| 로그 관리 | 수동 설정 | 자동 생성 ✓ |
| 통합성 | 별도 스크립트 | 파이프라인 통합 ✓ |
| 학습 곡선 | 새로운 문법 | 이미 사용 중 ✓ |

---

## 📊 Snakefile 주요 기능

### 1. 자동 비교군 탐지

```python
# DE 분석 output 디렉토리에서 자동으로 비교군 탐지
COMPARISONS = get_comparisons()
# → ['H2O2_vs_Control', 'GABA_vs_Control', ...]
```

### 2. Rule 기반 워크플로우

```python
rule convert_de_to_gsea:
    input:  final_de_results.csv
    output: {comparison}_DE_results.xlsx
    shell:  python3 convert_de_to_gsea.py ...
```

### 3. 유용한 유틸리티 Rules

- `show_config`: 현재 설정 표시
- `list_comparisons`: 감지된 비교군 목록
- `validate_inputs`: 입력 파일 검증
- `clean_converted_files`: 변환 파일 정리

### 4. 병렬 처리 지원

```bash
# 4개 비교군을 4개 코어로 동시 처리
snakemake -s bridge/Snakefile --cores 4
```

### 5. 자동 로그 생성

```
output/H2O2_Neuron/logs/
├── bridge_convert_H2O2_vs_Control.log
└── bridge_convert_GABA_vs_Control.log
```

---

## 🔄 전체 워크플로우

### 통합 워크플로우 예시

```bash
# Step 1: DE 분석 실행
snakemake --configfile config_H2O2_Neuron.yml --cores 4

# Step 2: 결과 변환 (Bridge)
snakemake -s bridge/Snakefile --cores 4

# Step 3: GSEA 분석 실행
cd ../RNA-Seq_GO_GSEA_analysis
snakemake -s workflow/Snakefile_GO --cores 4
```

### 한 줄로 실행 (향후 가능)

```bash
# 메인 Snakefile에 bridge 포함 시
snakemake --cores 4  # DE 분석 + 변환까지 한 번에!
```

---

## 📖 문서 구조

### 읽는 순서 (권장)

1. **`SNAKEMAKE_GUIDE.md`** ⭐
   - Snakemake 사용법
   - Rule 설명
   - 예제 명령어
   - 문제 해결

2. **`INDEX.md`**
   - 빠른 네비게이션
   - 파일 설명
   - 학습 경로

3. **`WORKFLOW_DIAGRAM.md`**
   - 시각적 다이어그램
   - 데이터 흐름
   - 아키텍처

4. **`README.md`** (필요시)
   - Python 스크립트 직접 사용
   - 세부 옵션 설명

---

## 🎯 실전 사용 시나리오

### 시나리오 1: 일반적인 사용

```bash
# DE 분석 완료 후
cd /path/to/RNA-Seq_DE_GO_analysis

# 모든 결과 변환
snakemake -s bridge/Snakefile --cores 1

# 출력 확인
ls -lh ../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline/
```

### 시나리오 2: 디버깅

```bash
# 설정 확인
snakemake -s bridge/Snakefile show_config

# 비교군 확인
snakemake -s bridge/Snakefile list_comparisons

# 입력 검증
snakemake -s bridge/Snakefile validate_inputs

# Dry-run
snakemake -s bridge/Snakefile -n
```

### 시나리오 3: 특정 비교군만

```bash
# H2O2_vs_Control만 재변환
snakemake -s bridge/Snakefile \
  --forcerun convert_de_to_gsea \
  --config comparison=H2O2_vs_Control \
  --cores 1
```

### 시나리오 4: 대량 병렬 처리

```bash
# 20개 비교군을 8코어로 병렬 처리
snakemake -s bridge/Snakefile --cores 8
```

---

## 🔧 커스터마이징

### 기본 설정 변경

Snakefile 상단의 설정 변경:

```python
# 기본 경로 수정
DE_PIPELINE_ROOT = "/custom/path/to/de_pipeline"
GSEA_PIPELINE_ROOT = "/custom/path/to/gsea_pipeline"

# 실험 디렉토리 변경
EXPERIMENT = "My_Experiment"
```

### 명령줄에서 설정

```bash
snakemake -s bridge/Snakefile \
  --config \
    experiment=Custom_Experiment \
    de_output_dir=/custom/de/output \
    gsea_input_dir=/custom/gsea/input \
  --cores 1
```

---

## 📝 주요 변경사항

### README.md 업데이트

- ✅ Snakemake 방식을 "방법 1 (권장)"으로 강조
- ✅ Python 스크립트는 "방법 2 (대체)"로 유지
- ✅ Bash 스크립트는 참고용으로 유지
- ✅ SNAKEMAKE_GUIDE.md 링크 추가

### INDEX.md 업데이트

- ✅ Snakefile을 Core Files 최상단에 배치
- ✅ Learning Path에 Snakemake 사용자 섹션 추가
- ✅ Quick Commands에 Snakemake 명령어 우선 표시

---

## 🎓 학습 리소스

### 초급 (Snakemake 처음 사용)

1. `SNAKEMAKE_GUIDE.md` - 빠른 시작 섹션
2. 실습: `snakemake -s bridge/Snakefile show_config`
3. 실습: `snakemake -s bridge/Snakefile -n`
4. 실습: `snakemake -s bridge/Snakefile --cores 1`

### 중급 (Snakemake 익숙)

1. `SNAKEMAKE_GUIDE.md` - 고급 사용법 섹션
2. Rule 커스터마이징
3. 병렬 처리 최적화
4. 메인 파이프라인에 통합

### 고급 (파이프라인 개발자)

1. Snakefile 코드 분석
2. 새로운 Rule 추가
3. GSEA 파이프라인과 완전 통합
4. 클러스터 환경 최적화

---

## ✅ 체크리스트

### Snakemake 방식 사용 시

- [x] Snakefile 생성 완료
- [x] SNAKEMAKE_GUIDE.md 작성 완료
- [x] INDEX.md 업데이트 완료
- [x] README.md 업데이트 완료
- [ ] 실제 데이터로 테스트
- [ ] 메인 Snakefile에 통합 (선택사항)

### 사용자 액션 필요

- [ ] `snakemake -s bridge/Snakefile show_config` 실행
- [ ] `snakemake -s bridge/Snakefile list_comparisons` 실행
- [ ] `snakemake -s bridge/Snakefile -n` 실행 (dry-run)
- [ ] `snakemake -s bridge/Snakefile --cores 1` 실행 (실제 변환)
- [ ] 변환된 Excel 파일 검증
- [ ] GSEA 파이프라인에서 사용

---

## 🎉 결론

**Bash 스크립트에서 Snakemake 기반 워크플로우로 성공적으로 전환했습니다!**

### 주요 개선사항

✅ **통합성**: 두 Snakemake 파이프라인이 자연스럽게 연결  
✅ **자동화**: 의존성 추적 및 자동 재실행  
✅ **효율성**: 병렬 처리로 처리 시간 단축  
✅ **일관성**: 모든 파이프라인이 동일한 도구 사용  
✅ **유지보수**: Rule 기반으로 관리 용이  

### 다음 단계

1. 실제 데이터로 테스트
2. 성능 모니터링
3. 필요시 Rule 최적화
4. 메인 파이프라인 통합 고려

---

**문서 위치**:
- **Snakemake 가이드**: `bridge/SNAKEMAKE_GUIDE.md` ⭐
- **빠른 네비게이션**: `bridge/INDEX.md`
- **Python 스크립트**: `bridge/README.md`
- **워크플로우 다이어그램**: `bridge/WORKFLOW_DIAGRAM.md`

**이제 Snakemake만으로 두 파이프라인을 완벽하게 연결할 수 있습니다!** 🚀
