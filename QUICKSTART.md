# Quick Start Guide - CLI 빠른 시작 가이드

RNA-Seq 분석 파이프라인을 CLI로 실행하기 위한 빠른 시작 가이드입니다.

## 📋 사전 준비

### 필수 요구사항
- R (버전 4.2 이상)
- Rscript 실행 가능
- 충분한 메모리 (최소 4GB 권장)

### 선택 사항 (워크플로우 관리 도구)
- Snakemake (Python 기반)
- Nextflow (Java 기반)

## 🚀 5분 안에 시작하기

### 1단계: 프로젝트 클론
```bash
git clone https://github.com/ibs-CMG-NGS/RNA-Seq_DE_GO_analysis
cd RNA-Seq_DE_GO_analysis
```

### 2단계: Conda 환경 설정 (최초 1회)
```bash
# Conda 환경 생성 (모든 의존성 자동 설치)
conda env create -f environment.yml

# 환경 활성화
conda activate rna-seq-de-go-analysis
```

### 3단계: 파이프라인 실행
```bash
# 방법 A: 전체 파이프라인 실행
Rscript run_pipeline.R

# 방법 B: Shell 스크립트 사용
./run_pipeline.sh
```

### 4단계: 결과 확인
```bash
# 결과는 output/YYYY-MM-DD_HH-MM-SS/ 디렉토리에 저장됩니다
ls -lh output/*/
```

## 📝 설정 파일 커스터마이징

분석을 시작하기 전에 `config.yml`을 편집하여 자신의 데이터에 맞게 수정하세요:

```bash
# config.yml 편집
nano config.yml
# 또는
vim config.yml
```

주요 설정 항목:
- `count_data_path`: 카운트 데이터 파일 경로
- `metadata_path`: 메타데이터 파일 경로
- `species`: "human" 또는 "mouse"
- `de_analysis.method`: "DESeq2", "edgeR", "limma-voom" 중 선택

## 🔧 고급 사용법

### 특정 분석 단계만 실행
```bash
# DE 분석만
Rscript run_pipeline.R --step de

# 플롯 생성만
Rscript run_pipeline.R --step plots

# Enrichment 분석만
Rscript run_pipeline.R --step enrichment

# GO 플롯만
Rscript run_pipeline.R --step go_plots
```

### 커스텀 출력 디렉토리 사용
```bash
Rscript run_pipeline.R --output my_results
```

### 여러 실험 배치 실행
```bash
#!/bin/bash
for config in configs/*.yml; do
    name=$(basename $config .yml)
    Rscript run_pipeline.R \
        --config $config \
        --output results/$name
done
```

## 🐍 Snakemake로 실행

### 설치 (conda 환경 권장)
```bash
conda install -c bioconda snakemake
```

### 실행
```bash
# Dry-run
snakemake --cores 4 --dry-run

# 실제 실행
snakemake --cores 4
```

## 🌊 Nextflow로 실행

### 설치
```bash
# Java 필요
curl -s https://get.nextflow.io | bash
```

### 실행
```bash
nextflow run nextflow_pipeline.nf
```

## 📊 예상 결과물

파이프라인 실행 후 다음 파일들이 생성됩니다:

```
output/YYYY-MM-DD_HH-MM-SS/
├── final_de_results.csv         # DE 분석 결과
├── final_de_results.xlsx        # Excel 형식 결과
├── pca_plot.png                 # PCA 플롯
├── volcano_plot.png             # Volcano 플롯
├── go_enrichment_*.csv          # GO 분석 결과
├── go_dotplot_*.png             # GO Dot plot
├── kegg_enrichment_*.csv        # KEGG 분석 결과
├── kegg_dotplot_*.png           # KEGG Dot plot
├── go_barplot_*.png             # GO Bar plot
└── config_used.yml              # 사용된 설정 파일
```

## ⚠️ 문제 해결

### 1. "Rscript: command not found"
```bash
# R 설치 확인
which R

# PATH에 R 추가 (예시)
export PATH="/usr/local/bin:$PATH"
```

### 2. "package 'optparse' is not available"
```bash
# Conda 환경이 활성화되었는지 확인
conda activate rna-seq-de-go-analysis

# 또는 환경 재생성
conda env remove -n rna-seq-de-go-analysis
conda env create -f environment.yml
```

### 3. 메모리 부족 오류
```bash
# R 메모리 한도 증가
export R_MAX_VSIZE=16Gb
Rscript run_pipeline.R
```

### 4. 권한 오류
```bash
# 실행 권한 부여
chmod +x run_pipeline.R run_pipeline.sh
```

## 🧪 테스트 실행

설치가 제대로 되었는지 확인:
```bash
./test_cli.sh
```

## 📚 추가 자료

- 상세 CLI 가이드: [CLI_GUIDE.md](CLI_GUIDE.md)
- 전체 README: [README.md](README.md)
- 설정 파일 예제: [config.yml](config.yml)

## 💡 팁

1. **첫 실행 시간**: Conda 환경 생성은 5-10분 정도 소요될 수 있으며, 첫 분석은 10-20분 정도 소요됩니다.

2. **병렬 처리**: Snakemake나 Nextflow를 사용하면 여러 샘플을 병렬로 처리할 수 있어 시간을 절약할 수 있습니다.

3. **재현성**: 각 실행마다 사용된 `config.yml`이 결과 디렉토리에 복사되므로 나중에 분석을 재현할 수 있습니다.

4. **디버깅**: 문제가 발생하면 `--verbose` 옵션을 사용하여 상세 로그를 확인하세요:
   ```bash
   Rscript run_pipeline.R --verbose
   ```

5. **단계별 실행**: 각 단계를 따로 실행하면 문제가 발생한 부분을 쉽게 파악할 수 있습니다.

## 🎯 다음 단계

분석이 완료되면:
1. 결과 파일들을 검토하세요
2. Volcano plot에서 유의미한 유전자를 확인하세요
3. GO/KEGG 분석 결과로 생물학적 의미를 해석하세요
4. 필요하면 `config.yml`을 조정하여 재분석하세요

---

궁금한 점이 있으면 GitHub Issues에 질문해 주세요!
