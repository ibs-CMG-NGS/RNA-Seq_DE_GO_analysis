# CLI 사용 가이드 (Command Line Interface Guide)

이 문서는 RNA-Seq 분석 파이프라인을 커맨드 라인에서 실행하는 방법을 상세히 설명합니다.

## 목차
- [기본 사용법](#기본-사용법)
- [옵션 상세 설명](#옵션-상세-설명)
- [사용 예제](#사용-예제)
- [Snakemake 통합](#snakemake-통합)
- [문제 해결](#문제-해결)

## 기본 사용법

### 최초 실행 (의존성 설치 포함)
```bash
Rscript run_pipeline.R --install-deps
```

### 표준 실행
```bash
# R 스크립트로 직접 실행
Rscript run_pipeline.R

# 또는 Shell 스크립트 사용
./run_pipeline.sh
```

## 옵션 상세 설명

### `--config, -c <파일경로>`
- 사용할 설정 파일의 경로를 지정합니다
- 기본값: `config.yml`
- 예제: `--config experiments/exp1_config.yml`

### `--output, -o <디렉토리경로>`
- 분석 결과를 저장할 디렉토리를 지정합니다
- 기본값: `output/YYYY-MM-DD_HH-MM-SS` (타임스탬프 자동 생성)
- 예제: `--output results/experiment_01`

### `--step, -s <단계명>`
- 실행할 분석 단계를 선택합니다
- 가능한 값:
  - `all`: 모든 단계 실행 (기본값)
  - `de`: Differential Expression 분석만
  - `plots`: PCA, Volcano plot 생성만
  - `enrichment`: GO/KEGG enrichment 분석만
  - `go_plots`: GO bar plot 생성만
- 예제: `--step de`

### `--install-deps`
- 분석에 필요한 모든 R 패키지를 자동으로 설치합니다
- 최초 실행 시 또는 새로운 환경에서 실행할 때 사용합니다
- 예제: `--install-deps`

### `--verbose, -v`
- 상세한 로그를 출력합니다 (디버깅 용도)
- 예제: `--verbose`

### `--help`
- 도움말을 표시합니다
- 예제: `--help`

## 사용 예제

### 1. 기본 전체 파이프라인 실행
```bash
Rscript run_pipeline.R
```
- config.yml 사용
- output/YYYY-MM-DD_HH-MM-SS 폴더에 결과 저장
- 모든 분석 단계 실행

### 2. 커스텀 설정 파일로 실행
```bash
Rscript run_pipeline.R --config my_experiment.yml
```

### 3. 특정 출력 디렉토리 지정
```bash
Rscript run_pipeline.R --output results/mouse_liver_analysis
```

### 4. DE 분석만 실행
```bash
Rscript run_pipeline.R --step de --output results/de_only
```

### 5. 이전 DE 분석 결과로 enrichment만 실행
```bash
# 먼저 DE 분석 실행
Rscript run_pipeline.R --step de --output results/my_analysis

# 동일한 출력 디렉토리에 enrichment 추가
Rscript run_pipeline.R --step enrichment --output results/my_analysis
```

### 6. Verbose 모드로 디버깅
```bash
Rscript run_pipeline.R --verbose --step de
```

### 7. 여러 옵션 조합
```bash
Rscript run_pipeline.R \
  --config experiments/exp1.yml \
  --output results/exp1_results \
  --step all \
  --verbose
```

## Snakemake 통합

프로젝트에 포함된 `Snakefile`을 사용하여 파이프라인을 실행할 수 있습니다.

### 기본 실행
```bash
# Dry-run (실행 계획 확인)
snakemake --cores 4 --dry-run

# 실제 실행
snakemake --cores 4
```

### 특정 단계까지만 실행
```bash
# DE 분석만
snakemake --cores 4 de_analysis

# Enrichment까지
snakemake --cores 4 enrichment
```

### 워크플로우 시각화
```bash
# DAG (Directed Acyclic Graph) 생성
snakemake --dag | dot -Tpng > workflow_diagram.png

# 규칙 그래프 생성
snakemake --rulegraph | dot -Tpng > rules_diagram.png
```

### 로그 확인
```bash
# 로그는 logs/ 디렉토리에 저장됩니다
cat logs/*_de_analysis.log
cat logs/*_enrichment.log
```

### 결과 정리
```bash
snakemake clean
```

## Nextflow 통합

Nextflow를 선호하는 경우 `nextflow_pipeline.nf`를 사용할 수 있습니다.

### 기본 실행
```bash
# 전체 파이프라인 실행
nextflow run nextflow_pipeline.nf

# 리포트 생성
nextflow run nextflow_pipeline.nf -with-report report.html
```

### 파라미터 지정
```bash
# 커스텀 config 사용
nextflow run nextflow_pipeline.nf --config my_config.yml

# 출력 디렉토리 지정
nextflow run nextflow_pipeline.nf --outdir results/exp1

# 의존성 설치 포함
nextflow run nextflow_pipeline.nf --install_deps true
```

### 재시작 기능
```bash
# 실패한 작업만 재실행
nextflow run nextflow_pipeline.nf -resume
```

## Snakemake 고급 설정

### 특정 config 사용
Snakefile을 수정하여 CONFIGFILE 변수를 변경하거나, 명령줄에서 지정:
```bash
snakemake --config configfile=my_config.yml --cores 4
```

### 재시작 기능
실패한 작업만 재실행:
```bash
snakemake --cores 4 --rerun-incomplete
```

### 클러스터에서 실행
```bash
snakemake --cluster "qsub -cwd -V" --jobs 10
```

## 문제 해결

### 1. "command not found: Rscript"
**원인**: R이 설치되지 않았거나 PATH에 없습니다.

**해결**:
```bash
# R 설치 확인
which R

# PATH 설정 (필요시)
export PATH="/usr/local/bin:$PATH"
```

### 2. "cannot open file 'run_pipeline.R'"
**원인**: 스크립트가 실행되는 위치 문제

**해결**:
```bash
# 프로젝트 루트 디렉토리로 이동
cd /path/to/RNA-Seq_DE_GO_analysis
./run_pipeline.sh
```

### 3. "package 'optparse' is not available"
**원인**: 필수 패키지가 설치되지 않음

**해결**:
```bash
Rscript run_pipeline.R --install-deps
```

### 4. "Config file not found"
**원인**: config.yml 파일을 찾을 수 없음

**해결**:
```bash
# 절대 경로 사용
Rscript run_pipeline.R --config /full/path/to/config.yml

# 또는 프로젝트 루트에서 실행
cd /path/to/RNA-Seq_DE_GO_analysis
Rscript run_pipeline.R
```

### 5. Snakemake 실행 오류
**원인**: Python 환경 또는 Snakemake 미설치

**해결**:
```bash
# Snakemake 설치 (conda 사용)
conda install -c bioconda snakemake

# 또는 pip 사용
pip install snakemake
```

## 성능 최적화

### 병렬 처리
Snakemake를 사용하면 독립적인 작업을 병렬로 실행할 수 있습니다:
```bash
# 4개 코어 사용
snakemake --cores 4

# 모든 가용 코어 사용
snakemake --cores all
```

### 메모리 관리
대용량 데이터 처리 시:
```bash
# R 메모리 제한 증가 (예: 16GB)
export R_MAX_VSIZE=16Gb
Rscript run_pipeline.R
```

## 배치 작업 예제

여러 실험을 순차적으로 실행:
```bash
#!/bin/bash
# batch_analysis.sh

experiments=("exp1" "exp2" "exp3")

for exp in "${experiments[@]}"; do
  echo "Running analysis for $exp"
  Rscript run_pipeline.R \
    --config "configs/${exp}_config.yml" \
    --output "results/${exp}" \
    2>&1 | tee "logs/${exp}.log"
done
```

## 추가 리소스

- [R 공식 문서](https://www.r-project.org/)
- [Snakemake 튜토리얼](https://snakemake.readthedocs.io/)
- [optparse R 패키지](https://cran.r-project.org/web/packages/optparse/)

## 문의 및 기여

버그 리포트나 기능 제안은 GitHub Issues를 통해 제출해 주세요.
