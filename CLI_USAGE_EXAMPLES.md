# CLI 사용 예제 모음

이 문서는 RNA-Seq 분석 파이프라인의 다양한 CLI 사용 예제를 제공합니다.

## 📋 목차
- [기본 사용 예제](#기본-사용-예제)
- [분석 단계별 실행](#분석-단계별-실행)
- [배치 처리 예제](#배치-처리-예제)
- [워크플로우 관리 도구](#워크플로우-관리-도구)
- [HPC 클러스터에서 실행](#hpc-클러스터에서-실행)

---

## 기본 사용 예제

### 예제 1: 표준 실행
```bash
# config.yml을 사용하여 전체 파이프라인 실행
Rscript run_pipeline.R

# 또는 shell 스크립트 사용
./run_pipeline.sh
```

### 예제 2: 도움말 보기
```bash
Rscript run_pipeline.R --help
```

### 예제 3: Verbose 모드로 실행
```bash
# 상세 로그 출력
Rscript run_pipeline.R --verbose
```

---

## 분석 단계별 실행

### 예제 4: DE 분석만 실행
```bash
Rscript run_pipeline.R --step de
```

### 예제 5: 플롯 생성만 실행
```bash
# 이전에 DE 분석이 완료되어 있어야 함
Rscript run_pipeline.R --step plots --output results/previous_run
```

### 예제 6: Enrichment 분석만 실행
```bash
Rscript run_pipeline.R --step enrichment --output results/previous_run
```

### 예제 7: GO 플롯만 생성
```bash
Rscript run_pipeline.R --step go_plots --output results/previous_run
```

---

## 커스텀 설정 사용

### 예제 8: 커스텀 config 파일 사용
```bash
Rscript run_pipeline.R --config experiments/liver_vs_kidney.yml
```

### 예제 9: 특정 출력 디렉토리 지정
```bash
Rscript run_pipeline.R --output results/mouse_experiment_001
```

### 예제 10: 여러 옵션 조합
```bash
Rscript run_pipeline.R \
    --config custom_configs/exp1.yml \
    --output results/exp1_output \
    --step all \
    --verbose
```

---

## 배치 처리 예제

### 예제 11: 여러 실험 순차 실행
```bash
#!/bin/bash
# batch_run.sh

configs=("exp1.yml" "exp2.yml" "exp3.yml")

for config in "${configs[@]}"; do
    echo "Processing $config..."
    Rscript run_pipeline.R \
        --config "configs/$config" \
        --output "results/$(basename $config .yml)"
done
```

### 예제 12: 병렬 실행 (GNU Parallel 사용)
```bash
# 4개 작업을 병렬로 실행
ls configs/*.yml | parallel -j 4 \
    'Rscript run_pipeline.R --config {} --output results/{/.}'
```

### 예제 13: 조건부 실행
```bash
#!/bin/bash
# conditional_run.sh

if [ ! -f "results/exp1/final_de_results.csv" ]; then
    echo "Running DE analysis..."
    Rscript run_pipeline.R --step de --output results/exp1
fi

if [ ! -f "results/exp1/pca_plot.png" ]; then
    echo "Generating plots..."
    Rscript run_pipeline.R --step plots --output results/exp1
fi
```

---

## 워크플로우 관리 도구

### Snakemake 예제

#### 예제 14: Snakemake 기본 실행
```bash
# Dry-run
snakemake --cores 4 --dry-run

# 실제 실행
snakemake --cores 4
```

#### 예제 15: 특정 타겟만 실행
```bash
# DE 분석만
snakemake --cores 4 de_analysis

# Enrichment까지
snakemake --cores 4 enrichment
```

#### 예제 16: 클러스터에서 Snakemake 실행
```bash
# SLURM 클러스터
snakemake --cluster "sbatch -p normal -t 2:00:00" --jobs 10

# SGE 클러스터
snakemake --cluster "qsub -cwd -V" --jobs 10
```

#### 예제 17: Snakemake 워크플로우 시각화
```bash
# DAG 생성
snakemake --dag | dot -Tpng > workflow.png

# 규칙 그래프
snakemake --rulegraph | dot -Tpng > rules.png
```

### Nextflow 예제

#### 예제 18: Nextflow 기본 실행
```bash
nextflow run nextflow_pipeline.nf
```

#### 예제 19: Nextflow 파라미터 지정
```bash
nextflow run nextflow_pipeline.nf \
    --config my_config.yml \
    --outdir results/nextflow_run
```

#### 예제 20: Nextflow 리포트 생성
```bash
nextflow run nextflow_pipeline.nf \
    -with-report report.html \
    -with-timeline timeline.html \
    -with-dag flowchart.html
```

#### 예제 21: Nextflow 재시작
```bash
# 실패한 작업만 재실행
nextflow run nextflow_pipeline.nf -resume
```

---

## HPC 클러스터에서 실행

### 예제 22: SLURM 배치 작업
```bash
#!/bin/bash
#SBATCH --job-name=rnaseq_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=rnaseq_%j.log

# R 모듈 로드 (클러스터에 따라 다름)
module load R/4.3.0

# 작업 디렉토리로 이동
cd $SLURM_SUBMIT_DIR

# 파이프라인 실행
Rscript run_pipeline.R --verbose
```

실행:
```bash
sbatch slurm_job.sh
```

### 예제 23: SGE 배치 작업
```bash
#!/bin/bash
#$ -N rnaseq_analysis
#$ -cwd
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=04:00:00
#$ -o rnaseq_$JOB_ID.log

# R 모듈 로드
module load R

# 파이프라인 실행
Rscript run_pipeline.R --verbose
```

실행:
```bash
qsub sge_job.sh
```

### 예제 24: PBS 배치 작업
```bash
#!/bin/bash
#PBS -N rnaseq_analysis
#PBS -l nodes=1:ppn=4
#PBS -l mem=16gb
#PBS -l walltime=04:00:00
#PBS -o rnaseq.log
#PBS -e rnaseq.err

cd $PBS_O_WORKDIR

module load R

Rscript run_pipeline.R --verbose
```

실행:
```bash
qsub pbs_job.sh
```

---

## 고급 사용 예제

### 예제 25: 메모리 제한 증가
```bash
# R 메모리 제한을 16GB로 설정
export R_MAX_VSIZE=16Gb
Rscript run_pipeline.R
```

### 예제 26: 로그 파일 저장
```bash
# 로그를 파일로 저장
Rscript run_pipeline.R --verbose 2>&1 | tee analysis.log
```

### 예제 27: 타임스탬프와 함께 로그 저장
```bash
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
Rscript run_pipeline.R --verbose 2>&1 | tee logs/run_${TIMESTAMP}.log
```

### 예제 28: 에러 발생 시 알림
```bash
#!/bin/bash
Rscript run_pipeline.R || {
    echo "Analysis failed!" | mail -s "RNA-Seq Pipeline Error" user@example.com
    exit 1
}
```

### 예제 29: 성공 시 후처리
```bash
#!/bin/bash
if Rscript run_pipeline.R; then
    echo "Analysis completed successfully"
    # 결과 압축
    tar -czf results_$(date +%Y%m%d).tar.gz output/
    # 원격 서버로 전송
    scp results_*.tar.gz user@server:/path/to/backup/
fi
```

---

## 디버깅 예제

### 예제 30: 단계별 디버깅
```bash
# 1단계: DE 분석
Rscript run_pipeline.R --step de --verbose --output debug/step1

# 결과 확인 후 다음 단계
Rscript run_pipeline.R --step plots --verbose --output debug/step1

# 계속...
Rscript run_pipeline.R --step enrichment --verbose --output debug/step1
Rscript run_pipeline.R --step go_plots --verbose --output debug/step1
```

### 예제 31: R 인터랙티브 디버깅
```bash
# R 대화형 모드로 스크립트 실행
R --vanilla

# R 콘솔에서
> source("run_pipeline.R")
> # 에러 발생 시 traceback() 사용
> traceback()
```

### 예제 32: 의존성 확인
```bash
# CLI 테스트 실행
./test_cli.sh

# Conda 환경 확인
conda env list
conda list -n rna-seq-de-go-analysis | grep -E "r-base|bioconductor"
```

---

## 성능 최적화 예제

### 예제 33: 병렬 처리 (Snakemake)
```bash
# 모든 가용 코어 사용
snakemake --cores all

# 특정 개수의 코어 사용
snakemake --cores 8
```

### 예제 34: 재시작 기능 활용
```bash
# Snakemake - 완료된 작업 건너뛰기
snakemake --cores 4 --rerun-incomplete

# Nextflow - 이전 실행 재개
nextflow run nextflow_pipeline.nf -resume
```

---

## 모니터링 예제

### 예제 35: 진행 상황 모니터링
```bash
# 백그라운드에서 실행
Rscript run_pipeline.R &

# 프로세스 확인
ps aux | grep run_pipeline

# 로그 실시간 확인
tail -f output/*/analysis.log
```

### 예제 36: 리소스 사용량 모니터링
```bash
# htop으로 실시간 모니터링
htop -p $(pgrep -f run_pipeline)

# 또는 시스템 모니터링
watch -n 5 'ps aux | grep Rscript | grep -v grep'
```

---

## 문제 해결 예제

### 예제 37: 권한 문제 해결
```bash
# 실행 권한 부여
chmod +x run_pipeline.R run_pipeline.sh test_cli.sh

# 출력 디렉토리 권한 확인
ls -ld output/
```

### 예제 38: 경로 문제 해결
```bash
# 절대 경로 사용
Rscript /full/path/to/run_pipeline.R \
    --config /full/path/to/config.yml \
    --output /full/path/to/results
```

### 예제 39: 환경 변수 설정
```bash
# R 경로 추가
export PATH="/usr/local/bin:$PATH"

# R 라이브러리 경로 추가
export R_LIBS_USER="$HOME/R/library"

# 실행
Rscript run_pipeline.R
```

---

## 요약

이 문서는 39가지 다양한 CLI 사용 예제를 제공합니다:
- 기본 실행: 예제 1-3
- 단계별 실행: 예제 4-7
- 커스텀 설정: 예제 8-10
- 배치 처리: 예제 11-13
- 워크플로우 도구: 예제 14-21
- HPC 클러스터: 예제 22-24
- 고급 사용: 예제 25-29
- 디버깅: 예제 30-32
- 성능 최적화: 예제 33-34
- 모니터링: 예제 35-36
- 문제 해결: 예제 37-39

더 많은 정보는 [CLI_GUIDE.md](CLI_GUIDE.md)와 [QUICKSTART.md](QUICKSTART.md)를 참조하세요.
