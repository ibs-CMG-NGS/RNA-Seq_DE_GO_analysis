#!/bin/bash
# 파일 경로: run_pipeline.sh
# RNA-Seq 분석 파이프라인 실행 스크립트

# 스크립트가 위치한 디렉토리로 이동
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# run_pipeline.R 실행
Rscript run_pipeline.R "$@"
