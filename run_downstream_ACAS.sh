#!/bin/bash
# ACAS downstream 분석 실행 스크립트
# final_de_results.csv가 이미 있는 상태에서 enrichment부터 실행

set -e

echo "=========================================="
echo "ACAS Downstream 분석 시작"
echo "=========================================="
echo ""

# 환경 활성화
echo "Conda 환경 활성화 중..."
conda activate rna-seq-de-go-analysis

# Snakemake 실행
echo "Snakemake 실행 중..."
snakemake --use-conda --cores 4 \
    output/ACAS/pairwise/SVad_vs_Control/final_go_results.xlsx

echo ""
echo "=========================================="
echo "✅ 분석 완료!"
echo "=========================================="
echo ""
echo "결과 파일 위치:"
echo "  - GO enrichment: output/ACAS/pairwise/SVad_vs_Control/go_enrichment_*.csv"
echo "  - KEGG enrichment: output/ACAS/pairwise/SVad_vs_Control/kegg_enrichment_*.csv"
echo "  - GO summary table: output/ACAS/pairwise/SVad_vs_Control/final_go_results.xlsx"
echo ""
