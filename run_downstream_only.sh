#!/bin/bash
# 외부 DE 결과로 downstream 분석 실행 예시
# Usage: bash run_downstream_only.sh

set -e  # 에러 발생시 중단

echo "=========================================="
echo "외부 DE 결과로 Downstream 분석 시작"
echo "=========================================="
echo ""

# ============================================
# 1. 설정
# ============================================
# 수정 필요: 사용자의 파일 경로로 변경하세요
DE_RESULT_FILE="your_de_results.csv"          # 외부 DE 결과 파일
COMPARISON="H2O2_vs_Control"                  # 비교쌍 이름
OUTPUT_DIR="output/H2O2_Neuron"               # 출력 디렉토리
CONFIG_FILE="config_H2O2_Neuron.yml"          # Config 파일

# ============================================
# 2. 사전 확인
# ============================================
echo "Step 1: 파일 존재 확인..."

if [ ! -f "$DE_RESULT_FILE" ]; then
    echo "❌ 오류: DE 결과 파일을 찾을 수 없습니다: $DE_RESULT_FILE"
    echo "   파일 경로를 확인하고 스크립트 상단의 DE_RESULT_FILE을 수정하세요."
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ 오류: Config 파일을 찾을 수 없습니다: $CONFIG_FILE"
    exit 1
fi

echo "✓ DE 결과 파일: $DE_RESULT_FILE"
echo "✓ Config 파일: $CONFIG_FILE"
echo ""

# ============================================
# 3. DE 결과 파일 준비
# ============================================
echo "Step 2: DE 결과 파일 준비..."

# Python 스크립트가 있는지 확인
if [ -f "prepare_external_de_results.py" ]; then
    echo "자동 스크립트를 사용하여 파일 준비 중..."
    python prepare_external_de_results.py \
        --input "$DE_RESULT_FILE" \
        --comparison "$COMPARISON" \
        --output-dir "$OUTPUT_DIR" \
        --config "$CONFIG_FILE"
    
    if [ $? -ne 0 ]; then
        echo "❌ 오류: 파일 준비 실패"
        exit 1
    fi
else
    echo "수동으로 파일 배치 중..."
    # 디렉토리 생성
    mkdir -p "$OUTPUT_DIR/pairwise/$COMPARISON"
    
    # 파일 복사
    cp "$DE_RESULT_FILE" "$OUTPUT_DIR/pairwise/$COMPARISON/final_de_results.csv"
    cp "$CONFIG_FILE" "$OUTPUT_DIR/pairwise/$COMPARISON/config_used.yml"
    
    echo "✓ 파일 배치 완료: $OUTPUT_DIR/pairwise/$COMPARISON/"
fi

echo ""

# ============================================
# 4. Dry-run으로 실행 계획 확인
# ============================================
echo "Step 3: Snakemake 실행 계획 확인 (dry-run)..."
echo ""

snakemake --dryrun --printshellcmds \
    "$OUTPUT_DIR/pairwise/$COMPARISON/final_go_results.xlsx"

if [ $? -ne 0 ]; then
    echo ""
    echo "❌ Dry-run 실패. 설정을 확인하세요."
    exit 1
fi

echo ""
echo "=========================================="
echo "Dry-run 성공! 실제 분석을 시작하시겠습니까?"
echo "=========================================="
read -p "계속하려면 'y'를 입력하세요 (y/n): " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "분석 취소됨"
    exit 0
fi

# ============================================
# 5. 실제 분석 실행
# ============================================
echo ""
echo "Step 4: Downstream 분석 실행 중..."
echo "=========================================="
echo ""

# 코어 수 설정 (사용 가능한 코어 수에 맞게 조정)
CORES=4

# Snakemake 실행
snakemake --cores $CORES --use-conda \
    "$OUTPUT_DIR/pairwise/$COMPARISON/final_go_results.xlsx"

if [ $? -ne 0 ]; then
    echo ""
    echo "❌ 분석 실패"
    exit 1
fi

# ============================================
# 6. 결과 확인
# ============================================
echo ""
echo "=========================================="
echo "✅ 분석 완료!"
echo "=========================================="
echo ""
echo "생성된 파일:"
echo ""

# 주요 결과 파일 나열
RESULT_DIR="$OUTPUT_DIR/pairwise/$COMPARISON"

if [ -f "$RESULT_DIR/final_go_results.xlsx" ]; then
    echo "✓ GO Summary Table: $RESULT_DIR/final_go_results.xlsx"
fi

echo ""
echo "Enrichment CSV 파일:"
ls -1 "$RESULT_DIR"/go_enrichment_*.csv 2>/dev/null | head -5
echo "   ... (더 많은 파일)"

echo ""
echo "Visualization PNG 파일:"
ls -1 "$RESULT_DIR"/*.png 2>/dev/null | head -5
echo "   ... (더 많은 파일)"

echo ""
echo "전체 결과 디렉토리: $RESULT_DIR"
echo ""

# ============================================
# 7. 다음 단계 안내
# ============================================
echo "=========================================="
echo "다음 단계:"
echo "=========================================="
echo ""
echo "1. 결과 확인:"
echo "   cd $RESULT_DIR"
echo "   ls -lh"
echo ""
echo "2. Excel 파일 열기:"
echo "   open $RESULT_DIR/final_go_results.xlsx"
echo ""
echo "3. 추가 비교쌍 분석:"
echo "   - 스크립트 상단의 COMPARISON 변수 수정"
echo "   - 다시 실행: bash run_downstream_only.sh"
echo ""
echo "4. QC plots 생성 (선택사항):"
echo "   config에서 generate_pairwise_qc: true 설정 후"
echo "   snakemake --cores $CORES --use-conda"
echo ""
