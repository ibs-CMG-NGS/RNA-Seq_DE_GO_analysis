# 외부 DE 결과로 Downstream 분석 실행 - PowerShell 버전
# Usage: .\run_downstream_only.ps1

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "외부 DE 결과로 Downstream 분석 시작" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""

# ============================================
# 1. 설정
# ============================================
# 수정 필요: 사용자의 파일 경로로 변경하세요
$DE_RESULT_FILE = "your_de_results.csv"          # 외부 DE 결과 파일
$COMPARISON = "H2O2_vs_Control"                  # 비교쌍 이름
$OUTPUT_DIR = "output/H2O2_Neuron"               # 출력 디렉토리
$CONFIG_FILE = "config_H2O2_Neuron.yml"          # Config 파일

# ============================================
# 2. 사전 확인
# ============================================
Write-Host "Step 1: 파일 존재 확인..." -ForegroundColor Yellow

if (-Not (Test-Path $DE_RESULT_FILE)) {
    Write-Host "❌ 오류: DE 결과 파일을 찾을 수 없습니다: $DE_RESULT_FILE" -ForegroundColor Red
    Write-Host "   파일 경로를 확인하고 스크립트 상단의 DE_RESULT_FILE을 수정하세요." -ForegroundColor Red
    exit 1
}

if (-Not (Test-Path $CONFIG_FILE)) {
    Write-Host "❌ 오류: Config 파일을 찾을 수 없습니다: $CONFIG_FILE" -ForegroundColor Red
    exit 1
}

Write-Host "✓ DE 결과 파일: $DE_RESULT_FILE" -ForegroundColor Green
Write-Host "✓ Config 파일: $CONFIG_FILE" -ForegroundColor Green
Write-Host ""

# ============================================
# 3. DE 결과 파일 준비
# ============================================
Write-Host "Step 2: DE 결과 파일 준비..." -ForegroundColor Yellow

# Python 스크립트가 있는지 확인
if (Test-Path "prepare_external_de_results.py") {
    Write-Host "자동 스크립트를 사용하여 파일 준비 중..." -ForegroundColor Cyan
    
    python prepare_external_de_results.py `
        --input $DE_RESULT_FILE `
        --comparison $COMPARISON `
        --output-dir $OUTPUT_DIR `
        --config $CONFIG_FILE
    
    if ($LASTEXITCODE -ne 0) {
        Write-Host "❌ 오류: 파일 준비 실패" -ForegroundColor Red
        exit 1
    }
} else {
    Write-Host "수동으로 파일 배치 중..." -ForegroundColor Cyan
    
    # 디렉토리 생성
    $targetDir = "$OUTPUT_DIR/pairwise/$COMPARISON"
    New-Item -ItemType Directory -Force -Path $targetDir | Out-Null
    
    # 파일 복사
    Copy-Item $DE_RESULT_FILE "$targetDir/final_de_results.csv"
    Copy-Item $CONFIG_FILE "$targetDir/config_used.yml"
    
    Write-Host "✓ 파일 배치 완료: $targetDir/" -ForegroundColor Green
}

Write-Host ""

# ============================================
# 4. Dry-run으로 실행 계획 확인
# ============================================
Write-Host "Step 3: Snakemake 실행 계획 확인 (dry-run)..." -ForegroundColor Yellow
Write-Host ""

$targetFile = "$OUTPUT_DIR/pairwise/$COMPARISON/final_go_results.xlsx"
snakemake --dryrun --printshellcmds $targetFile

if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "❌ Dry-run 실패. 설정을 확인하세요." -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Dry-run 성공! 실제 분석을 시작하시겠습니까?" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan

$response = Read-Host "계속하려면 'y'를 입력하세요 (y/n)"

if ($response -ne 'y' -and $response -ne 'Y') {
    Write-Host "분석 취소됨" -ForegroundColor Yellow
    exit 0
}

# ============================================
# 5. 실제 분석 실행
# ============================================
Write-Host ""
Write-Host "Step 4: Downstream 분석 실행 중..." -ForegroundColor Yellow
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""

# 코어 수 설정 (사용 가능한 코어 수에 맞게 조정)
$CORES = 4

# Snakemake 실행
snakemake --cores $CORES --use-conda $targetFile

if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "❌ 분석 실패" -ForegroundColor Red
    exit 1
}

# ============================================
# 6. 결과 확인
# ============================================
Write-Host ""
Write-Host "==========================================" -ForegroundColor Green
Write-Host "✅ 분석 완료!" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host ""
Write-Host "생성된 파일:" -ForegroundColor Yellow
Write-Host ""

# 주요 결과 파일 나열
$RESULT_DIR = "$OUTPUT_DIR/pairwise/$COMPARISON"

if (Test-Path "$RESULT_DIR/final_go_results.xlsx") {
    Write-Host "✓ GO Summary Table: $RESULT_DIR/final_go_results.xlsx" -ForegroundColor Green
}

Write-Host ""
Write-Host "Enrichment CSV 파일:" -ForegroundColor Yellow
Get-ChildItem "$RESULT_DIR/go_enrichment_*.csv" -ErrorAction SilentlyContinue | Select-Object -First 5 | ForEach-Object { Write-Host "  $_" }
Write-Host "   ... (더 많은 파일)"

Write-Host ""
Write-Host "Visualization PNG 파일:" -ForegroundColor Yellow
Get-ChildItem "$RESULT_DIR/*.png" -ErrorAction SilentlyContinue | Select-Object -First 5 | ForEach-Object { Write-Host "  $_" }
Write-Host "   ... (더 많은 파일)"

Write-Host ""
Write-Host "전체 결과 디렉토리: $RESULT_DIR" -ForegroundColor Cyan
Write-Host ""

# ============================================
# 7. 다음 단계 안내
# ============================================
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "다음 단계:" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "1. 결과 확인:"
Write-Host "   cd $RESULT_DIR"
Write-Host "   ls"
Write-Host ""
Write-Host "2. Excel 파일 열기:"
Write-Host "   start $RESULT_DIR/final_go_results.xlsx"
Write-Host ""
Write-Host "3. 추가 비교쌍 분석:"
Write-Host "   - 스크립트 상단의 `$COMPARISON 변수 수정"
Write-Host "   - 다시 실행: .\run_downstream_only.ps1"
Write-Host ""
Write-Host "4. QC plots 생성 (선택사항):"
Write-Host "   config에서 generate_pairwise_qc: true 설정 후"
Write-Host "   snakemake --cores $CORES --use-conda"
Write-Host ""
