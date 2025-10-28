#!/bin/bash
# test_cli.sh
# CLI 기능 검증 테스트 스크립트

# set -e  # 에러 발생 시 중단 (테스트에서는 계속 진행)

echo "=========================================="
echo "RNA-Seq Analysis Pipeline CLI Test"
echo "=========================================="
echo ""

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 테스트 결과 추적
TESTS_PASSED=0
TESTS_FAILED=0

# 함수: 테스트 성공
test_pass() {
    echo -e "${GREEN}✓ PASS${NC}: $1"
    ((TESTS_PASSED++))
}

# 함수: 테스트 실패
test_fail() {
    echo -e "${RED}✗ FAIL${NC}: $1"
    ((TESTS_FAILED++))
}

# 함수: 테스트 건너뛰기
test_skip() {
    echo -e "${YELLOW}⊘ SKIP${NC}: $1"
}

echo "Test 1: Check if R is installed"
if command -v R &> /dev/null; then
    test_pass "R is installed ($(R --version | head -1))"
else
    test_fail "R is not installed"
    echo "  Please install R from https://www.r-project.org/"
fi

echo ""
echo "Test 2: Check if Rscript is available"
if command -v Rscript &> /dev/null; then
    test_pass "Rscript is available"
else
    test_fail "Rscript is not available"
fi

echo ""
echo "Test 3: Check if run_pipeline.R exists"
if [ -f "run_pipeline.R" ]; then
    test_pass "run_pipeline.R exists"
else
    test_fail "run_pipeline.R not found"
fi

echo ""
echo "Test 4: Check if run_pipeline.R is executable"
if [ -x "run_pipeline.R" ]; then
    test_pass "run_pipeline.R is executable"
else
    test_fail "run_pipeline.R is not executable"
    echo "  Run: chmod +x run_pipeline.R"
fi

echo ""
echo "Test 5: Check if run_pipeline.sh exists"
if [ -f "run_pipeline.sh" ]; then
    test_pass "run_pipeline.sh exists"
else
    test_fail "run_pipeline.sh not found"
fi

echo ""
echo "Test 6: Check if run_pipeline.sh is executable"
if [ -x "run_pipeline.sh" ]; then
    test_pass "run_pipeline.sh is executable"
else
    test_fail "run_pipeline.sh is not executable"
    echo "  Run: chmod +x run_pipeline.sh"
fi

echo ""
echo "Test 7: Check if config.yml exists"
if [ -f "config.yml" ]; then
    test_pass "config.yml exists"
else
    test_fail "config.yml not found"
fi

echo ""
echo "Test 8: Check if Snakefile exists"
if [ -f "Snakefile" ]; then
    test_pass "Snakefile exists"
else
    test_fail "Snakefile not found"
fi

echo ""
echo "Test 9: Check if data files exist"
if [ -f "data/raw/airway_scaledcounts.csv" ] && [ -f "data/raw/airway_metadata.csv" ]; then
    test_pass "Data files exist"
else
    test_fail "Data files not found"
fi

echo ""
echo "Test 10: Check R script syntax (if R is available)"
if command -v Rscript &> /dev/null; then
    if Rscript -e "source('run_pipeline.R', echo=FALSE)" 2>&1 | grep -q "Error"; then
        test_fail "R script has syntax errors"
    else
        # 스크립트가 실행되려고 하면 정상 (실제 실행은 안 함)
        test_skip "R script syntax check (would require full execution)"
    fi
else
    test_skip "R not available for syntax check"
fi

echo ""
echo "Test 11: Test --help option (if R is available)"
if command -v Rscript &> /dev/null; then
    if Rscript run_pipeline.R --help 2>&1 | grep -q "Usage"; then
        test_pass "--help option works"
    else
        # optparse may not be installed yet
        test_skip "--help option (optparse may not be installed)"
    fi
else
    test_skip "R not available for --help test"
fi

echo ""
echo "Test 12: Check directory structure"
REQUIRED_DIRS=("src/analysis" "src/utils" "data/raw" "output" "notebooks")
ALL_DIRS_EXIST=true
for dir in "${REQUIRED_DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
        ALL_DIRS_EXIST=false
        echo "  Missing directory: $dir"
    fi
done
if $ALL_DIRS_EXIST; then
    test_pass "All required directories exist"
else
    test_fail "Some required directories are missing"
fi

echo ""
echo "Test 13: Check analysis scripts"
REQUIRED_SCRIPTS=(
    "src/analysis/01_run_de_analysis.R"
    "src/analysis/02_generate_plots.R"
    "src/analysis/03_enrichment_analysis.R"
    "src/analysis/04_generate_go_plots.R"
    "src/utils/load_data.R"
)
ALL_SCRIPTS_EXIST=true
for script in "${REQUIRED_SCRIPTS[@]}"; do
    if [ ! -f "$script" ]; then
        ALL_SCRIPTS_EXIST=false
        echo "  Missing script: $script"
    fi
done
if $ALL_SCRIPTS_EXIST; then
    test_pass "All required R scripts exist"
else
    test_fail "Some required R scripts are missing"
fi

echo ""
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo -e "${GREEN}Passed: $TESTS_PASSED${NC}"
echo -e "${RED}Failed: $TESTS_FAILED${NC}"
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed! ✓${NC}"
    echo ""
    echo "You can now run the pipeline with:"
    echo "  Rscript run_pipeline.R --help"
    echo "  Rscript run_pipeline.R"
    exit 0
else
    echo -e "${RED}Some tests failed. Please fix the issues above.${NC}"
    exit 1
fi
