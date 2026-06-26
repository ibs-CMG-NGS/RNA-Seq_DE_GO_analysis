#!/bin/bash
# run_batch.sh — DE-GO 파이프라인 일괄 실행
#
# 사용법:
#   bash run_batch.sh [옵션]
#
# 옵션:
#   -f, --force            조건 변경 후 완료된 작업도 강제 재실행 (--forceall)
#   -n, --dry-run          실제 실행 없이 계획만 출력
#   -j, --cores N          프로젝트당 사용할 코어 수 (기본값: 4)
#   -i, --include A,B,...  지정한 프로젝트만 실행 (쉼표 구분, 부분 일치)
#   -e, --exclude A,B,...  지정한 프로젝트 제외 (쉼표 구분, 부분 일치)
#   -l, --list             실행 대상 프로젝트 목록만 출력
#   -h, --help             도움말
#
# 예시:
#   bash run_batch.sh                              # 모든 프로젝트 순차 실행
#   bash run_batch.sh --force                      # 전체 강제 재실행 (조건 변경 시)
#   bash run_batch.sh --include mouse-h2o2,kkj     # 이름에 해당 문자열 포함된 것만 실행
#   bash run_batch.sh --exclude 2026-hiy           # 해당 프로젝트 제외
#   bash run_batch.sh --dry-run --include mouse    # 실행 계획만 확인
#   bash run_batch.sh --list                       # 실행 가능한 프로젝트 목록 확인

set -uo pipefail

# ============================================================
# 경로 설정
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIGS_DIR="$SCRIPT_DIR/configs"
BATCH_LOG_DIR="$SCRIPT_DIR/logs/batch"
CONDA_ENV="rna-seq-de-go-analysis"

# ============================================================
# 옵션 파싱
# ============================================================
FORCEALL=false
DRYRUN=false
LIST_ONLY=false
CORES=4
INCLUDE_PATTERNS=()
EXCLUDE_PATTERNS=()

usage() {
    sed -n '2,26p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--force)    FORCEALL=true ;;
        -n|--dry-run)  DRYRUN=true ;;
        -l|--list)     LIST_ONLY=true ;;
        -j|--cores)    CORES="$2"; shift ;;
        -i|--include)  IFS=',' read -ra INCLUDE_PATTERNS <<< "$2"; shift ;;
        -e|--exclude)  IFS=',' read -ra EXCLUDE_PATTERNS <<< "$2"; shift ;;
        -h|--help)     usage ;;
        *) echo "알 수 없는 옵션: $1"; echo ""; usage ;;
    esac
    shift
done

# ============================================================
# 실행 대상 수집
# ============================================================
declare -a TARGETS=()   # "프로젝트명:config경로" 형식

for config_file in "$CONFIGS_DIR"/config_*.yml; do
    proj=$(basename "$config_file" .yml | sed 's/config_//')

    # --include 필터 (부분 일치, 하나라도 매칭되면 포함)
    if [[ ${#INCLUDE_PATTERNS[@]} -gt 0 ]]; then
        matched=false
        for pat in "${INCLUDE_PATTERNS[@]}"; do
            [[ "$proj" == *"$pat"* ]] && matched=true && break
        done
        $matched || continue
    fi

    # --exclude 필터 (부분 일치, 하나라도 매칭되면 제외)
    skip=false
    for pat in "${EXCLUDE_PATTERNS[@]}"; do
        [[ "$proj" == *"$pat"* ]] && skip=true && break
    done
    $skip && continue

    # 데이터 파일 존재 확인
    count_path=$(grep "^count_data_path:" "$config_file" 2>/dev/null | awk '{print $2}')
    if [[ -z "$count_path" || ! -f "$SCRIPT_DIR/$count_path" ]]; then
        echo "[SKIP] $proj — 데이터 파일 없음: ${count_path:-'count_data_path 미설정'}"
        continue
    fi

    TARGETS+=("$proj:$config_file")
done

if [[ ${#TARGETS[@]} -eq 0 ]]; then
    echo "실행할 프로젝트가 없습니다."
    exit 1
fi

# --list 모드: 목록만 출력하고 종료
if $LIST_ONLY; then
    echo "실행 대상 프로젝트 (${#TARGETS[@]}개):"
    for entry in "${TARGETS[@]}"; do
        proj="${entry%%:*}"
        config="${entry##*:}"
        outdir=$(grep "^output_dir:" "$config" | awk '{print $2}')
        echo "  - $proj  →  $outdir"
    done
    exit 0
fi

# ============================================================
# 배치 실행 준비
# ============================================================
mkdir -p "$BATCH_LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
BATCH_LOG="$BATCH_LOG_DIR/batch_${TIMESTAMP}.log"

SNAKE_OPTS="--cores $CORES --use-conda"
$FORCEALL && SNAKE_OPTS="$SNAKE_OPTS --forceall"
$DRYRUN   && SNAKE_OPTS="$SNAKE_OPTS --dryrun"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  DE-GO Batch Pipeline"
echo "  시작    : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  프로젝트: ${#TARGETS[@]}개"
echo "  코어    : $CORES"
$FORCEALL && echo "  모드    : 강제 재실행 (--forceall)"
$DRYRUN   && echo "  모드    : Dry-run (실제 실행 없음)"
echo "  배치 로그: $BATCH_LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

{
    echo "=== DE-GO Batch Run ==="
    echo "Start   : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Projects: ${#TARGETS[@]}"
    echo "Options : $SNAKE_OPTS"
    echo ""
} > "$BATCH_LOG"

declare -a SUCCESS=()
declare -a FAILED=()
TOTAL=${#TARGETS[@]}
IDX=0

# ============================================================
# 프로젝트별 순차 실행
# ============================================================
for entry in "${TARGETS[@]}"; do
    proj="${entry%%:*}"
    config_file="${entry##*:}"
    IDX=$((IDX + 1))
    proj_log="$BATCH_LOG_DIR/${proj}_${TIMESTAMP}.log"

    echo "▶ [$IDX/$TOTAL] $proj"
    echo "  Config : $(basename "$config_file")"

    start_ts=$(date +%s)

    set +e
    conda run -n "$CONDA_ENV" \
        snakemake --configfile "$config_file" $SNAKE_OPTS \
        > "$proj_log" 2>&1
    exit_code=$?
    set -e

    elapsed=$(( $(date +%s) - start_ts ))
    elapsed_fmt=$(printf '%02d:%02d' $((elapsed/60)) $((elapsed%60)))

    if [[ $exit_code -eq 0 ]]; then
        echo "  ✓ 완료 (${elapsed_fmt})"
        SUCCESS+=("$proj")
        echo "[$IDX/$TOTAL] $proj  →  SUCCESS  (${elapsed_fmt})  log: $proj_log" >> "$BATCH_LOG"
    else
        echo "  ✗ 실패 (${elapsed_fmt})  →  로그: $proj_log"
        FAILED+=("$proj")
        echo "[$IDX/$TOTAL] $proj  →  FAILED   (${elapsed_fmt})  log: $proj_log" >> "$BATCH_LOG"
        # 실패한 마지막 10줄을 배치 로그에 첨부
        echo "    --- 오류 마지막 10줄 ---" >> "$BATCH_LOG"
        tail -10 "$proj_log" | sed 's/^/    /' >> "$BATCH_LOG"
        echo "" >> "$BATCH_LOG"
    fi
    echo ""
done

# ============================================================
# 최종 요약
# ============================================================
{
    echo ""
    echo "=== 결과 요약 ==="
    echo "완료: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "성공: ${#SUCCESS[@]} / $TOTAL"
    [[ ${#FAILED[@]} -gt 0 ]] && echo "실패: ${#FAILED[@]} / $TOTAL"
} >> "$BATCH_LOG"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  완료: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
printf "  ✓ 성공 (%d)\n" "${#SUCCESS[@]}"
for p in "${SUCCESS[@]}"; do echo "      $p"; done

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo ""
    printf "  ✗ 실패 (%d)\n" "${#FAILED[@]}"
    for p in "${FAILED[@]}"; do echo "      $p"; done
fi

echo ""
echo "  배치 로그: $BATCH_LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

[[ ${#FAILED[@]} -gt 0 ]] && exit 1
exit 0
