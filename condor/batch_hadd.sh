#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 OUT_FILE \"IN_FILES\" BATCH_SIZE NJOBS [-z|--zombie-check]" >&2
    echo "  Note: IN_FILES must be quoted to prevent shell expansion" >&2
    exit 1
}

if [[ $# -lt 4 || $# -gt 5 ]]; then
    usage
fi

OUT_FILE="$1"
IN_FILES="$2"
BATCH_SIZE="${3:-10}"
NJOBS="${4:-4}"

ZOMBIE_CHECK=false
if [[ $# -eq 5 ]]; then
    if [[ "$5" == "-z" || "$5" == "--zombie-check" ]]; then
        ZOMBIE_CHECK=true
    else
        usage
    fi
fi

echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"

source "$(dirname "${BASH_SOURCE[0]}")/draw_bar.sh"
pick_bar_color

MY_TMPDIR="$(dirname "$OUT_FILE")/hadd_tmp_$$"
LOG_FILE="$(dirname "$OUT_FILE")/hadd_$$.log"
eos mkdir -p "$MY_TMPDIR"

cleanup() {
    if (( $? == 0 )); then
        eos rm -r "$MY_TMPDIR"
        rm -f "$LOG_FILE"
    else
        if [[ -s "$LOG_FILE" ]]; then
            echo "Failure — log written to ${LOG_FILE}" >&2
        fi
        echo "Leaving ${MY_TMPDIR} intact for inspection" >&2
    fi
}
trap cleanup EXIT

# format_elapsed SECONDS → "Xh Ym Zs" (omits leading zero components)
format_elapsed() {
    local secs="$1"
    local h=$(( secs / 3600 ))
    local m=$(( (secs % 3600) / 60 ))
    local s=$(( secs % 60 ))
    if (( h > 0 )); then
        printf "%dh %dm %ds" "$h" "$m" "$s"
    elif (( m > 0 )); then
        printf "%dm %ds" "$m" "$s"
    else
        printf "%ds" "$s"
    fi
}

# normalise IN_FILES to a directory + filename pattern
IN_DIR="$(dirname "$IN_FILES")"
IN_PAT="$(basename "$IN_FILES")"

if [[ -d "$IN_FILES" ]]; then
    IN_DIR="$IN_FILES"
    IN_PAT="*.root"
elif [[ "$IN_PAT" == *root && "$IN_PAT" != *.root ]]; then
    IN_PAT="${IN_PAT/%root/.root}"
fi

mapfile -t FILES < <(find "$IN_DIR" -type f -name "$IN_PAT" | sort)

if (( ${#FILES[@]} == 0 )); then
    echo "No files matched: $IN_FILES" >&2
    echo "Pattern must be a valid path with files and wildcard(s), e.g. /eos/user/username/mydir/type1_*.root" >&2
    exit 1
fi

TOTAL_INPUT=${#FILES[@]}
echo "Found ${TOTAL_INPUT} files"

zombie_files=()

if $ZOMBIE_CHECK; then
    CHECK_JOBS=$(( NJOBS > 4 ? NJOBS : 4 ))
    echo "Checking for zombies (${CHECK_JOBS} parallel)..."
    ZOMBIE_START=$(date +%s)

    valid_files=()
    pids=()
    tmpdir_check=$(mktemp -d)
    checked=0
    declare -A collected

    start_bar_timer
    draw_bar "$BAR_COLOR" "Zombie check:" 0 "$TOTAL_INPUT"

    collect_check_results() {
        for p in "${pids[@]}"; do
            if [[ -z "${collected[$p]+x}" ]] && [[ -f "${tmpdir_check}/${p}" ]]; then
                local r
                r=$(< "${tmpdir_check}/${p}")
                if [[ "$r" == valid:* ]]; then
                    valid_files+=("${r#valid:}")
                else
                    zombie_files+=("${r#zombie:}")
                fi
                rm -f "${tmpdir_check}/${p}"
                collected[$p]=1
                checked=$(( checked + 1 ))
                draw_bar "$BAR_COLOR" "Zombie check:" "$checked" "$TOTAL_INPUT"
            fi
        done
    }

    for f in "${FILES[@]}"; do
        while (( $(jobs -rp | wc -l) >= CHECK_JOBS )); do
            collect_check_results
            sleep 0.1
        done

        (
            if rootls "$f" &>/dev/null; then
                echo "valid:$f" > "${tmpdir_check}/$BASHPID"
            else
                echo "zombie:$f" > "${tmpdir_check}/$BASHPID"
            fi
        ) &
        pids+=($!)

        collect_check_results
    done

    for pid in "${pids[@]}"; do
        wait "$pid"
        collect_check_results
    done

    rm -rf "$tmpdir_check"
    draw_bar "$BAR_COLOR" "Zombie check:" "$TOTAL_INPUT" "$TOTAL_INPUT"
    printf "\n\n"

    ZOMBIE_END=$(date +%s)
    ZOMBIE_ELAPSED=$(( ZOMBIE_END - ZOMBIE_START ))

    if (( ${#zombie_files[@]} > 0 )); then
        for zf in "${zombie_files[@]}"; do
            echo "  [ZOMBIE] $zf"
        done
        echo "Warning: skipping ${#zombie_files[@]} zombie/corrupt file(s)" >&2
    fi

    echo "Zombie check complete in $(format_elapsed $ZOMBIE_ELAPSED)"

    if (( ${#valid_files[@]} == 0 )); then
        echo "No valid ROOT files remain after zombie check." >&2
        exit 1
    fi

    FILES=("${valid_files[@]}")
    echo "Proceeding with ${#FILES[@]} valid files"
fi

TOTAL_FILES=${#FILES[@]}
echo "Temporary directory: ${MY_TMPDIR}"

level=0
current_files=("${FILES[@]}")
declare -a level_summary

while (( ${#current_files[@]} > 1 )); do
    n_batches=$(( (${#current_files[@]} + BATCH_SIZE - 1) / BATCH_SIZE ))
    echo "Level ${level}: ${#current_files[@]} files → ${n_batches} batches"
    LEVEL_START=$(date +%s)

    next_files=()
    batch=0
    completed=0
    pids=()
    tmpdir_hadd=$(mktemp -d)
    declare -A collected_hadd

    start_bar_timer
    draw_bar "$BAR_COLOR" "Level ${level}:" 0 "$n_batches"

    collect_hadd_results() {
        for p in "${pids[@]}"; do
            if [[ -z "${collected_hadd[$p]+x}" ]] && [[ -f "${tmpdir_hadd}/${p}" ]]; then
                if [[ $(< "${tmpdir_hadd}/${p}") != "0" ]]; then
                    printf "\n"
                    echo "A hadd batch failed at level ${level} — check ${LOG_FILE}" >&2
                    exit 1
                fi
                rm -f "${tmpdir_hadd}/${p}"
                collected_hadd[$p]=1
                completed=$(( completed + 1 ))
                draw_bar "$BAR_COLOR" "Level ${level}:" "$completed" "$n_batches"
            fi
        done
    }

    for ((i=0; i<${#current_files[@]}; i+=BATCH_SIZE)); do
        while (( $(jobs -rp | wc -l) >= NJOBS )); do
            collect_hadd_results
            sleep 0.5
        done

        chunk=("${current_files[@]:i:BATCH_SIZE}")
        partial="${MY_TMPDIR}/level${level}_batch${batch}.root"
        next_files+=("$partial")

        (
            hadd -k "$partial" "${chunk[@]}" >>"$LOG_FILE" 2>&1
            echo $? > "${tmpdir_hadd}/$BASHPID"
        ) &
        pids+=($!)
        batch=$(( batch + 1 ))
    done

    for pid in "${pids[@]}"; do
        wait "$pid"
        collect_hadd_results
    done

    draw_bar "$BAR_COLOR" "Level ${level}:" "$n_batches" "$n_batches"
    printf "\n\n"

    LEVEL_END=$(date +%s)
    LEVEL_ELAPSED=$(( LEVEL_END - LEVEL_START ))
    echo "Level ${level} complete in $(format_elapsed $LEVEL_ELAPSED)"

    rm -rf "$tmpdir_hadd"
    unset collected_hadd
    level_summary+=("level ${level}: ${batch} batches in $(format_elapsed $LEVEL_ELAPSED)")
    current_files=("${next_files[@]}")
    level=$(( level + 1 ))
done

echo "Writing final output: ${OUT_FILE}"
mv "${current_files[0]}" "$OUT_FILE"

echo ""
echo "=== Merge Summary ==="
echo "Input files found:  ${TOTAL_INPUT}"
echo "Zombies skipped:    ${#zombie_files[@]}"
echo "Files merged:       ${TOTAL_FILES}"
echo "Max batch size:     ${BATCH_SIZE}"
echo "Total levels:       ${level}"
for entry in "${level_summary[@]}"; do
    echo "  ${entry}"
done
echo "Output: ${OUT_FILE}"
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "====================="