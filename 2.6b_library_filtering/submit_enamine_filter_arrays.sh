#!/usr/bin/env bash
set -euo pipefail

# Submit Enamine chunk filtering as one or more LSF job arrays.
#
# Defaults are suitable for a first test over chunks 00000-00199.
#
# Example test:
#   bash submit_enamine_filter_arrays.sh \
#       --start 0 \
#       --end 199 \
#       --output-root /pi/summer.thyme-umw/2.6b_library_filtering/output
#
# Full run:
#   bash submit_enamine_filter_arrays.sh \
#       --start 0 \
#       --end 53083 \
#       --output-root /pi/summer.thyme-umw/2.6b_library_filtering/output \
#       --block-size 1000 \
#       --max-concurrent 1000

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/run_enamine_filter_array_job.sh"

START_CHUNK=0
END_CHUNK=199
OUTPUT_ROOT=""
FILTER_SCRIPT="${SCRIPT_DIR}/filter_enamine_chunk_similarity_panfs_safe_v2.py"
BLOCK_SIZE=1000
MAX_CONCURRENT=1000
SIMILARITY_CUTOFF=0.75
MAX_LIGANDS=50000
PROGRESS_EVERY=1000

usage() {
    cat <<EOF
Usage: $0 --output-root PATH [options]

Required:
  --output-root PATH       Root for <superchunk>/<chunk>/ output directories.

Range:
  --start N                First global chunk. Default: 0
  --end N                  Last global chunk, inclusive. Default: 199

Submission:
  --block-size N           Array elements per bsub submission. Default: 1000
  --max-concurrent N       Maximum running elements per array. Default: 1000

Filtering:
  --filter-script PATH     Filtering Python script.
                           Default: ${FILTER_SCRIPT}
  --similarity-cutoff X    Default: 0.75
  --max-ligands N          Default: 50000
  --progress-every N       Default: 1000

LSF resources are fixed at:
  queue=long, runtime=1:00, cores=1, memory=4000 MB
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --start)
            START_CHUNK="$2"; shift 2 ;;
        --end)
            END_CHUNK="$2"; shift 2 ;;
        --output-root)
            OUTPUT_ROOT="$2"; shift 2 ;;
        --filter-script)
            FILTER_SCRIPT="$2"; shift 2 ;;
        --block-size)
            BLOCK_SIZE="$2"; shift 2 ;;
        --max-concurrent)
            MAX_CONCURRENT="$2"; shift 2 ;;
        --similarity-cutoff)
            SIMILARITY_CUTOFF="$2"; shift 2 ;;
        --max-ligands)
            MAX_LIGANDS="$2"; shift 2 ;;
        --progress-every)
            PROGRESS_EVERY="$2"; shift 2 ;;
        -h|--help)
            usage
            exit 0 ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage >&2
            exit 2 ;;
    esac
done

if [[ -z "${OUTPUT_ROOT}" ]]; then
    echo "ERROR: --output-root is required." >&2
    usage >&2
    exit 2
fi

for value in START_CHUNK END_CHUNK BLOCK_SIZE MAX_CONCURRENT MAX_LIGANDS PROGRESS_EVERY; do
    if ! [[ "${!value}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: ${value} must be a nonnegative integer." >&2
        exit 2
    fi
done

if (( START_CHUNK > END_CHUNK )); then
    echo "ERROR: --start cannot be greater than --end." >&2
    exit 2
fi

if (( END_CHUNK > 53083 )); then
    echo "ERROR: maximum valid chunk is 53083." >&2
    exit 2
fi

if (( BLOCK_SIZE < 1 || MAX_CONCURRENT < 1 )); then
    echo "ERROR: --block-size and --max-concurrent must be at least 1." >&2
    exit 2
fi

if [[ ! -f "${WORKER_SCRIPT}" ]]; then
    echo "ERROR: worker script not found: ${WORKER_SCRIPT}" >&2
    exit 2
fi

if [[ ! -f "${FILTER_SCRIPT}" ]]; then
    echo "ERROR: filter script not found: ${FILTER_SCRIPT}" >&2
    exit 2
fi

mkdir -p "${OUTPUT_ROOT}"
OUTPUT_ROOT="$(cd "${OUTPUT_ROOT}" && pwd)"
FILTER_SCRIPT="$(cd "$(dirname "${FILTER_SCRIPT}")" && pwd)/$(basename "${FILTER_SCRIPT}")"
WORKER_SCRIPT="$(cd "$(dirname "${WORKER_SCRIPT}")" && pwd)/$(basename "${WORKER_SCRIPT}")"

LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

TOTAL=$((END_CHUNK - START_CHUNK + 1))
echo "Submitting ${TOTAL} chunks: ${START_CHUNK}-${END_CHUNK}"
echo "Output root: ${OUTPUT_ROOT}"
echo "Block size: ${BLOCK_SIZE}"
echo "Maximum concurrent elements per array: ${MAX_CONCURRENT}"

BLOCK_START="${START_CHUNK}"
SUBMISSION_COUNT=0

while (( BLOCK_START <= END_CHUNK )); do
    BLOCK_END=$((BLOCK_START + BLOCK_SIZE - 1))
    if (( BLOCK_END > END_CHUNK )); then
        BLOCK_END="${END_CHUNK}"
    fi

    ELEMENT_COUNT=$((BLOCK_END - BLOCK_START + 1))
    ARRAY_NAME="enamine_filter_${BLOCK_START}_${BLOCK_END}[1-${ELEMENT_COUNT}]%${MAX_CONCURRENT}"

    echo "Submitting block ${BLOCK_START}-${BLOCK_END} (${ELEMENT_COUNT} elements)"

    bsub \
        -J "${ARRAY_NAME}" \
        -q long \
        -W 1:00 \
        -n 1 \
        -M 4000 \
        -R "rusage[mem=4000] span[hosts=1]" \
        -o "${LOG_DIR}/%J.%I.out" \
        -e "${LOG_DIR}/%J.%I.err" \
        /bin/bash "${WORKER_SCRIPT}" \
            "${BLOCK_START}" \
            "${OUTPUT_ROOT}" \
            "${FILTER_SCRIPT}" \
            "${SIMILARITY_CUTOFF}" \
            "${MAX_LIGANDS}" \
            "${PROGRESS_EVERY}"

    SUBMISSION_COUNT=$((SUBMISSION_COUNT + 1))
    BLOCK_START=$((BLOCK_END + 1))
done

echo "Submitted ${SUBMISSION_COUNT} LSF array block(s)."
