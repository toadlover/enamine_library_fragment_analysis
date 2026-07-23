#!/usr/bin/env bash
set -euo pipefail

# LSF array worker for one Enamine chunk.
#
# Arguments:
#   1. First chunk represented by this array block
#   2. Output root
#   3. Path to filter_enamine_chunk_similarity_panfs_safe_v2.py
#   4. Similarity cutoff
#   5. Maximum ligands
#   6. Progress interval
#
# The actual chunk is:
#   BLOCK_START + LSB_JOBINDEX - 1

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 BLOCK_START OUTPUT_ROOT FILTER_SCRIPT SIMILARITY_CUTOFF MAX_LIGANDS PROGRESS_EVERY" >&2
    exit 2
fi

BLOCK_START="$1"
OUTPUT_ROOT="$2"
FILTER_SCRIPT="$3"
SIMILARITY_CUTOFF="$4"
MAX_LIGANDS="$5"
PROGRESS_EVERY="$6"

if [[ -z "${LSB_JOBINDEX:-}" ]]; then
    echo "ERROR: LSB_JOBINDEX is not set. This script must run as an LSF array element." >&2
    exit 2
fi

CHUNK=$((BLOCK_START + LSB_JOBINDEX - 1))
SUPERCHUNK=$((CHUNK / 100))
CHUNK_PADDED=$(printf "%05d" "${CHUNK}")

CHUNK_OUTPUT_DIR="${OUTPUT_ROOT}/${SUPERCHUNK}/${CHUNK_PADDED}"
mkdir -p "${CHUNK_OUTPUT_DIR}"

FILTERED_SDF="${CHUNK_OUTPUT_DIR}/chunk_${CHUNK_PADDED}_filtered.sdf"
REPORT_CSV="${CHUNK_OUTPUT_DIR}/chunk_${CHUNK_PADDED}_similarity_report.csv"
REPORT_ARCHIVE="${REPORT_CSV}.tar.gz"
SUMMARY_JSON="${CHUNK_OUTPUT_DIR}/chunk_${CHUNK_PADDED}_summary.json"
SUCCESS_MARKER="${CHUNK_OUTPUT_DIR}/chunk_${CHUNK_PADDED}.complete"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting chunk ${CHUNK_PADDED}"
echo "Superchunk: ${SUPERCHUNK}"
echo "Output directory: ${CHUNK_OUTPUT_DIR}"

# Permit safe re-runs. A completed chunk is skipped.
if [[ -f "${SUCCESS_MARKER}" && -f "${REPORT_ARCHIVE}" && -f "${SUMMARY_JSON}" ]]; then
    echo "Chunk ${CHUNK_PADDED} is already complete; skipping."
    exit 0
fi

# Remove stale products from a prior incomplete attempt.
rm -f \
    "${SUCCESS_MARKER}" \
    "${REPORT_ARCHIVE}" \
    "${FILTERED_SDF}" \
    "${REPORT_CSV}" \
    "${REPORT_CSV}.partial" \
    "${FILTERED_SDF}.partial" \
    "${SUMMARY_JSON}.partial"

python "${FILTER_SCRIPT}" \
    "${CHUNK_PADDED}" \
    "${CHUNK_OUTPUT_DIR}" \
    --similarity-cutoff "${SIMILARITY_CUTOFF}" \
    --max-ligands "${MAX_LIGANDS}" \
    --progress-every "${PROGRESS_EVERY}" \
    --quiet-rdkit

if [[ ! -f "${REPORT_CSV}" ]]; then
    echo "ERROR: expected report was not produced: ${REPORT_CSV}" >&2
    exit 1
fi

if [[ ! -f "${SUMMARY_JSON}" ]]; then
    echo "ERROR: expected summary was not produced: ${SUMMARY_JSON}" >&2
    exit 1
fi

# Create the tar.gz first. Delete the CSV only after tar succeeds.
tar -czf "${REPORT_ARCHIVE}.partial" \
    -C "${CHUNK_OUTPUT_DIR}" \
    "$(basename "${REPORT_CSV}")"

mv "${REPORT_ARCHIVE}.partial" "${REPORT_ARCHIVE}"
rm -f "${REPORT_CSV}"

# The filtered SDF is intentionally not retained.
rm -f "${FILTERED_SDF}"

# Mark the chunk complete only after all required post-processing succeeds.
{
    echo "chunk=${CHUNK_PADDED}"
    echo "superchunk=${SUPERCHUNK}"
    echo "completed=$(date --iso-8601=seconds)"
    echo "report_archive=${REPORT_ARCHIVE}"
    echo "summary_json=${SUMMARY_JSON}"
} > "${SUCCESS_MARKER}.partial"

mv "${SUCCESS_MARKER}.partial" "${SUCCESS_MARKER}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed chunk ${CHUNK_PADDED}"
echo "Retained:"
echo "  ${REPORT_ARCHIVE}"
echo "  ${SUMMARY_JSON}"
echo "  ${SUCCESS_MARKER}"
