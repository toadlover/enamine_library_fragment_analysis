#!/bin/bash
set -Eeuo pipefail

# Six arrays cover chunks 00000-53085.
#
# Because scratch usage, not RAM, is the main constraint, this defaults to only
# four running elements per array: 6 x 4 = at most 24 concurrent jobs.
#
# Override after measuring representative peak workspace sizes:
#   MAX_PER_ARRAY=8 ./submit_all_manifest_arrays.sh
#
# Safe concurrency estimate:
#   total_concurrency <= usable_scratch_GiB / peak_workspace_GiB_per_job
#
# Leave additional headroom for filesystem metadata, logs, transfer-failure
# outputs, and variation among chunks.

readonly MAX_PER_ARRAY="${MAX_PER_ARRAY:-4}"
readonly JOB_FILE="${JOB_FILE:-submit_filter_manifest.job}"

if ! [[ "${MAX_PER_ARRAY}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MAX_PER_ARRAY must be a positive integer" >&2
    exit 2
fi

mkdir -p logs

echo "Submitting six arrays with %${MAX_PER_ARRAY} each"
echo "Maximum simultaneous jobs: $((6 * MAX_PER_ARRAY))"

bsub -J "manifest_00[1-10000]%${MAX_PER_ARRAY}" -env "all,OFFSET=0"     < "${JOB_FILE}"
bsub -J "manifest_01[1-10000]%${MAX_PER_ARRAY}" -env "all,OFFSET=10000" < "${JOB_FILE}"
bsub -J "manifest_02[1-10000]%${MAX_PER_ARRAY}" -env "all,OFFSET=20000" < "${JOB_FILE}"
bsub -J "manifest_03[1-10000]%${MAX_PER_ARRAY}" -env "all,OFFSET=30000" < "${JOB_FILE}"
bsub -J "manifest_04[1-10000]%${MAX_PER_ARRAY}" -env "all,OFFSET=40000" < "${JOB_FILE}"
bsub -J "manifest_05[1-3086]%${MAX_PER_ARRAY}"  -env "all,OFFSET=50000" < "${JOB_FILE}"
