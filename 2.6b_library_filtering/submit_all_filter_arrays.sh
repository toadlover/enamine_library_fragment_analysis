#!/bin/bash
# Submit all 53,086 chunks as six arrays while keeping total array concurrency
# near 100 jobs: 6 arrays x 17 concurrent elements = at most 102 jobs.

set -Eeuo pipefail

SCRIPT_DIR="/pi/summer.thyme-umw/2.6b_library_filtering/enamine_library_fragment_analysis/2.6b_library_filtering"
cd "${SCRIPT_DIR}"
mkdir -p logs

bsub -J 'library_00[1-10000]%17' -env 'all,OFFSET=0'     < submit_filter.job
bsub -J 'library_01[1-10000]%17' -env 'all,OFFSET=10000' < submit_filter.job
bsub -J 'library_02[1-10000]%17' -env 'all,OFFSET=20000' < submit_filter.job
bsub -J 'library_03[1-10000]%17' -env 'all,OFFSET=30000' < submit_filter.job
bsub -J 'library_04[1-10000]%17' -env 'all,OFFSET=40000' < submit_filter.job
bsub -J 'library_05[1-3086]%17'  -env 'all,OFFSET=50000' < submit_filter.job
