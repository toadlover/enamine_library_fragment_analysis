#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# VERSION 2: EXPLICIT PARTIAL-NODE RESOURCE REQUESTS
#
# QOS constraints observed:
#   MaxJobsPU = 4
#   MaxTRESPU = cpu=64
#
# This launcher requests 16 CPUs per job, allowing up to four running jobs:
#   4 jobs x 16 CPUs = 64 CPUs
#
# Each job starts 8 independent extraction workers internally:
#   4 jobs x 8 workers = up to 32 simultaneous extraction batches
###############################################################################

INPUT_ROOT="/work/umassmed/thymelab_umassmed/compressed_unmodified_64b_library"
FEATURE_ROOT="/work/umassmed/thymelab_umassmed/65b_library_extracted_features/features"

NODE_PIPELINE="/work/umassmed/thymelab_umassmed/enamine_library_fragment_analysis/65b_fragment_processing/enamine_recipe_feature_pipeline_with_query_slurm_node_v2.py"

ORIGINAL_PIPELINE="/work/umassmed/thymelab_umassmed/enamine_library_fragment_analysis/65b_fragment_processing/enamine_recipe_feature_pipeline_with_query.py"

PYTHON="python"

###############################################################################
# SLURM SETTINGS
###############################################################################

PARTITION="cpu"
WALLTIME="24:00:00"

# Explicit allocation. Do not use --exclusive or --mem=0.
CPUS_PER_TASK=16
SLURM_MEMORY="128G"

# Independent extraction subprocesses inside each Slurm job.
LOCAL_WORKERS=8

###############################################################################
# PIPELINE SETTINGS
###############################################################################

BATCH_SIZE=3000000
PATTERN="*.cxsmiles.bz2"

# 1 = submit all jobs now; 0 = generate .slurm files only.
SUBMIT_JOBS=1

LOG_DIR="${FEATURE_ROOT}/logs"
JOB_SCRIPT_DIR="${FEATURE_ROOT}/full_node_slurm_scripts_v2"

mkdir -p "${LOG_DIR}" "${JOB_SCRIPT_DIR}"

job_count=0

for library_type in M S; do
    type_input_dir="${INPUT_ROOT}/${library_type}"
    [[ -d "${type_input_dir}" ]] || continue

    for h_dir in "${type_input_dir}"/H*; do
        [[ -d "${h_dir}" ]] || continue

        h_name="$(basename "${h_dir}")"
        output_root="${FEATURE_ROOT}/${library_type}/${h_name}"
        job_name="node_${library_type}_${h_name}"
        slurm_file="${JOB_SCRIPT_DIR}/${job_name}.slurm"

        mkdir -p "${output_root}"

        cat > "${slurm_file}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${job_name}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --mem=${SLURM_MEMORY}
#SBATCH --time=${WALLTIME}
#SBATCH --output=${LOG_DIR}/${job_name}_%j.out
#SBATCH --error=${LOG_DIR}/${job_name}_%j.err

set -euo pipefail

echo "Job ID: \${SLURM_JOB_ID}"
echo "Node: \$(hostname)"
echo "CPUs per task: \${SLURM_CPUS_PER_TASK:-unknown}"
echo "CPUs on node allocation: \${SLURM_CPUS_ON_NODE:-unknown}"
echo "Memory per node: \${SLURM_MEM_PER_NODE:-unknown}"
echo "Input: ${h_dir}"
echo "Output: ${output_root}"
echo "Local extraction workers: ${LOCAL_WORKERS}"
echo "Started: \$(date --iso-8601=seconds)"

"${PYTHON}" "${NODE_PIPELINE}" \
    --pipeline-script "${ORIGINAL_PIPELINE}" \
    node "${h_dir}" \
    --output-root "${output_root}" \
    --pattern "${PATTERN}" \
    --batch-size "${BATCH_SIZE}" \
    --local-workers "${LOCAL_WORKERS}" \
    --python-executable "${PYTHON}" \
    --reuse-manifest \
    --skip-completed

echo "Finished: \$(date --iso-8601=seconds)"
EOF

        chmod +x "${slurm_file}"
        echo "Created: ${slurm_file}"

        if [[ "${SUBMIT_JOBS}" -eq 1 ]]; then
            sbatch "${slurm_file}"
        fi

        job_count=$((job_count + 1))
    done
done

echo
echo "Generated ${job_count} H-folder jobs."
echo "Scripts: ${JOB_SCRIPT_DIR}"
echo "Logs: ${LOG_DIR}"
