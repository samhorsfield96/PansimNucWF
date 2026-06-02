#!/bin/bash
# =============================================================================
# Slurm array job – PansimNucWF
#
# Submit (adjust --array upper bound to N_SAMPLES-1):
#   sbatch --array=0-99 run_lhs_sweep.slurm
#
# Or set the array range here and submit without the flag:
#   sbatch run_lhs_sweep.slurm
# =============================================================================

#SBATCH --job-name=pansimnucWF_lhs
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
##SBATCH --array=0-99%5              # <-- Set to 0-(N_SAMPLES-1)
#SBATCH --array=0              # Testing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        # Should match threads= in your .conf
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=normal.168h   # Uncomment and set if required
##SBATCH --mail-user=samuel.horsfield@unine.ch       # Uncomment and set if required
##SBATCH --mail-type=END,FAIL

# =============================================================================
# User configuration – edit these paths before submitting
# =============================================================================

# snakemake directory
SNAKEMAKE_DIR="/home/sam/Software/PansimNucWF"

# Base output directory for PansimNuc results (--outdir-base from lhs script)
# Each run writes to $OUTDIR_BASE/run_<id>/
OUTDIR_BASE="/scratch/sam_simulations/lhs_parameter_sweep/exon_doubexp_selection_one_split_gen_500"

# =============================================================================

set -euo pipefail

RUN_ID=${SLURM_ARRAY_TASK_ID}
RUN_CONFIG="${SNAKEMAKE_DIR}/config/config.yaml"
RUN_INPUT_DIR="${OUTDIR_BASE}/run_${RUN_ID}"
RUN_OUTDIR="${RUN_INPUT_DIR}/PansimNucWF"

# Sanity checks
if [[ ! -f "$RUN_CONFIG" ]]; then
    echo "ERROR: Config file not found: $RUN_CONFIG" >&2
    exit 1
fi

mkdir -p "${RUN_OUTDIR}"

echo "========================================================"
echo "Job ID   : ${SLURM_JOB_ID:-local}_${RUN_ID}"
echo "Run ID   : ${RUN_ID}"
echo "Config   : ${RUN_CONFIG}"
echo "Output   : ${RUN_OUTDIR}"
echo "Started  : $(date -Is)"
echo "========================================================"

cd ${SNAKEMAKE_DIR}

pixi run snakemake -s ${SNAKEMAKE_DIR}/Snakefile --directory ${RUN_OUTDIR} --configfile ${RUN_CONFIG} --cores $SLURM_CPUS_PER_TASK --use-conda --conda-prefix ${SNAKEMAKE_DIR}/.snakemake/conda --config reference=${RUN_INPUT_DIR}/root.fasta input_dir=${RUN_INPUT_DIR} output_dir=${RUN_OUTDIR} simulated="True" --unlock
pixi run snakemake -s ${SNAKEMAKE_DIR}/Snakefile --directory ${RUN_OUTDIR} --configfile ${RUN_CONFIG} --cores $SLURM_CPUS_PER_TASK --use-conda --conda-prefix ${SNAKEMAKE_DIR}/.snakemake/conda --config reference=${RUN_INPUT_DIR}/root.fasta input_dir=${RUN_INPUT_DIR} output_dir=${RUN_OUTDIR} simulated="True" > "${RUN_OUTDIR}/pansimnucWF.log" 2>&1

EXIT_CODE=$?

echo "Finished : $(date -Is)"
echo "Exit code: ${EXIT_CODE}"

exit ${EXIT_CODE}
