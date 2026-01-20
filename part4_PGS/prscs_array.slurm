#!/bin/bash

#SBATCH --job-name=PRScs_time
#SBATCH --partition=CPU
#SBATCH --array=1-22%5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

CHR=${SLURM_ARRAY_TASK_ID}      # 1-22
TIMEPOINT="time"                # e.g. 7years, 8years

PRScs_py="/path/to/PRScs/PRScs.py"
REF_DIR="/path/to/ldblk_ukbb_eur"

# Genotype
BIM_PREFIX="/path/to/PRScs_input/PRScs_input_chr${CHR}"

# GWAS summary stats for PRS-CS
SST_FILE="/path/to/sumstats/PRSCS_sumstats_${TIMEPOINT}.txt"

N_GWAS=samplesize

phi_list=(1e-2)

# Fixed seed
SEED=seed_ID${CHR}

# Output root directory
OUT_ROOT="/path/to/PRS_result/PRScs/${TIMEPOINT}"

for PHI in "${phi_list[@]}"; do

  PHI_LABEL=$(echo "${PHI}" | sed 's/-/_/g')

  OUT_DIR="${OUT_ROOT}/phi${PHI_LABEL}"
  mkdir -p "${OUT_DIR}"

  echo "[$(date +'%F %T')] Running PRS-CS: chr=${CHR}, phi=${PHI}, seed=${SEED}"

  python "${PRScs_py}" \
    --ref_dir    "${REF_DIR}" \
    --bim_prefix "${BIM_PREFIX}" \
    --sst_file   "${SST_FILE}" \
    --n_gwas     "${N_GWAS}" \
    --chrom      "${CHR}" \
    --out_dir    "${OUT_DIR}" \
    --phi        "${PHI}" \
    --seed       "${SEED}"

  echo "[$(date +'%F %T')] Finished: chr=${CHR}, phi=${PHI}"
done
