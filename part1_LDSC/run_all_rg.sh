#!/bin/bash

# Path to LDSC script and input data
LDSC_PATH="path/to/ldsc.py"
INPUT_DIR="path/to/input_data/"
REF_LD_DIR="path/to/reference_LD/"

# Output directory for results
OUT_DIR="./rg_results"

mkdir -p "$OUT_DIR"

# ========== input files ==========
FILES=(
  "childhood_bmi_nature_metab_2022_birth.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_6weeks.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_3months.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_6months.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_8months.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_1year.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_1.5years.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_2years.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_3years.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_5years.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_7years.sumstats.sumstats.gz"
  "childhood_bmi_nature_metab_2022_8years.sumstats.sumstats.gz"
  "adult_Locke_UKB_BMI.sumstats.sumstats.gz"
)

# Loop through the files to run LDSC
for (( i=0; i<${#FILES[@]}; i++ )); do
  for (( j=i+1; j<${#FILES[@]}; j++ )); do
    
    FILE1=${FILES[$i]}
    FILE2=${FILES[$j]}

    # Extract names without extensions
    NAME1=$(basename "$FILE1" .sumstats.sumstats.gz | sed 's/.sumstats$//')
    NAME2=$(basename "$FILE2" .sumstats.sumstats.gz | sed 's/.sumstats$//')

    OUT_FILE="${OUT_DIR}/rg_${NAME1}_vs_${NAME2}"

    echo "Analyzing: $NAME1 vs $NAME2"
    
    # Run LDSC command
    python "$LDSC_PATH" \
      --rg "${INPUT_DIR}${FILE1},${INPUT_DIR}${FILE2}" \
      --ref-ld-chr "$REF_LD_DIR" \
      --w-ld-chr "$REF_LD_DIR" \
      --out "$OUT_FILE"
    
    echo "Completed: $OUT_FILE.log"
    echo "--------------------------------------"
  done
done

echo "Done"
