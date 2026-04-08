#!/usr/bin/env bash
set -euo pipefail

# Batch driver for the AR(2) case study.
# Runs 50 replicates for each of 4 scenarios, with K in {5,10,20}.
#
# Directory layout:
# case_study_pilot/
#   iid/rep_001/
#   case1/rep_001/
#   case2/rep_001/
#   case3/rep_001/
#   ...

N_REPS=50
T_VAL=100000
P_VAL=10
K_LIST="5,10,20"

SIM_SCRIPT="simulate_data.r"
FIT_SCRIPT="fit_posteriors.r"
WASS_SCRIPT="wasserstein_average.r"
STAN_FILE="model_ar2_errors.stan"

OUT_ROOT="case_study"

# MCMC settings: tune these as needed
CHAINS=4
PARALLEL_CHAINS=4
ITER_WARMUP=500
ITER_SAMPLING=500
ADAPT_DELTA=0.9
MAX_TREEDEPTH=12

mkdir -p "${OUT_ROOT}"

run_one_rep() {
  local scenario="$1"
  local phi1="$2"
  local phi2="$3"
  local rep="$4"

  local rep_id
  rep_id=$(printf "rep_%03d" "${rep}")

  local out_dir="${OUT_ROOT}/${scenario}/${rep_id}"
  local data_file="${out_dir}/sim_data.rds"
  local seed=$((20270000 + 1000 * rep))

  mkdir -p "${out_dir}"

  echo "=================================================="
  echo "Scenario: ${scenario} | Replicate: ${rep_id}"
  echo "phi1=${phi1}, phi2=${phi2}, seed=${seed}"
  echo "Output dir: ${out_dir}"
  echo "=================================================="

  Rscript "${SIM_SCRIPT}" \
    --T "${T_VAL}" \
    --p "${P_VAL}" \
    --phi1 "${phi1}" \
    --phi2 "${phi2}" \
    --seed "${seed}" \
    --output "${data_file}"

  Rscript "${FIT_SCRIPT}" \
    --data "${data_file}" \
    --stan-file "${STAN_FILE}" \
    --output-dir "${out_dir}" \
    --K "${K_LIST}" \
    --chains "${CHAINS}" \
    --parallel-chains "${PARALLEL_CHAINS}" \
    --iter-warmup "${ITER_WARMUP}" \
    --iter-sampling "${ITER_SAMPLING}" \
    --seed "${seed}" \
    --adapt-delta "${ADAPT_DELTA}" \
    --max-treedepth "${MAX_TREEDEPTH}"

  for K in 5 10 20; do
    Rscript "${WASS_SCRIPT}" \
      --input "${out_dir}/subset_posteriors_K${K}.rds" \
      --output "${out_dir}/wasserstein_beta_K${K}.rds"
  done
}

# Scenarios from the paper excerpt
for rep in $(seq 14 "${N_REPS}"); do
  run_one_rep "iid"   "0.0" "0.0" "${rep}"
  run_one_rep "case1" "0.3" "0.1" "${rep}"
  run_one_rep "case2" "0.8" "0.0" "${rep}"
  run_one_rep "case3" "0.4" "0.4" "${rep}"
done

echo
echo "Batch run complete."
echo "Outputs written under: ${OUT_ROOT}"