#!/usr/bin/env bash
# Rebuild every demonstration output of the LEAD-chi analysis code from scratch, then write SHA256SUMS.
# Synthetic outcomes only. About 20-30 min on the Mac (4 R processes in parallel), plus the guard tests.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p synthetic demo_output sim
Rscript make_synthetic.R --out synthetic
for sc in null delta-0.6 null_dominant; do
  Rscript run_all.R --input "synthetic/SYNTHETIC_roi_table_${sc}.csv" --gates "synthetic/SYNTHETIC_gates_${sc}.json" \
    --multiverse "synthetic/SYNTHETIC_multiverse_${sc}.csv" --out "demo_output/${sc}" --label "${sc}" > "demo_output/${sc}.log" 2>&1 &
done
Rscript sim_typeI_power.R --out sim > sim/sim.log 2>&1 &
wait
for sc in null delta-0.6 null_dominant; do
  python3 verify_statsmodels.py --input "synthetic/SYNTHETIC_roi_table_${sc}.csv" --gates "synthetic/SYNTHETIC_gates_${sc}.json" \
    --results "demo_output/${sc}/results.json" --out "demo_output/${sc}/verify_statsmodels.json" > "demo_output/${sc}_verify.log" 2>&1
done
Rscript tests/test_guards.R > demo_output/test_guards.log 2>&1
find . -type f ! -name SHA256SUMS ! -name '.DS_Store' ! -path './__pycache__/*' | LC_ALL=C sort | xargs shasum -a 256 > SHA256SUMS
echo "done; see SHA256SUMS"
