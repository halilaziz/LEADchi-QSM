#!/usr/bin/env bash
# run_tests_SIMULATED.sh - end-to-end test of analyze_kit2.py on SIMULATED ratings only (never real ratings).
# Writes everything into analysis/test_run_SIMULATED/ (must not exist yet).
set -u
cd "$(dirname "$0")"
T=test_run_SIMULATED
[ -e "$T" ] && { echo "$T exists; nothing was done"; exit 1; }
mkdir -p "$T/inputs"
echo "ALL FILES IN THIS FOLDER ARE SIMULATED TEST DATA, NOT REAL RATINGS. Built $(date -u +%FT%TZ) by run_tests_SIMULATED.sh." > "$T/README_SIMULATED.txt"
python3 simulate_kit2_ratings.py --outdir $T/inputs --seed 777 --blank 1
python3 simulate_kit2_ratings.py --outdir $T/inputs --seed 778 --blank 0
python3 - <<'EOF'
import csv
K = "/path/to/Penn_LEAD/LEADchi/human_qc_kit/"
for key, idc, out in ((K + "_key_DO_NOT_SEND/kit_key.csv", "kit_id", "SIMULATED_kit1_identity_AI.csv"),
                      (K + "kit2_2026-10-28/_key_DO_NOT_SEND/kit2_key.csv", "kit2_id", "SIMULATED_kit2_identity_AI.csv")):
    with open("test_run_SIMULATED/inputs/" + out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["kit_id", "grade", "failure_codes", "placement", "note", "timestamp", "rater"])
        for r in csv.DictReader(open(key)):
            w.writerow([r[idc], r["ai_grade"], "", r["ai_placement"], "SIMULATED identity", "", "SIM"])
EOF
I=$T/inputs
{
echo "=== SIMULATED tests of analyze_kit2.py ($(date -u +%FT%TZ)); script sha256 $(shasum -a 256 analyze_kit2.py | cut -c1-64)"
echo "--- 1. simulated ratings, one kit-2 item blank (seed 777): runs; complete_human_qc must be false"
python3 analyze_kit2.py --simulated --kit1-ratings $I/SIMULATED_kit1_ratings_seed777.csv --kit2-ratings $I/SIMULATED_kit2_ratings_seed777_blank1.csv --out $T/run_seed777_blank1
echo "--- 2. complete simulated ratings (seed 778)"
python3 analyze_kit2.py --simulated --kit1-ratings $I/SIMULATED_kit1_ratings_seed778.csv --kit2-ratings $I/SIMULATED_kit2_ratings_seed778_blank0.csv --out $T/run_seed778_complete
echo "--- 3. identity (human = AI everywhere): K1_H_replace and H_all180 must equal AI_frozen (n 80, beta 4.96280154409223); retest 1.0"
python3 analyze_kit2.py --simulated --kit1-ratings $I/SIMULATED_kit1_identity_AI.csv --kit2-ratings $I/SIMULATED_kit2_identity_AI.csv --out $T/run_identity
echo "--- 4. cross-check: kit-1 analyze_human_qc.py on the seed-777 kit-1 file; its H_replace must equal K1_H_replace of test 1"
python3 -B ../../analysis/analyze_human_qc.py --simulated --ratings $I/SIMULATED_kit1_ratings_seed777.csv --out $T/crosscheck_kit1_analyzer_seed777 --nextrap 0 --nboot 200 > $T/crosscheck_kit1_analyzer_seed777_stdout.txt 2>&1
python3 - <<'EOF'
import json, pandas as pd
T = "test_run_SIMULATED/"
a = json.load(open(T + "crosscheck_kit1_analyzer_seed777/SIMULATED_m1_refit.json"))["scenarios"]["H_replace"]
b = json.load(open(T + "run_seed777_blank1/SIMULATED_kit2_m1_refit.json"))["scenarios"]["K1_H_replace"]
for k in ("n", "beta_PS_ppb", "se_hc3", "p_two"):
    print(k, a[k], b[k], "EQUAL" if a[k] == b[k] else "DIFFERENT")
print("first_usable_keys identical:", a["first_usable_keys"] == b["first_usable_keys"])
s1 = pd.read_csv(T + "crosscheck_kit1_analyzer_seed777/SIMULATED_scenario_usable.csv")
s2 = pd.read_csv(T + "run_seed777_blank1/SIMULATED_kit2_scenario_usable.csv")
m = s1.merge(s2, on="sess_key")
print("usable/first_usable identical over", len(m), "rows:",
      bool((m["usable__H_replace"] == m["usable__K1_H_replace"]).all() and
           (m["first_usable__H_replace"] == m["first_usable__K1_H_replace"]).all()))
EOF
echo "--- 5. kit-1 CSV given as the kit-2 file: must stop and write nothing"
python3 analyze_kit2.py --simulated --kit1-ratings $I/SIMULATED_kit1_ratings_seed778.csv --kit2-ratings $I/SIMULATED_kit1_ratings_seed778.csv --out $T/run_wrongfile
echo "exit $?; output folder created: $(test -d $T/run_wrongfile && echo yes || echo no)"
echo "--- 6. --simulated with the default (real) kit-1 CSV: must stop before parsing any rating file"
python3 analyze_kit2.py --simulated --kit2-ratings $I/SIMULATED_kit2_ratings_seed778_blank0.csv --out $T/run_realk1
echo "exit $?; output folder created: $(test -d $T/run_realk1 && echo yes || echo no)"
echo "--- 7. output folder that already holds results: must stop"
python3 analyze_kit2.py --simulated --kit1-ratings $I/SIMULATED_kit1_ratings_seed778.csv --kit2-ratings $I/SIMULATED_kit2_ratings_seed778_blank0.csv --out $T/run_seed778_complete
echo "exit $?"
} > $T/SIMULATED_test_log.txt 2>&1
cat $T/SIMULATED_test_log.txt
