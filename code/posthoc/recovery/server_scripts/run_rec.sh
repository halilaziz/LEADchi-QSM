#!/bin/bash
# LEAD-chi POST HOC (2026-10-07): QSMxT re-run of an injected BIDS tree with the registered PRIMARY settings
# (V-SHARP 12->1.5 mm, MEDI lambda from amendment 1, Mask 3, same config TOML). No R2* (magnitude unchanged).
# Usage: run_rec.sh <variant>   reads $LX/posthoc/recovery/bids_<variant>, writes $LX/qsm/<variant>
source /path/to/Penn_LEAD/LEADchi/scripts/pipeline/env.sh
export LX_NICE=19 LX_MAXJOBS=16 LX_THREADS=2 OMP_NUM_THREADS=2 RAYON_NUM_THREADS=2
set -euo pipefail
V=$1; BIDS=$LX/posthoc/recovery/bids_$V; OUT=$LX/qsm/$V
case "$OUT" in $LX/qsm/rec_*|$LX/qsm/e2_*) ;; *) echo "bad variant $V"; exit 2 ;; esac
[ -d "$BIDS" ] || { echo "$BIDS missing"; exit 2; }
[ -e "$OUT" ] && { echo "$OUT exists"; exit 2; }
LAMBDA=$(head -1 $LX/prereg/amendment1_lambda_qsmxt.txt)
ARGS=(--phase-offset-removal true --phase-offset-sigma 4 4 4 --unwrapping-algorithm romeo
      --no-romeo-individual --romeo-template 2 --b0-estimation weighted-avg --b0-weight-type phase-snr
      --qsm-reference none --output-space acquired --obliquity-threshold=-1
      --config $LX/prereg/qsm_primary.toml --n-procs 16 --mem-limit-gb 200 --use-custom-masks mask3
      --bf-algorithm vsharp --vsharp-max-radius 12 --vsharp-min-radius 1.5 --vsharp-threshold 0.05
      --qsm-algorithm medi --medi-lambda "$LAMBDA" --medi-merit false --medi-data-weighting 1
      --medi-percentage 0.3 --medi-max-iter 30 --medi-cg-max-iter 10 --medi-tol 0.1)
mkdir -p "$OUT"
echo "$(date -u +%FT%TZ) start $V lambda $LAMBDA" >> $LX/logs/posthoc.log
nice -n 19 "$QSMXT_BIN" run "$BIDS" "$OUT" "${ARGS[@]}" > $LX/logs/posthoc_qsmxt_$V.log 2>&1
echo "$(date -u +%FT%TZ) qsmxt $V rc=$?" >> $LX/logs/posthoc.log
