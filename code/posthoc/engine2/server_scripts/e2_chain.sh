#!/bin/bash
# LEAD-chi POST HOC (2026-10-07): engine-2 (QSM.jl) full run on all coded sessions and on the -14.24 ppb injected tree,
# then the registered plausibility + ROI scripts with LX pointed to the shadow tree (nothing written to frozen folders).
source /path/to/Penn_LEAD/LEADchi/scripts/pipeline/env.sh
source $LX/tools/julia/julia_env.sh
LXR=$LX; D=$LXR/posthoc/engine2; SH=$LXR/posthoc/shadow; mkdir -p $D/logs
st() { echo "$(date -u +%FT%TZ) $*" >> $D/e2_chain.log; }
st "e2 chain start"
$PY $D/scripts/e2_io.py prepare --tsv $D/sessions_all.tsv --work $D/work_all >> $D/e2_chain.log 2>&1
$PY $D/scripts/e2_io.py prepare --tsv $D/sessions_rec_m14.tsv --work $D/work_rec_m14 --bids $LXR/posthoc/recovery/bids_rec_m14 >> $D/e2_chain.log 2>&1
for set in all rec_m14; do
  for w in 1 2 3 4 5 6; do
    nice -n 19 env JULIA_NUM_THREADS=4 $JULIA_BIN $D/scripts/e2_qsmjl.jl $D/sessions_$set.tsv $w 6 $D/work_$set > $D/logs/${set}_w$w.log 2>&1 &
  done
  wait; st "julia $set done: $(ls $D/work_$set/*/done 2>/dev/null | wc -l) sessions; errors $(cat $D/logs/${set}_w*.log | grep -c error)"
done
export LX=$SH
$PY $D/scripts/e2_io.py finalize --tsv $D/sessions_all.tsv --work $D/work_all --prefix e2_lap_vsharp >> $D/e2_chain.log 2>&1; st "finalize all rc=$?"
$PY $D/scripts/e2_io.py finalize --tsv $D/sessions_rec_m14.tsv --work $D/work_rec_m14 --prefix e2rec_m14 >> $D/e2_chain.log 2>&1; st "finalize rec rc=$?"
for v in e2_lap_vsharp_rts e2_lap_vsharp_tv e2rec_m14_rts e2rec_m14_tv; do
  nice -n 19 $PY $LX_SCRIPTS/qsm_plausibility.py --stage indices --variant $v --jobs 16 > $D/logs/idx_$v.log 2>&1; st "idx $v rc=$?"
  nice -n 19 $PY $LX_SCRIPTS/qsm_plausibility.py --stage flags --variant $v > $D/logs/flg_$v.log 2>&1; st "flg $v rc=$?"
  nice -n 19 $PY $LX_SCRIPTS/extract_rois.py --stage roi --variant $v --jobs 16 > $D/logs/roi_$v.log 2>&1; st "roi $v rc=$?"
done
st "e2 chain done"
