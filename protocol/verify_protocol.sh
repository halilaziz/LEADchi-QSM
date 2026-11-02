#!/usr/bin/env bash
# Verify the internally time-stamped LEAD-chi protocol files and the frozen analysis code.
# Usage (from the repository root or from protocol/):  bash protocol/verify_protocol.sh
set -u
cd "$(dirname "$0")"
fail=0
echo "== 1. Original per-document freeze records (*.sha256; the 'frozen <UTC>' line is a time note, not a checksum)"
for f in LEADchi_amendment_*.sha256 LEADchi_prechi_freeze.sha256; do
  shasum -a 256 -c "$f" 2>/dev/null || fail=1
done
echo "== 2. Protocol v1.0 FINAL against the original freeze manifest (MANIFEST.sha256)"
grep ' \./LEADchi_prereg_v1.0_FINAL.md$' MANIFEST.sha256 | shasum -a 256 -c - || fail=1
echo "== 3. Protocol FINAL and amendment 0 against the pre-chi freeze (LEADchi_prechi_freeze.md)"
grep -E '^[0-9a-f]{64}  (LEADchi_prereg_v1.0_FINAL.md|LEADchi_amendment_0_route.md)$' LEADchi_prechi_freeze.md | shasum -a 256 -c - || fail=1
echo "== 4. Repository checksum list for every file in protocol/ (SHA256SUMS, made 2026-11-02 from the byte-identical copies)"
shasum -a 256 -c SHA256SUMS || fail=1
echo "== 5. Frozen registered analysis code (code/analysis) against MANIFEST.sha256"
# secondary.R and verify_statsmodels.py were patched after the freeze (deviation C4); their registered
# versions are kept in code/analysis/deviations/. R/synthetic.R was sanitized (one path) and is skipped.
grep ' \./code/analysis/' MANIFEST.sha256 \
  | grep -v -E 'demo_output|/SHA256SUMS$|/R/synthetic.R$' \
  | sed -e 's# \./code/analysis/R/secondary.R$# ./deviations/secondary_v1.0_registered.R#' \
        -e 's# \./code/analysis/verify_statsmodels.py$# ./deviations/verify_statsmodels_v1.0_registered.py#' \
        -e 's# \./code/analysis/# ./#' \
  | (cd ../code/analysis && shasum -a 256 -c -) || fail=1
echo "== 6. Data freeze manifest against the unlock record"
want=$(grep -o '"freeze": "[0-9a-f]*"' ../results/confirmatory/UNLOCK_RECORD.json | grep -o '[0-9a-f]\{64\}')
have=$(shasum -a 256 freeze_manifest.json | cut -c1-64)
if [ "$want" = "$have" ]; then echo "freeze_manifest.json: OK"; else echo "freeze_manifest.json: FAILED"; fail=1; fi
if [ $fail -eq 0 ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED"; fi
exit $fail
