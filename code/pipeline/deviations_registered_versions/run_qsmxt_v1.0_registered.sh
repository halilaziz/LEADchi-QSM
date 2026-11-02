#!/usr/bin/env bash
# run_qsmxt.sh - LEAD-chi QSM, R2* and T2* with the pinned QSMxT binary (task P5b, 2026-09-25).
# Written and hashed before the OSF registration; NOT executed before the approvals. Runs on server only.
#
# Registered pipeline (registration Annex B.2; SINGLE PASS on Mask 3, no --two-pass): phase-offset removal
# (MCPC-3D-S type, sigma 4 4 4 voxels) -> ROMEO in template mode (template echo 2, 1-indexed) -> weighted echo
# combination (weighted-avg, phase-snr weights w = |S|*TE) -> V-SHARP 12 -> 1.5 mm, TSVD 0.05 -> MEDI (lambda
# from the per-engine L-curve, amendment 1; SMV off via the TOML; MERIT off; data weighting 1; edge fraction 0.3;
# max-iter 30, CG max-iter 10, tol 0.1; no CSF term) -> no internal reference (two-step whole-brain referencing on
# the support S is done in qsm_plausibility.py / extract_rois.py). Oblique resampling disabled
# (--obliquity-threshold -1; chi is never resampled). Split-echo {2,4}: QSM.rs converts the ppm field to radians
# with the FIRST echo time of the input set (pipeline/inversion.rs, v0.38.0), so TE1 = 15 ms doubles the radian
# scale; lambda for split24 is therefore 2 x lambda (quadratic data term, L1 regulariser); split13 (TE1 = 7.5 ms)
# uses lambda unchanged.
# B0: QSMxT `run` has no B0 flag; it reads MagneticFieldStrength from each coded sidecar, which blind_setup.py
# wrote as round(f0/42.576, 4) T (E3 finding A1). Phase is already float32 radians (E3 finding A3).
# R2*/T2* (ARLO, --do-r2starmap/--do-t2starmap) are produced by the primary and pilot runs only.
#
# Usage: run_qsmxt.sh --stage STAGE [--variant NAME] [--dry-run]
#   install      A1: download the pinned tarball, check its SHA-256, unpack into $LX/tools/qsmxt, capture
#                `qsmxt --version`, `qsmxt run --help`, `qsmxt init`, `qsmxt quality-map --help`
#   config       write prereg/qsm_primary.toml (defaults with MEDI SMV off) and prereg/pipeline_choices.json
#   subsets      symlink BIDS trees inside blind/: bids_pilot, bids_lcurve, bids_gate20, bids_split13, bids_split24
#   negate       sign-check branch only (pilot_checks.py sign = "flip"): write blind/bids_all_neg with every phase
#                echo multiplied by -1 (headers unchanged) and record prereg/phase_sign.txt; then rebuild the subsets
#                and recompute every run (the pilot chi is discarded)
#   mask3        build the custom Mask 3 per coded session (MASK_MODE custom3: SynthStrip Mask 1, QSMxT phase-quality
#                Mask 2 = Otsu threshold within Mask 1, Mask 3 = fill-holes(Mask 1 AND Mask 2) minus the first and
#                last 2 slab partitions)
#   plan         `qsmxt validate` + `qsmxt run --dry` for every variant (creates nothing in the output dirs)
#   run          run --variant pilot | lcurve | g2 | g2alt | primary | split | multiverse | fallback | <single variant>
# Every real run requires prereg/qsmxt_flags_verified.txt (user sign-off after reading the captured help, before
# posting) and aborts if any flag used here is absent from prereg/qsmxt_run_help.txt.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=env.sh
source "$HERE/env.sh"

STAGE="" VARIANT="" DRY=0
usage() { sed -n '2,27p' "$0"; }
while [ $# -gt 0 ]; do
  case "$1" in
    --stage) STAGE=$2; shift 2 ;;
    --variant) VARIANT=$2; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[ "$DRY" = 1 ] && export LX_DRY=1
[ -n "$STAGE" ] || { usage; exit 2; }

# ---- pins (E3 section B; protocol sec. 12.1) ---------------------------------------------------------------
QSMXT_TARBALL_NAME=qsmxt-$QSMXT_TAG-x86_64-unknown-linux-musl.tar.gz
QSMXT_URL=https://github.com/QSMxT/QSMxT/releases/download/$QSMXT_TAG/$QSMXT_TARBALL_NAME
QSMXT_TARBALL_SHA256=6593438c029749b5e3fb5c78ce485f6f5e282545b1892832cf30539b0326ca49
HELP_TXT=$LX/prereg/qsmxt_run_help.txt
QMAP_HELP_TXT=$LX/prereg/qsmxt_quality_map_help.txt
DEFAULTS_TOML=$LX/prereg/qsmxt_${QSMXT_TAG}_defaults.toml
CONFIG_TOML=$LX/prereg/qsm_primary.toml
CHOICES=$LX/prereg/pipeline_choices.json
VERIFIED=$LX/prereg/qsmxt_flags_verified.txt
QMAP_CMD_FILE=$LX/prereg/qsmxt_quality_map_cmd.txt   # one line with {MAG1..4} {PHASE1..4} {OUT} placeholders
LAMBDA_FILE=$LX/prereg/amendment1_lambda_qsmxt.txt   # written after the L-curve (amendment 1)
LAMBDA_DEFAULT=0.001                                 # Cornell lambda 1000 on the QSM.rs scale (E3 A7, inferred)
LAMBDA_GRID=(0.004 0.002 0.001 0.0005)               # Cornell {250, 500, 1000, 2000} (E3 A7); names match pilot_checks.py
MASK_MODE=${MASK_MODE:-custom3}                      # custom3 (registered) | preset (registered technical fallback)
N_PARTITIONS=96                                      # slab partitions (matrix 120 x 160 x 96); edge trim 2 each side
EDGE_TRIM=2
MEM_GB=${LX_MEMGB:-200}
BLIND=$LX/blind
SIGN_FILE=$LX/prereg/phase_sign.txt
BIDS_ALL=$BLIND/bids_all
[ -s "$SIGN_FILE" ] && [ "$(head -1 "$SIGN_FILE")" = negated ] && BIDS_ALL=$BLIND/bids_all_neg
BF_LEVELS=(vsharp pdf lbv)
INV_LEVELS=(medi ilsqr rts)

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

current_lambda() {
  if [ -s "$LAMBDA_FILE" ]; then head -1 "$LAMBDA_FILE" | tr -d '[:space:]'
  elif [ "$DRY" = 1 ]; then echo "<lambda-from-amendment-1>"
  elif [ "$STAGE" = plan ]; then echo "$LAMBDA_DEFAULT"          # plan (--dry of QSMxT) only; never a real run
  else die "$LAMBDA_FILE missing: the L-curve rule (amendment 1) must fix lambda before this run"; fi
}

# ---- variant table -----------------------------------------------------------------------------------------
# sets V_BIDS V_OUT V_BF V_INV V_LAMBDA V_R2S
variant_spec() {
  local v=$1
  V_R2S=0; V_BF=vsharp; V_INV=medi; V_BIDS=$BIDS_ALL
  case "$v" in
    pilot)     V_BIDS=$BLIND/bids_pilot;  V_OUT=$LX/qsm/pilot;   V_LAMBDA=$LAMBDA_DEFAULT; V_R2S=1 ;;
    lcurve_*)  V_BIDS=$BLIND/bids_lcurve; V_OUT=$LX/qsm/$v; V_LAMBDA=${v#lcurve_} ;;
    g2)        V_BIDS=$BLIND/bids_gate20; V_OUT=$LX/qsm/g2; V_LAMBDA=$(current_lambda) ;;
    g2_mv_*_*) V_BIDS=$BLIND/bids_gate20; V_BF=$(echo "$v" | cut -d_ -f3); V_INV=$(echo "$v" | cut -d_ -f4)
               V_OUT=$LX/qsm/$v; V_LAMBDA=$(current_lambda) ;;
    primary)   V_OUT=$LX/qsm/primary; V_LAMBDA=$(current_lambda); V_R2S=1 ;;
    split13)   V_BIDS=$BLIND/bids_split13; V_OUT=$LX/qsm/split13; V_LAMBDA=$(current_lambda) ;;
    split24)   V_BIDS=$BLIND/bids_split24; V_OUT=$LX/qsm/split24
               V_LAMBDA=$(current_lambda); case "$V_LAMBDA" in \<*) ;; *) V_LAMBDA=$(awk -v l="$V_LAMBDA" 'BEGIN{printf "%.6g", 2*l}') ;; esac ;;
    mv_*_*)    V_BF=$(echo "$v" | cut -d_ -f2); V_INV=$(echo "$v" | cut -d_ -f3)
               V_OUT=$LX/qsm/$v; V_LAMBDA=$(current_lambda)
               [ "$v" = mv_vsharp_medi ] && die "mv_vsharp_medi is the primary run" ;;
    *) die "unknown variant $v" ;;
  esac
  return 0
}

choice_get() {   # choice_get KEY : one value from prereg/pipeline_choices.json (written by --stage config)
  if [ ! -s "$CHOICES" ]; then [ "$DRY" = 1 ] && { echo "<$1-from-config>"; return 0; }; die "$CHOICES missing: run --stage config"; fi
  "$SYSPY" -c 'import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])' "$CHOICES" "$1"
}

build_args() {   # fills the global array ARGS for the current V_* values
  ARGS=(--phase-offset-removal true --phase-offset-sigma 4 4 4 --unwrapping-algorithm romeo
        --no-romeo-individual --romeo-template 2
        --b0-estimation weighted-avg --b0-weight-type phase-snr
        --qsm-reference none --output-space acquired --obliquity-threshold -1
        --config "$CONFIG_TOML" --n-procs "$LX_MAXJOBS" --mem-limit-gb "$MEM_GB")
  case "$MASK_MODE" in
    custom3) ARGS+=(--use-custom-masks mask3) ;;
    preset) ARGS+=(--mask-preset bet-and-phase) ;;
    *) die "MASK_MODE must be custom3 or preset" ;;
  esac
  case "$V_BF" in
    vsharp) ARGS+=(--bf-algorithm vsharp --vsharp-max-radius 12 --vsharp-min-radius 1.5 --vsharp-threshold 0.05) ;;
    pdf)    ARGS+=(--bf-algorithm pdf --pdf-tol "$(choice_get pdf_tol)") ;;
    lbv)    ARGS+=(--bf-algorithm lbv --lbv-tol "$(choice_get lbv_tol)") ;;
    *) die "background-removal level $V_BF not registered" ;;
  esac
  case "$V_INV" in
    medi)  ARGS+=(--qsm-algorithm medi --medi-lambda "$V_LAMBDA" --medi-merit false --medi-data-weighting 1
                  --medi-percentage 0.3 --medi-max-iter 30 --medi-cg-max-iter 10 --medi-tol 0.1) ;;
    ilsqr) ARGS+=(--qsm-algorithm ilsqr) ;;
    rts)   ARGS+=(--qsm-algorithm rts) ;;
    *) die "inversion level $V_INV not registered" ;;
  esac
  [ "$V_R2S" = 1 ] && ARGS+=(--do-r2starmap --do-t2starmap)
  return 0
}

check_flags() {   # every --flag in ARGS must appear in the captured help of the pinned binary
  local a missing=0
  if [ ! -s "$HELP_TXT" ]; then
    [ "$DRY" = 1 ] && { log "flag check skipped (no $HELP_TXT yet; dry-run)"; return 0; }
    die "$HELP_TXT missing: run --stage install first"
  fi
  for a in "${ARGS[@]}"; do
    case "$a" in --*) grep -qE -- "(^|[[:space:],])${a}([[:space:]=,<]|\$)" "$HELP_TXT" \
                        || { echo "flag not in qsmxt run --help: $a" >&2; missing=1; } ;; esac
  done
  [ "$missing" = 0 ] || die "flag check failed; resolve against the help text and log it in amendment 1"
}

require_verified() {
  [ "$DRY" = 1 ] && return 0
  [ -s "$VERIFIED" ] || die "$VERIFIED missing: the user signs off the captured help/init output before posting"
  [ -s "$CONFIG_TOML" ] || die "$CONFIG_TOML missing: run --stage config"
  "$QSMXT_BIN" --version 2>&1 | grep -q "${QSMXT_TAG#v}" || die "installed qsmxt is not $QSMXT_TAG"
}

masks_present() {   # custom masks must exist for every session of the BIDS tree (QSMxT would fall back silently)
  local bids=$1 n=0 m=0 d c k
  [ "$MASK_MODE" = preset ] && return 0
  for d in "$bids"/sub-*/ses-*; do
    [ -d "$d" ] || continue
    c=$(basename "$(dirname "$d")"); k=$(basename "$d"); n=$((n + 1))
    [ -s "$BLIND/bids_all/derivatives/mask3/$c/$k/anat/${c}_${k}_mask.nii.gz" ] && m=$((m + 1))
  done
  log "custom Mask 3 present for $m/$n sessions in $(basename "$bids")"
  [ "$DRY" = 1 ] || [ "$m" = "$n" ] || die "custom masks missing; QSMxT would silently fall back (E3 A6)"
}

run_variant() {
  local v=$1
  variant_spec "$v"; build_args; check_flags; require_verified; masks_present "$V_BIDS"
  lx_assert_inside "$V_OUT"
  if [ -e "$V_OUT/derivatives/qsmxt" ] && [ "$DRY" != 1 ]; then die "$V_OUT exists; outputs are never overwritten"; fi
  lx_stamp "qsmxt_$v"
  log "variant $v: $V_BF + $V_INV (lambda $V_LAMBDA), R2* $V_R2S, input $(basename "$V_BIDS"), output $V_OUT"
  if [ "$DRY" = 1 ]; then
    lx_run "$QSMXT_BIN" run "$V_BIDS" "$V_OUT" "${ARGS[@]}"
  else
    mkdir -p "$V_OUT"
    lx_run "$QSMXT_BIN" run "$V_BIDS" "$V_OUT" "${ARGS[@]}" 2>&1 | tee "$LX/logs/qsmxt_$v.log"
    grep -i "b0.dir\|b0 direction\|B0_dir" "$LX/logs/qsmxt_$v.log" > "$LX/prereg/b0dir_$v.log" || true   # archived (Annex B.2)
    local cfg="$V_OUT/derivatives/qsmxt/pipeline_config.toml"
    [ -s "$cfg" ] && echo "$(sha256sum "$cfg" | cut -d' ' -f1)  qsm/$v/derivatives/qsmxt/pipeline_config.toml" \
      >> "$LX/prereg/frozen.sha256"
  fi
  lx_home_check "qsmxt_$v"
}

# ---- stages ------------------------------------------------------------------------------------------------
stage_install() {    # approval A1 (download of the pinned release binary; ~19 MB)
  local dst=$LX/tools/qsmxt
  lx_assert_inside "$dst"
  lx_stamp install
  lx_run curl -fL -o "$dst/$QSMXT_TARBALL_NAME" "$QSMXT_URL"
  if [ "$DRY" != 1 ]; then
    local got; got=$(sha256sum "$dst/$QSMXT_TARBALL_NAME" | cut -d' ' -f1)
    [ "$got" = "$QSMXT_TARBALL_SHA256" ] || die "tarball SHA-256 $got differs from the pinned value"
  fi
  lx_run tar -xzf "$dst/$QSMXT_TARBALL_NAME" -C "$dst"
  if [ "$DRY" = 1 ]; then
    echo "[DRY-RUN] capture qsmxt --version, run --help, init and quality-map --help into $LX/prereg/"
  else
    [ -x "$QSMXT_BIN" ] || die "qsmxt binary not found at $QSMXT_BIN after unpacking (adjust for a subfolder)"
    { echo "qsmxt $("$QSMXT_BIN" --version 2>&1) tarball_sha256 $QSMXT_TARBALL_SHA256"; } >> "$LX/prereg/software_versions.txt"
    "$QSMXT_BIN" run --help > "$HELP_TXT" 2>&1
    "$QSMXT_BIN" init > "$DEFAULTS_TOML" 2>&1
    "$QSMXT_BIN" quality-map --help > "$QMAP_HELP_TXT" 2>&1 || true
    sha256sum "$HELP_TXT" "$DEFAULTS_TOML" "$QMAP_HELP_TXT" >> "$LX/prereg/frozen.sha256"
  fi
  lx_home_check install
}

stage_config() {
  lx_require "$DEFAULTS_TOML"
  if [ "$DRY" = 1 ]; then
    echo "[DRY-RUN] write $CONFIG_TOML = $DEFAULTS_TOML with the MEDI SMV key set to false; write $CHOICES"
    return 0
  fi
  [ -e "$CONFIG_TOML" ] && die "$CONFIG_TOML exists (frozen); edit only through an amendment"
  "$SYSPY" - "$DEFAULTS_TOML" "$CONFIG_TOML" <<'PYEOF'
import re, sys
src, dst = sys.argv[1], sys.argv[2]
lines = open(src, encoding="utf-8").read().splitlines()
sec, hits = "", []
for i, l in enumerate(lines):
    m = re.match(r"\s*\[+([^\]]+)\]+\s*$", l)
    if m:
        sec = m.group(1).lower()
        continue
    if re.match(r"\s*(medi[_-]smv|use[_-]smv)\s*=\s*(true|false)\s*(#.*)?$", l, re.I) or \
       ("medi" in sec and re.match(r"\s*smv\s*=\s*(true|false)\s*(#.*)?$", l, re.I)):
        hits.append(i)
if len(hits) != 1:
    sys.exit(f"MEDI SMV key found {len(hits)} times in {src}; set it by hand, hash the file and log it (E3 F2)")
lines[hits[0]] = re.sub(r"=\s*(true|false)", "= false", lines[hits[0]], flags=re.I)
open(dst, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print(f"MEDI SMV set to false at line {hits[0] + 1}")
PYEOF
  "$SYSPY" - "$CHOICES" "$MASK_MODE" "$QSMXT_TAG" "$(sha256sum "$CONFIG_TOML" | cut -d' ' -f1)" "$DEFAULTS_TOML" <<'PYEOF'
import json, re, sys
p, mask, tag, cfg, dflt = sys.argv[1:6]
txt = open(dflt, encoding="utf-8").read()
def one(pat, what):
    hits = re.findall(pat, txt, flags=re.I | re.M)
    if len(hits) != 1:
        sys.exit(f"{what}: found {len(hits)} candidate keys in {dflt}; set it by hand, hash the file and log it in amendment 1")
    return float(hits[0])
pdf_tol = one(r"^\s*pdf[_-]?tol(?:erance)?\s*=\s*([0-9.eE+-]+)", "PDF tolerance")
lbv_tol = one(r"^\s*lbv[_-]?tol(?:erance)?\s*=\s*([0-9.eE+-]+)", "LBV tolerance")
json.dump({"mask_mode": mask, "mask_source": "qsmxt_mask" if mask == "preset" else "custom_mask3",
           "pdf_tol": pdf_tol, "lbv_tol": lbv_tol,
           "qsmxt_tag": tag, "qsm_primary_toml_sha256": cfg,
           "reference": "two-step whole-brain reference on the support S (qsm_plausibility.py, extract_rois.py)",
           "variants": ["primary", "split13", "split24"] + [f"mv_{b}_{i}" for b in ("vsharp", "pdf", "lbv")
                        for i in ("medi", "ilsqr", "rts") if (b, i) != ("vsharp", "medi")],
           "fallback_variant": "mv_pdf_ilsqr"}, open(p, "w"), indent=2)
PYEOF
  sha256sum "$CONFIG_TOML" "$CHOICES" >> "$LX/prereg/frozen.sha256"
}

mk_subset() {   # mk_subset NAME LISTFILE "new:orig new:orig ..."
  local name=$1 list=$2 map=$3 dst=$BLIND/$1 c k pair new orig p ext
  lx_require "$list" "$BIDS_ALL/dataset_description.json"
  lx_assert_inside "$dst"
  [ -s "$list" ] || { log "(dry-run) $list not present yet; subset $name skipped"; return 0; }
  [ -e "$dst" ] && [ "$DRY" != 1 ] && die "$dst exists"
  lx_run mkdir -p "$dst"
  lx_run cp "$BIDS_ALL/dataset_description.json" "$dst/"
  lx_run ln -s "../$(basename "$BIDS_ALL")/derivatives" "$dst/derivatives"
  while IFS=$'\t' read -r c k; do
    [ -n "$c" ] || continue
    lx_run mkdir -p "$dst/sub-$c/ses-$k/anat"
    for pair in $map; do
      new=${pair%%:*}; orig=${pair##*:}
      for p in mag phase; do for ext in nii.gz json; do
        lx_run ln -s "../../../../$(basename "$BIDS_ALL")/sub-$c/ses-$k/anat/sub-${c}_ses-${k}_echo-${orig}_part-${p}_MEGRE.$ext" \
                     "$dst/sub-$c/ses-$k/anat/sub-${c}_ses-${k}_echo-${new}_part-${p}_MEGRE.$ext"
      done; done
    done
  done < "$list"
  log "subset $name built from $(basename "$list")"
}

stage_subsets() {
  mk_subset bids_pilot   "$BLIND/pilot10.txt"          "1:1 2:2 3:3 4:4"
  mk_subset bids_lcurve  "$BLIND/lcurve2.txt"          "1:1 2:2 3:3 4:4"
  mk_subset bids_gate20  "$BLIND/gate20.txt"           "1:1 2:2 3:3 4:4"
  mk_subset bids_split13 "$BLIND/processing_order.txt" "1:1 2:3"
  mk_subset bids_split24 "$BLIND/processing_order.txt" "1:2 2:4"
}

mask3_one() {   # mask3_one CODE SES : Mask 1 (SynthStrip on RSS magnitude) [+ Mask 2] -> filled Mask 3
  local c=$1 k=$2 lab=sub-$1_ses-$2 src=$BLIND/bids_all/sub-$1/ses-$2/anat
  local w=$LX/work/mask3/$lab out=$BLIND/bids_all/derivatives/mask3/sub-$1/ses-$2/anat/${lab}_mask.nii.gz
  lx_assert_inside "$out"
  [ -s "$out" ] && [ "$DRY" != 1 ] && { log "$lab mask exists"; return 0; }
  lx_run mkdir -p "$w" "$(dirname "$out")"
  lx_run "$PY" - "$src" "$lab" "$w/rss.nii.gz" <<'PYEOF'
import sys, nibabel as nib, numpy as np
src, lab, out = sys.argv[1:4]
imgs = [nib.load(f"{src}/{lab}_echo-{n}_part-mag_MEGRE.nii.gz") for n in (1, 2, 3, 4)]
rss = np.sqrt(sum(np.asarray(i.get_fdata(dtype=np.float32)) ** 2 for i in imgs)).astype(np.float32)
nib.save(nib.Nifti1Image(rss, imgs[0].affine), out)
PYEOF
  lx_apx "$QSIPREP_SIF" mri_synthstrip -i "$w/rss.nii.gz" -m "$w/mask1.nii.gz"
  local q=""
  if [ "$MASK_MODE" = custom3 ]; then
    [ -s "$QMAP_CMD_FILE" ] || { [ "$DRY" = 1 ] || die "$QMAP_CMD_FILE missing (fill from $QMAP_HELP_TXT at A1)"; }
    local cmd="{QMAP-COMMAND-FROM-$QMAP_CMD_FILE}"
    [ -s "$QMAP_CMD_FILE" ] && cmd=$(head -1 "$QMAP_CMD_FILE")
    local n
    for n in 1 2 3 4; do
      cmd=${cmd//\{MAG$n\}/$src/${lab}_echo-${n}_part-mag_MEGRE.nii.gz}
      cmd=${cmd//\{PHASE$n\}/$src/${lab}_echo-${n}_part-phase_MEGRE.nii.gz}
    done
    cmd=${cmd//\{OUT\}/$w/quality.nii.gz}
    case "$cmd" in *VERIFY*|*\{*\}*) [ "$DRY" = 1 ] || die "quality-map command template unresolved: $cmd" ;; esac
    # shellcheck disable=SC2086
    lx_run $cmd
    q=$w/quality.nii.gz
  fi
  lx_run "$PY" - "$w/mask1.nii.gz" "${q:-none}" "$src/${lab}_echo-1_part-mag_MEGRE.nii.gz" "$out" "$HERE" "$N_PARTITIONS" "$EDGE_TRIM" <<'PYEOF'
import sys, nibabel as nib, numpy as np
from scipy import ndimage
m1p, qp, refp, out, here, npart, ntrim = sys.argv[1:8]
sys.path.insert(0, here)
import lx_common as L
ref = nib.load(refp)
m1 = np.asarray(nib.load(m1p).dataobj) > 0
if m1.shape != ref.shape[:3]:
    sys.exit("SynthStrip mask grid differs from the MEGRE grid")
if qp == "none":
    sys.exit("Mask 2 needs the phase-quality map (MASK_MODE custom3)")
q = np.asarray(nib.load(qp).get_fdata(dtype=np.float32))
m2 = q > L.otsu_threshold(q[m1])           # Mask 2: phase quality above its Otsu threshold within Mask 1 (Annex B.2)
m3 = ndimage.binary_fill_holes(m1 & m2)    # Mask 3 = fill-holes(Mask 1 AND Mask 2) ...
m3 = L.trim_partition_edges(m3, int(npart), int(ntrim))   # ... minus the first and last 2 slab partitions
hdr = ref.header.copy(); hdr.set_data_dtype(np.uint8); hdr.set_slope_inter(1, 0)
nib.save(nib.Nifti1Image(m3.astype(np.uint8), ref.affine, hdr), out)
PYEOF
}

stage_mask3() {
  [ "$MASK_MODE" = preset ] && { log "MASK_MODE=preset: QSMxT builds the masks; nothing to do"; return 0; }
  lx_require "$BLIND/processing_order.txt" "$QSIPREP_SIF"
  lx_stamp mask3
  export -f mask3_one lx_run lx_apx lx_assert_inside log die
  export MASK_MODE BLIND QMAP_CMD_FILE QMAP_HELP_TXT DRY HERE N_PARTITIONS EDGE_TRIM
  if [ "$DRY" = 1 ]; then
    head -1 "$BLIND/processing_order.txt" 2>/dev/null | while IFS=$'\t' read -r c k; do mask3_one "$c" "$k"; done
    echo "[DRY-RUN] ... repeated for every coded session in processing_order.txt with $LX_MAXJOBS parallel jobs"
  else
    tr '\t' ' ' < "$BLIND/processing_order.txt" | xargs -P "$LX_MAXJOBS" -L 1 bash -c 'mask3_one "$0" "$1"'
  fi
  lx_home_check mask3
}

all_variants() {
  local b i out=(primary split13 split24)
  for b in "${BF_LEVELS[@]}"; do for i in "${INV_LEVELS[@]}"; do
    [ "$b-$i" = vsharp-medi ] || out+=("mv_${b}_${i}")
  done; done
  echo "${out[@]}"
}

stage_negate() {   # sign-check branch (Annex B.4): phase * -1 for every session, headers unchanged
  [ -e "$BLIND/bids_all_neg" ] && [ "$DRY" != 1 ] && die "$BLIND/bids_all_neg exists"
  lx_stamp negate
  lx_run "$PY" - "$BLIND/bids_all" "$BLIND/bids_all_neg" "$HERE" <<'PYEOF'
import glob, os, shutil, sys
src, dst, here = sys.argv[1:4]
sys.path.insert(0, here)
import lx_common as L
import nibabel as nib, numpy as np
for f in sorted(glob.glob(os.path.join(src, "sub-*", "ses-*", "anat", "*_MEGRE.*"))):
    o = os.path.join(dst, os.path.relpath(f, src)); L.ensure_parent(o)
    if f.endswith(".json") or "_part-mag_" in f:
        shutil.copyfile(f, o)
    else:
        img = nib.load(f); nib.save(nib.Nifti1Image((-np.asarray(img.dataobj, np.float32)).astype(np.float32), None, header=img.header), o)
shutil.copyfile(os.path.join(src, "dataset_description.json"), os.path.join(dst, "dataset_description.json"))
os.symlink("../bids_all/derivatives", os.path.join(dst, "derivatives"))
print("negated phase tree written")
PYEOF
  [ "$DRY" = 1 ] || echo negated > "$SIGN_FILE"
  log "phase sign negated for every session; rebuild the subsets (--stage subsets on a clean blind/bids_* set) and recompute all runs"
}

stage_plan() {
  local v
  [ -x "$QSMXT_BIN" ] || [ "$DRY" = 1 ] || die "qsmxt not installed"
  lx_run "$QSMXT_BIN" validate "$BIDS_ALL"
  for v in pilot $(all_variants); do
    variant_spec "$v"; build_args; check_flags
    if [ "$DRY" = 1 ]; then lx_run "$QSMXT_BIN" run "$V_BIDS" "$V_OUT" "${ARGS[@]}" --dry
    else "$QSMXT_BIN" run "$V_BIDS" "$V_OUT" "${ARGS[@]}" --dry > "$LX/logs/qsmxt_plan_$v.log" 2>&1
         log "plan for $v -> logs/qsmxt_plan_$v.log"; fi
  done
}

stage_run() {
  local v lam
  case "$VARIANT" in
    pilot)      run_variant pilot ;;
    lcurve)     for lam in "${LAMBDA_GRID[@]}"; do run_variant "lcurve_$lam"; done ;;
    g2)         run_variant g2 ;;
    g2alt)      run_variant g2_mv_vsharp_rts; run_variant g2_mv_pdf_ilsqr ;;   # only if approval A2 was declined
    primary)    run_variant primary ;;
    split)      run_variant split13; run_variant split24 ;;
    multiverse) for v in $(all_variants); do case "$v" in mv_*) run_variant "$v" ;; esac; done ;;
    fallback)   run_variant mv_pdf_ilsqr ;;
    "")         die "--variant required" ;;
    *)          run_variant "$VARIANT" ;;
  esac
}

case "$STAGE" in
  install) stage_install ;;
  config)  stage_config ;;
  subsets) stage_subsets ;;
  mask3)   stage_mask3 ;;
  negate)  stage_negate ;;
  plan)    stage_plan ;;
  run)     stage_run ;;
  *) die "unknown stage $STAGE" ;;
esac
