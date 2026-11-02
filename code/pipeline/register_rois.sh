#!/usr/bin/env bash
# register_rois.sh - T1w-only, chi-independent registration and atlas warps for LEAD-chi (task P5b, 2026-09-25).
# Written and hashed before the OSF registration; NOT executed before approvals A4/A5/A6. Runs on server only.
#
# Registered chain (protocol v1.1 sec. 12.2; E3 findings A10-A11):
#   MEGRE -> same-session T1w : N4-corrected echo-1 magnitude (registration only), ANTs 2.5.3 rigid (6 dof), Mattes
#                               MI, fixed-image mask = dilated SynthSeg brain.
#   T1w  -> MNI152NLin2009cAsym: REGISTERED ROUTE R for every session (registration Annex B.6): own
#          antsRegistrationSyN.sh (-t s) of the SynthSeg-masked T1w to the TemplateFlow 2009cAsym brain. The QSIPrep
#          transforms are not used (no inverse exists, and the coded copies carry per-session translation offsets).
#   Atlases: CIT168 v1.1 (2009cAsym) directly; SUIT dentate and JHU PLIC (FSL MNI152 = NLin6Asym) through the
#            label-free template-to-template SyN (E3 A11 option b). All transforms are composed into ONE
#            antsApplyTransforms call per atlas (single interpolation); probabilities Linear, labels GenericLabel.
#            chi is never resampled.
#   G5 inputs: TemplateFlow-2009c SynthSeg warped to the T1w grid (T1w->MNI Dice), SynthSeg of the coded echo-1
#            magnitude and the T1w SynthSeg on that grid with the MEGRE field of view (MEGRE<->T1w Dice).
#   Coverage: the MEGRE field of view on the T1w grid (brain-in-FOV fraction) and in 2009cAsym; stage `harmref`
#            (after every session is registered) builds the MNI region covered in >= 95 % of sessions (label-free) and
#            warps it back to each MEGRE grid (harm_megre.nii.gz; coverage-harmonised reference, sec. 24.4).
#
# Usage: register_rois.sh --stage templates|session|all|harmref [--code C --ses K] [--dry-run]
#   templates  once: template brains, template SynthSeg, NLin6->2009c SyN, atlas subsets (label-free)
#   session    one coded session -> $LX/reg/sub-C/ses-K/
#   all        every coded session in blind/processing_order.txt, <= LX_MAXJOBS parallel jobs
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=env.sh
source "$HERE/env.sh"
CFG=$HERE/atlas_config.json

STAGE="" CODE="" SES="" DRY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --stage) STAGE=$2; shift 2 ;;
    --code) CODE=$2; shift 2 ;;
    --ses) SES=$2; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) sed -n '2,26p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[ "$DRY" = 1 ] && export LX_DRY=1
log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

cfg() {   # cfg dotted.key -> value; relative paths resolved against $LX, ${VAR} expanded
  "$SYSPY" -c '
import json, os, sys
v = json.load(open(sys.argv[1]))
for k in sys.argv[2].split("."):
    v = v[int(k)] if isinstance(v, list) else v[k]
if isinstance(v, str):
    v = os.path.expandvars(v)
    print(v if os.path.isabs(v) or not sys.argv[3] else os.path.join(sys.argv[3], v))
else:
    print(json.dumps(v))
' "$CFG" "$1" "${2-$LX}"
}

ROIS=$LX/rois
T09=$(cfg templates.MNI152NLin2009cAsym.t1w);  M09=$(cfg templates.MNI152NLin2009cAsym.brain_mask)
T06=$(cfg templates.MNI152NLin6Asym.t1w);      M06=$(cfg templates.MNI152NLin6Asym.brain_mask)
B09=$ROIS/templates/tpl09_brain.nii.gz;        B06=$ROIS/templates/tpl06_brain.nii.gz
S09=$ROIS/templates/tpl09_synthseg.nii.gz
N6=$(cfg template_to_template.prefix)
CIT_SUB=$(cfg cit168.subset_file); DN_SUB=$(cfg dentate.subset_file); PLIC_SUB=$(cfg plic.subset_file)

ants_check() {
  [ "$DRY" = 1 ] && return 0
  antsRegistration --version 2>&1 | grep -q "2\.5\.3" || die "antsRegistration is not ANTs 2.5.3 (ANTSPATH=$ANTSPATH)"
}

brain_of() {   # brain_of T1 MASK OUT : template brain (mask if present and non-empty, else SynthStrip)
  local t=$1 m=$2 o=$3
  if [ -s "$m" ]; then lx_run ImageMath 3 "$o" m "$t" "$m"
  else lx_apx "$QSIPREP_SIF" mri_synthstrip -i "$t" -o "$o" -m "${o%.nii.gz}_mask.nii.gz"; fi
}

stage_templates() {   # approvals: A5 (template/atlas copies present), A6 (container), A4 (compute ~1 CPU-h)
  lx_require "$T09" "$T06" "$(cfg cit168.file)" "$QSIPREP_SIF"
  ants_check
  lx_stamp templates
  lx_run mkdir -p "$ROIS/templates" "$ROIS/xfm" "$(dirname "$CIT_SUB")" "$(dirname "$DN_SUB")" "$(dirname "$PLIC_SUB")"
  brain_of "$T09" "$M09" "$B09"
  brain_of "$T06" "$M06" "$B06"
  lx_apx "$QSIPREP_SIF" "$SYNTHSEG_BIN" --i "$T09" --o "$S09" --robust --threads "$LX_THREADS" --cpu
  lx_run antsRegistrationSyN.sh -d 3 -f "$B09" -m "$B06" -o "$N6" -t s -n "$LX_THREADS"
  # atlas subsets (label-free; checks volume counts and the pinned CIT168 SHA-256)
  lx_run "$PY" - "$CFG" "$LX" <<'PYEOF'
import hashlib, json, os, sys
import nibabel as nib, numpy as np
cfg, lx = json.load(open(sys.argv[1])), sys.argv[2]
P = lambda p: (lambda q: q if os.path.isabs(q) else os.path.join(lx, q))(os.path.expandvars(p))
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
c = cfg["cit168"]; f = P(c["file"])
if sha(f) != c["sha256"]:
    sys.exit(f"CIT168 SHA-256 mismatch for {f}")
img = nib.load(f)
if img.ndim != 4 or img.shape[3] != c["n_volumes"]:
    sys.exit(f"CIT168: expected 4D with {c['n_volumes']} volumes, got {img.shape}")
d = np.asarray(img.dataobj, dtype=np.float32)
sub = np.stack([d[..., c["volume_index_0based"][k]] for k in c["subset_order"]], axis=-1)
nib.save(nib.Nifti1Image(sub, img.affine), P(c["subset_file"]))
dn = cfg["dentate"]; fd = P(dn["file"])
if os.path.exists(fd):
    im = nib.load(fd); a = np.asarray(im.dataobj, dtype=np.float32)
    if a.ndim != 4:
        sys.exit(f"dentate probseg {fd} is not 4D; check the layout (atlas_config VERIFY)")
    sub = np.stack([a[..., dn["volume_index_0based"][k]] for k in dn["subset_order"]], axis=-1)
    print(f"dentate atlas {dn['atlas']} probseg sha256 {sha(fd)}")
else:
    fb = dn["dseg_fallback"]["file"]
    print(f"DEVIATION: dentate probseg absent; deterministic fallback {fb} (sha256 {sha(fb)})")
    im = nib.load(fb); lab = np.asarray(im.dataobj)
    sub = np.stack([(lab == dn["label_1based"][k]).astype(np.float32) for k in dn["subset_order"]], axis=-1)
nib.save(nib.Nifti1Image(sub.astype(np.float32), im.affine), P(dn["subset_file"]))
pl = cfg["plic"]; fp = P(pl["file"]); im = nib.load(fp); lab = np.asarray(im.dataobj)
m = np.isin(lab, list(pl["labels"].values())).astype(np.uint8)
if m.sum() == 0:
    sys.exit("PLIC labels not found in the JHU label image (atlas_config VERIFY)")
nib.save(nib.Nifti1Image(m, im.affine), P(pl["subset_file"]))
print(f"JHU labels sha256 {sha(fp)}; PLIC voxels {int(m.sum())}")
PYEOF
  lx_home_check templates
}

flag_of() {   # flag_of CODE SES COLUMN : read one label-free column of blind/session_flags.csv
  "$SYSPY" -c '
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if r["code"] == sys.argv[2] and r["ses_order"] == sys.argv[3]:
        print(r[sys.argv[4]]); break
' "$LX/blind/session_flags.csv" "$1" "$2" "$3"
}

stage_session() {
  [ -n "$CODE" ] && [ -n "$SES" ] || die "--code and --ses required"
  local lab=sub-${CODE}_ses-${SES}
  local A=$LX/blind/anat_all/sub-$CODE/ses-$SES/anat M=$LX/blind/bids_all/sub-$CODE/ses-$SES/anat
  local O=$LX/reg/sub-$CODE/ses-$SES
  local E1=$M/${lab}_echo-1_part-mag_MEGRE.nii.gz
  lx_assert_inside "$O"
  local route src frame
  route=$("$SYSPY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["xfm_route"])' "$LX/blind/setup_manifest.json" 2>/dev/null || echo "?")
  src=$(flag_of "$CODE" "$SES" anat_src 2>/dev/null || echo "?")
  if [ "$DRY" = 1 ]; then
    if [ "$route" = "?" ]; then route=R; fi
    if [ "$src" = "?" ]; then src=qsiprep; fi
  fi
  [ "$src" = none ] && { log "$lab: no T1w; not registered (first-usable criterion b)"; return 0; }
  [ "$route" = R ] || die "registered route is R for every session; setup_manifest says $route"
  if [ -e "$O/route.json" ] && [ "$DRY" != 1 ]; then log "$lab done"; return 0; fi
  lx_stamp "reg_$lab"
  lx_run mkdir -p "$O"
  local T1 ASEG
  if [ "$src" = qsiprep ]; then
    T1=$A/${lab}_space-ACPC_desc-preproc_T1w.nii.gz; ASEG=$A/${lab}_space-ACPC_desc-aseg_dseg.nii.gz; frame=acpc
  else
    T1=$A/${lab}_T1w.nii.gz; ASEG=$O/t1_synthseg.nii.gz; frame=t1w_raw
    lx_apx "$QSIPREP_SIF" "$SYNTHSEG_BIN" --i "$T1" --o "$ASEG" --threads "$LX_THREADS" --cpu
  fi
  lx_require "$E1" "$T1" "$ASEG"
  # --- T1w brain mask on the T1w grid (SynthSeg > 0, dilated 2, holes filled) and brain image
  lx_run antsApplyTransforms -d 3 -i "$ASEG" -r "$T1" -o "$O/aseg_t1grid.nii.gz" -n GenericLabel
  lx_run ThresholdImage 3 "$O/aseg_t1grid.nii.gz" "$O/brainmask.nii.gz" 1 100000 1 0
  lx_run ImageMath 3 "$O/brainmask.nii.gz" MD "$O/brainmask.nii.gz" 2
  lx_run ImageMath 3 "$O/brainmask.nii.gz" FillHoles "$O/brainmask.nii.gz"
  lx_run ImageMath 3 "$O/t1_brain.nii.gz" m "$T1" "$O/brainmask.nii.gz"
  # --- MEGRE -> T1w rigid (registration image only)
  lx_run N4BiasFieldCorrection -d 3 -i "$E1" -o "$O/megre_e1_n4.nii.gz" -s 2 -b "[200]" -c "[50x50x50x50,1e-7]"
  lx_run antsRegistration -d 3 --float 0 -o "[$O/megre2t1_]" \
      -r "[$T1,$O/megre_e1_n4.nii.gz,1]" -t "Rigid[0.1]" \
      -m "MI[$T1,$O/megre_e1_n4.nii.gz,1,32,Regular,0.25]" -c "[1000x500x250,1e-6,10]" -f 4x2x1 -s 2x1x0vox \
      -x "[$O/brainmask.nii.gz,NULL]" -v 0
  local R=$O/megre2t1_0GenericAffine.mat
  # --- T1w -> template chain (points in the T1w frame -> points in 2009cAsym), per route
  local -a CH_T1
  lx_run antsRegistrationSyN.sh -d 3 -f "$B09" -m "$O/t1_brain.nii.gz" -o "$O/t1_to_2009c_" -t s -n "$LX_THREADS"
  CH_T1=(-t "[$O/t1_to_2009c_0GenericAffine.mat,1]" -t "$O/t1_to_2009c_1InverseWarp.nii.gz")
  local -a CH_M=(-t "[$R,1]" "${CH_T1[@]}")
  local -a N6C=(-t "${N6}1Warp.nii.gz" -t "${N6}0GenericAffine.mat")
  # --- atlas warps (single interpolation each)
  lx_run antsApplyTransforms -d 3 -e 3 -i "$CIT_SUB" -r "$E1" -o "$O/cit_megre.nii.gz" -n Linear "${CH_M[@]}"
  lx_run antsApplyTransforms -d 3 -e 3 -i "$CIT_SUB" -r "$T1" -o "$O/cit_t1.nii.gz" -n Linear "${CH_T1[@]}"
  lx_run antsApplyTransforms -d 3 -e 3 -i "$DN_SUB" -r "$E1" -o "$O/dn_megre.nii.gz" -n Linear "${CH_M[@]}" "${N6C[@]}"
  lx_run antsApplyTransforms -d 3 -e 3 -i "$DN_SUB" -r "$T1" -o "$O/dn_t1.nii.gz" -n Linear "${CH_T1[@]}" "${N6C[@]}"
  lx_run antsApplyTransforms -d 3 -i "$PLIC_SUB" -r "$E1" -o "$O/plic_megre.nii.gz" -n GenericLabel "${CH_M[@]}" "${N6C[@]}"
  lx_run antsApplyTransforms -d 3 -i "$ASEG" -r "$E1" -o "$O/aseg_megre.nii.gz" -n GenericLabel -t "[$R,1]"
  lx_run antsApplyTransforms -d 3 -i "$T1" -r "$E1" -o "$O/t1_megre.nii.gz" -n Linear -t "[$R,1]"
  # --- G5 inputs
  lx_run antsApplyTransforms -d 3 -i "$S09" -r "$O/aseg_t1grid.nii.gz" -o "$O/tplseg_t1.nii.gz" -n GenericLabel "${CH_T1[@]}"
  lx_apx "$QSIPREP_SIF" "$SYNTHSEG_BIN" --i "$E1" --o "$O/megre_synthseg.nii.gz" --robust --threads "$LX_THREADS" --cpu
  lx_run antsApplyTransforms -d 3 -i "$ASEG" -r "$O/megre_synthseg.nii.gz" -o "$O/aseg_on_megreseg.nii.gz" -n GenericLabel -t "[$R,1]"
  lx_run ThresholdImage 3 "$E1" "$O/fov_megre.nii.gz" -1e30 1e30 1 0
  lx_run antsApplyTransforms -d 3 -i "$O/fov_megre.nii.gz" -r "$O/megre_synthseg.nii.gz" -o "$O/fov_on_megreseg.nii.gz" -n NearestNeighbor
  # --- coverage inputs: MEGRE field of view on the T1w grid and in 2009cAsym (forward chain: MEGRE -> T1w -> MNI)
  lx_run antsApplyTransforms -d 3 -i "$O/fov_megre.nii.gz" -r "$T1" -o "$O/fov_t1.nii.gz" -n NearestNeighbor -t "$R"
  lx_run antsApplyTransforms -d 3 -i "$O/fov_megre.nii.gz" -r "$B09" -o "$O/fov_2009c.nii.gz" -n NearestNeighbor \
      -t "$O/t1_to_2009c_1Warp.nii.gz" -t "$O/t1_to_2009c_0GenericAffine.mat" -t "$R"
  # --- provenance
  if [ "$DRY" = 1 ]; then
    echo "[DRY-RUN] write $O/route.json (route=$route, anat_src=$src, frame=$frame)"
  else
    "$SYSPY" - "$O/route.json" "$route" "$src" "$frame" "$(antsRegistration --version 2>&1 | head -1)" <<'PYEOF'
import datetime, json, sys
p, route, src, frame, ants = sys.argv[1:6]
json.dump({"route": route, "anat_src": src, "orient_frame": frame, "ants": ants,
           "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}, open(p, "w"), indent=2)
PYEOF
  fi
  lx_home_check "reg_$lab"
}

stage_all() {
  lx_require "$LX/blind/processing_order.txt" "$CIT_SUB" "$DN_SUB" "$PLIC_SUB" "$S09" "${N6}1Warp.nii.gz"
  ants_check
  if [ "$DRY" = 1 ]; then
    local c k
    IFS=$'\t' read -r c k < "$LX/blind/processing_order.txt" 2>/dev/null || { c=CODE; k=1; }
    CODE=$c SES=$k stage_session
    echo "[DRY-RUN] ... repeated for every coded session with $LX_MAXJOBS parallel jobs of $LX_THREADS threads"
    return 0
  fi
  tr '\t' ' ' < "$LX/blind/processing_order.txt" | \
    xargs -P "$LX_MAXJOBS" -L 1 bash -c '"$0" --stage session --code "$1" --ses "$2" > "$LX/logs/reg_sub-$1_ses-$2.log" 2>&1 || echo "FAILED sub-$1 ses-$2"' "$HERE/register_rois.sh"
}

stage_harmref() {   # label-free: MNI region inside the MEGRE FOV in >= 95 % of registered sessions -> each MEGRE grid
  lx_require "$LX/blind/processing_order.txt" "$B09"
  lx_stamp harmref
  lx_run "$PY" - "$LX" "$ROIS/templates/harm_region_2009c.nii.gz" <<'PYEOF'
import glob, os, sys
import nibabel as nib, numpy as np
lx, out = sys.argv[1:3]
fs = sorted(glob.glob(os.path.join(lx, "reg", "sub-*", "ses-*", "fov_2009c.nii.gz")))
if not fs:
    sys.exit("no fov_2009c.nii.gz found; run --stage all first")
ref = nib.load(fs[0]); acc = np.zeros(ref.shape[:3], np.float32)
for f in fs:
    acc += np.asarray(nib.load(f).dataobj) > 0
cov = acc / len(fs)
brain = np.asarray(nib.load(os.path.join(lx, "rois", "templates", "tpl09_brain.nii.gz")).dataobj) > 0
region = (cov >= 0.95) & brain
nib.save(nib.Nifti1Image(region.astype(np.uint8), ref.affine), out)
print(f"coverage-harmonised region: {int(region.sum())} template voxels from {len(fs)} sessions")
PYEOF
  local c k O E1
  while IFS=$'\t' read -r c k; do
    O=$LX/reg/sub-$c/ses-$k; E1=$LX/blind/bids_all/sub-$c/ses-$k/anat/sub-${c}_ses-${k}_echo-1_part-mag_MEGRE.nii.gz
    [ -e "$O/route.json" ] || continue
    lx_run antsApplyTransforms -d 3 -i "$ROIS/templates/harm_region_2009c.nii.gz" -r "$E1" -o "$O/harm_megre.nii.gz" \
        -n GenericLabel -t "[$O/megre2t1_0GenericAffine.mat,1]" -t "[$O/t1_to_2009c_0GenericAffine.mat,1]" -t "$O/t1_to_2009c_1InverseWarp.nii.gz"
  done < "$LX/blind/processing_order.txt"
  lx_home_check harmref
}

case "$STAGE" in
  templates) stage_templates ;;
  harmref)   stage_harmref ;;
  session)   ants_check; stage_session ;;
  all)       stage_all ;;
  *) sed -n '2,26p' "$0"; exit 2 ;;
esac
