#!/usr/bin/env python3
"""extract_rois.py - coded ROI tables for LEAD-chi (task P5b, 2026-09-25).

Written and hashed before the OSF registration; runs on server (AFQ_env python) after run_qsmxt.sh,
register_rois.sh and qsm_plausibility.py. Reads only the coded tree; writes code-keyed tables that the analyst
does not open (v1.2 sec. 2.5) and visible pass/fail flags. Column names follow P5's schema
(LEADchi/prereg/code/analysis/schema/roi_table_schema.csv); unlock_and_merge.py adds the label columns.

Stage `geom` (variant-independent, one row per coded session):
  roi_vol        native sum of w over CIT168 SN-VTA voxels with w >= 0.5, mm^3 (atlas only: no mask, no CSF buffer;
                 chi-independent; registered definition, sec. 12.2)
  roi_vol_SNonly as roi_vol for SNc + SNr (G3 fallback ROI)
  DN_R_vol       sum of p over right-dentate voxels with p >= 0.5, mm^3 (dentate atlas, D8)
  str_vol        volume of the eroded SynthSeg striatum used for chi_STR, mm^3
  slabcov_*      share of the ROI's T1w-grid weight that falls inside the MEGRE field of view (review < .95/.90)
  icv_seg_ml, v4_vol_ml   SynthSeg intracranial (all labels > 0) and 4th-ventricle volumes (coded side only)
  pitch/roll/yaw/b0_tilt, abs_pitch   head orientation relative to B0, ONE frame for all sessions: the rotation of
                 the composed MEGRE->T1w rigid and the affine stage of the route-R T1w->MNI152NLin2009cAsym registration
  motion_gsr     magnitude ghost-to-signal ratio on the echo-1 magnitude (phase-encode axis --pe-axis; label-free)
  fov_brain_frac share of the T1w SynthSeg brain inside the MEGRE field of view
  G5 Dice        T1w->MNI: SynthSeg brainstem + ventral DC vs template SynthSeg warped through the same chain;
                 MEGRE<->T1w: SynthSeg on echo-1 magnitude vs T1w SynthSeg, brainstem 16 / ventral DC 28+60 /
                 cerebellar WM 7+46, inside the MEGRE field of view; pass = every Dice >= 0.80 (gate G5)
  -> hidden tables_coded/session_geom.csv
Stage `g5`: visible flags_coded/g5_flags.csv (pass/fail, slab-review flag) and g5_summary.json, with the VDC threshold
  recalibrated on the pilot when prereg/g5_vdc_threshold.txt exists (pilot_checks.py --stage g5cal; amendment 1)
Stage `roi` (one reconstruction variant):
  chi_<ROI> (ppb) on the support S (Mask 3 AND the reconstruction mask AND finite non-zero chi), referenced to the
  two-step whole-brain reference (mean over S minus implausible voxels), implausible voxels excluded;
  chi_<ROI>_ref2 (CSF-only reference; NA when no valid CSF), r2s_<ROI> (primary R2*, clipped [0, 200] s^-1),
  in-support / implausible shares, chi_SNVTA_median, chi_SNVTA_impkeep, chi_SNVTA_harm (coverage-harmonised
  reference: S restricted to the MNI region covered in >= 95 % of sessions), reference-domain tissue fractions,
  SN-only columns for the G3 branch.
Stage `multiverse` (all 9 background-removal x inversion nodes at once):
  per session the INTERSECTION of the 9 supports; each node summarised on it (two-step reference, CSF-only
  reference, SN-VTA and SN-only weighted means and medians); a node value is NA where that node's own plausibility
  flags failed the ROI. -> hidden tables_coded/multiverse_<engine>.csv (G3 and the specification curve)
  ROIs other than SN-VTA and SN-only that failed a plausibility rule are written as NA (v1.2 sec. 5.3); SN-VTA
  and SN-only keep their values because failure removes the session through criterion (c') (build_usable.py).
  -> hidden tables_coded/roi_<engine>_<variant>.csv
--dry-run prints the plan only.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

CHI_ROIS = ["SNVTA", "SNonly", "VTAPBP", "RN", "STN", "GP", "PUT", "PUT_R", "CAU", "THAL", "STR", "DGM", "DN_R",
            "DN_L", "DNIP_R"]
KEEP_ON_FAIL = {"SNVTA", "SNonly"}
R2S_ROIS = ["SNVTA", "SNonly", "GP", "PUT", "CAU", "THAL", "RN", "DN_R", "STR"]


def reg_dir(lx, code, ses):
    return os.path.join(lx, "reg", f"sub-{code}", f"ses-{ses}")


# ------------------------------------------------------------------------------------------------ geom

def geom_one(job):
    import numpy as np
    lx, code, ses, cfg = job
    rd = reg_dir(lx, code, ses)
    row = {"code": code, "ses_order": ses, "registered": 0}
    if not os.path.exists(os.path.join(rd, "route.json")):
        return row
    row["registered"] = 1
    with open(os.path.join(rd, "route.json"), encoding="utf-8") as f:
        route = json.load(f)
    R = L.build_rois(rd, cfg)
    v = R["vox"]
    rois = R["rois"]
    row["roi_vol"] = float(rois["SNVTA"]["w"][rois["SNVTA"]["sup_atlas"]].sum() * v)
    row["roi_vol_SNonly"] = float(rois["SNonly"]["w"][rois["SNonly"]["sup_atlas"]].sum() * v)
    row["DN_R_vol"] = float(rois["DN_R"]["w"][rois["DN_R"]["sup_atlas"]].sum() * v)
    row["str_vol"] = float(rois["STR"]["sup"].sum() * v)
    # slab coverage: weight inside the MEGRE grid vs weight in the full T1w grid (unthresholded, mm^3)
    cimg, cit = L.load_nii(os.path.join(rd, "cit_megre.nii.gz"))
    ct_img, cit_t1 = L.load_nii(os.path.join(rd, "cit_t1.nii.gz"))
    dimg, dn = L.load_nii(os.path.join(rd, "dn_megre.nii.gz"))
    dt_img, dn_t1 = L.load_nii(os.path.join(rd, "dn_t1.nii.gz"))
    order = cfg["cit168"]["subset_order"]
    sn_idx = [order.index(k) for k in ("SNc", "SNr", "VTA", "PBP")]
    w_m = np.clip(cit[..., sn_idx].sum(-1), 0, 1).sum() * v
    w_t = np.clip(cit_t1[..., sn_idx].sum(-1), 0, 1).sum() * L.voxel_volume_mm3(ct_img)
    row["slabcov_SNVTA"] = float(w_m / w_t) if w_t > 0 else float("nan")
    di = cfg["dentate"]["subset_order"].index("DN_R")
    d_m = dn[..., di].sum() * L.voxel_volume_mm3(dimg)
    d_t = dn_t1[..., di].sum() * L.voxel_volume_mm3(dt_img)
    row["slabcov_DN_R"] = float(d_m / d_t) if d_t > 0 else float("nan")
    # SynthSeg volumes on the T1w grid
    a_img, a_t1 = L.load_labels(os.path.join(rd, "aseg_t1grid.nii.gz"))
    va = L.voxel_volume_mm3(a_img)
    row["icv_seg_ml"] = float((a_t1 > 0).sum() * va / 1000.0)
    row["v4_vol_ml"] = float((a_t1 == cfg["synthseg"]["V4"]).sum() * va / 1000.0)
    # orientation: one frame for every session (MNI152NLin2009cAsym via the route-R affine stage)
    chain = [L.read_itk_affine(os.path.join(rd, "t1_to_2009c_0GenericAffine.mat")),
             L.read_itk_affine(os.path.join(rd, "megre2t1_0GenericAffine.mat"))]
    row.update(L.orientation_angles(chain))
    row["orient_frame"] = "mni2009c_affine"
    # motion index (magnitude ghost-to-signal ratio) and brain-in-FOV fraction
    e1 = os.path.join(lx, "blind", "bids_all", f"sub-{code}", f"ses-{ses}", "anat", f"{L.sess_label(code, ses)}_echo-1_part-mag_MEGRE.nii.gz")
    if os.path.exists(e1):
        L.forbid_raw_inputs([e1])
        row["motion_gsr"] = L.ghost_to_signal(L.load_nii(e1)[1], int(cfg.get("pe_axis", 0)))
    fp = os.path.join(rd, "fov_t1.nii.gz")
    if os.path.exists(fp):
        _, fv = L.load_labels(fp)
        brain = a_t1 > 0
        row["fov_brain_frac"] = float((brain & (fv > 0)).sum() / max(1, brain.sum()))
    # G5
    ss = cfg["synthseg"]
    t_img, tpl = L.load_labels(os.path.join(rd, "tplseg_t1.nii.gz"))
    bsvdc = [ss["BS"], ss["L_VDC"], ss["R_VDC"]]
    row["dice_t1mni_bsvdc"] = L.dice(np.isin(a_t1, bsvdc), np.isin(tpl, bsvdc)) if L.same_grid(t_img, a_img) \
        else float("nan")
    s_img, sseg = L.load_labels(os.path.join(rd, "megre_synthseg.nii.gz"))
    _, aon = L.load_labels(os.path.join(rd, "aseg_on_megreseg.nii.gz"))
    _, fov = L.load_labels(os.path.join(rd, "fov_on_megreseg.nii.gz"))
    fov = fov > 0
    for nm, labs in (("bs", [ss["BS"]]), ("vdc", [ss["L_VDC"], ss["R_VDC"]]), ("cwm", [ss["L_CWM"], ss["R_CWM"]])):
        row[f"dice_megre_{nm}"] = L.dice(np.isin(sseg, labs) & fov, np.isin(aon, labs) & fov)
    # midbrain-level boundary distance (label-free QC aggregate for amendment 2): symmetric mean surface distance, mm,
    # between the brainstem + ventral-DC boundaries of the echo-1-magnitude SynthSeg and the T1w SynthSeg
    from scipy import ndimage
    a_ = np.isin(sseg, [ss["BS"], ss["L_VDC"], ss["R_VDC"]]) & fov
    b_ = np.isin(aon, [ss["BS"], ss["L_VDC"], ss["R_VDC"]]) & fov
    sa, sb_ = a_ & ~L.erode(a_, 1), b_ & ~L.erode(b_, 1)
    if sa.any() and sb_.any():
        z = s_img.header.get_zooms()[:3]
        da = ndimage.distance_transform_edt(~sa, sampling=z)
        db = ndimage.distance_transform_edt(~sb_, sampling=z)
        row["msd_bs_vdc_mm"] = float(0.5 * (da[sb_].mean() + db[sa].mean()))
    else:
        row["msd_bs_vdc_mm"] = float("nan")
    g5_apply(row, cfg["g5_dice_min"], cfg["g5_dice_min"])
    row["dice_min"] = min(row["dice_t1mni_bsvdc"], row["dice_megre_bs"], row["dice_megre_vdc"], row["dice_megre_cwm"])
    row["slab_review"] = int(not (row["slabcov_SNVTA"] >= 0.95 and row["slabcov_DN_R"] >= 0.90))
    return row


def g5_apply(row, dmin, vdc_min):
    """G5 per session: T1w->MNI brainstem+VDC Dice >= 0.80 and MEGRE<->T1w Dice >= 0.80 for brainstem and cerebellar
    WM and >= vdc_min for the ventral DC (vdc_min = 0.80 unless recalibrated on the pilot in amendment 1)."""
    ok = lambda x, t: x == x and x >= t
    row["g5_t1mni_pass"] = int(ok(row["dice_t1mni_bsvdc"], dmin))
    row["g5_megre_pass"] = int(ok(row["dice_megre_bs"], dmin) and ok(row["dice_megre_vdc"], vdc_min) and ok(row["dice_megre_cwm"], dmin))
    row["g5_pass"] = int(row["g5_t1mni_pass"] and row["g5_megre_pass"])
    row["slab_review"] = int(not (row.get("slabcov_SNVTA", 0) >= 0.95 and row.get("slabcov_DN_R", 0) >= 0.90))
    return row


GEOM_COLS = ["code", "ses_order", "registered", "roi_vol", "roi_vol_SNonly", "DN_R_vol", "str_vol", "slabcov_SNVTA",
             "slabcov_DN_R", "icv_seg_ml", "v4_vol_ml", "pitch_deg", "roll_deg", "yaw_deg", "b0_tilt_deg",
             "abs_pitch", "orient_frame", "motion_gsr", "fov_brain_frac", "dice_t1mni_bsvdc", "dice_megre_bs",
             "dice_megre_vdc", "dice_megre_cwm", "dice_min", "msd_bs_vdc_mm", "g5_t1mni_pass", "g5_megre_pass", "g5_pass",
             "slab_review"]


def stage_geom(args, cfg):
    lx = L.lx_root()
    sess = L.coded_sessions(lx)
    out_h = os.path.join(lx, "tables_coded", "session_geom.csv")
    L.say(f"stage geom: {len(sess)} coded sessions -> {out_h} (hidden)", args.dry_run)
    if args.dry_run:
        return
    if os.path.exists(out_h):
        sys.exit(f"{out_h} exists")
    rr = L.RunRecord(__file__, False)
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        rows = list(ex.map(geom_one, [(lx, c, s, cfg) for c, s in sess]))
    rows.sort(key=lambda r: (r["code"], r["ses_order"]))
    L.write_csv(out_h, rows, GEOM_COLS)
    rr.add_output(out_h)
    rr.finish()
    L.say(f"wrote {out_h} ({len(rows)} sessions; values hidden). Next: pilot_checks.py --stage g5cal, then --stage g5")


def stage_g5(args, cfg):
    """Visible G5 flags from the hidden geometry table with the registered (possibly pilot-recalibrated) thresholds."""
    lx = L.lx_root()
    gtab = os.path.join(lx, "tables_coded", "session_geom.csv")
    out_v = os.path.join(lx, "flags_coded", "g5_flags.csv")
    calp = os.path.join(lx, "prereg", "g5_vdc_threshold.txt")
    vdc = float(open(calp).read().strip()) if os.path.exists(calp) else cfg["g5_dice_min"]
    L.say(f"stage g5: {gtab} -> {out_v} (VDC threshold {vdc})", args.dry_run)
    if args.dry_run:
        return
    if os.path.exists(out_v):
        sys.exit(f"{out_v} exists")
    num = lambda v: float(v) if v not in ("", "NA", None) else float("nan")
    rows = []
    for r in L.read_csv(gtab):
        row = {k: (num(r[k]) if k.startswith(("dice_", "slabcov_")) else r[k]) for k in r}
        row["registered"] = int(r["registered"])
        if row["registered"]:
            g5_apply(row, cfg["g5_dice_min"], vdc)
        rows.append(row)
    L.write_csv(out_v, rows, ["code", "ses_order", "registered", "g5_t1mni_pass", "g5_megre_pass", "g5_pass", "slab_review"])
    summ = {"n": len(rows), "registered": sum(r["registered"] for r in rows), "g5_pass": sum(int(r.get("g5_pass", 0) or 0) for r in rows),
            "slab_review": sum(int(r.get("slab_review", 0) or 0) for r in rows), "vdc_threshold": vdc}
    summ["g5_fail_share"] = 1 - summ["g5_pass"] / max(1, summ["registered"])
    summ["g5_revise_chain"] = summ["g5_fail_share"] > 0.10     # > 10 % -> one blind revision from the registered menu
    L.write_json(os.path.join(lx, "flags_coded", "g5_summary.json"), summ)
    print(json.dumps(summ, indent=1))


# ------------------------------------------------------------------------------------------------ roi

def roi_one(job):
    import numpy as np
    lx, engine, variant, code, ses, cfg, choices, pattern, mpattern, plaus_row, geom_row, qphase = job
    lab = L.sess_label(code, ses)
    row = {"code": code, "ses_order": ses, "engine": engine, "variant": variant,
           "pipeline": "fallback_PDF_iLSQR" if variant == choices.get("fallback_variant") else
           ("primary" if variant == "primary" else variant)}
    row.update({k: geom_row.get(k) for k in ("roi_vol", "roi_vol_SNonly", "DN_R_vol", "str_vol", "abs_pitch",
                                             "pitch_deg", "roll_deg", "yaw_deg", "b0_tilt_deg", "orient_frame",
                                             "slabcov_SNVTA", "slabcov_DN_R", "icv_seg_ml", "v4_vol_ml", "g5_pass",
                                             "motion_gsr", "fov_brain_frac", "dice_min")})
    row["Q_phase"] = qphase
    vo = L.variant_outputs(lx, variant, code, ses, pattern, mpattern)
    rd = reg_dir(lx, code, ses)
    if plaus_row is None or plaus_row.get("status") != "ok" or len(vo["chi"]) != 1:
        row["status"] = (plaus_row or {}).get("status", "missing")
        return row
    row["status"] = "ok"
    mp = os.path.join(lx, "blind", "bids_all", "derivatives", "mask3", f"sub-{code}", f"ses-{ses}", "anat", f"{lab}_mask.nii.gz")
    L.forbid_raw_inputs([vo["chi"][0], mp] + vo["mask"])
    cimg, chi = L.load_nii(vo["chi"][0])
    _, m3 = L.load_nii(mp)
    qm = L.load_nii(vo["mask"][0])[1] if vo["mask"] else None
    R = L.build_rois(rd, cfg)
    if not L.same_grid(cimg, R["img"]):
        row["status"] = "grid_mismatch_reg"
        return row
    S = L.support_mask(m3, chi, qm)
    ref, imp, ref0 = L.two_step_reference(chi, S, cfg["implausible_voxel_abs_ppm"])   # identical to qsm_plausibility.py
    chi_ref = (chi - ref) * 1000.0
    ref2 = float(plaus_row["ref2_mean_ppb"]) if plaus_row.get("ref2_mean_ppb") not in (None, "", "NA") else float("nan")
    row["ref_ppb"] = ref * 1000.0
    row["ref2_type"] = plaus_row.get("ref2_type")
    row["ref2_minus_wb_ppb"] = ref2
    row.update(L.tissue_fractions(R["aseg"], S, imp, cfg))
    hp = os.path.join(rd, "harm_megre.nii.gz")
    harm_ref = float("nan")
    if os.path.exists(hp):
        _, hm = L.load_labels(hp)
        Sh = S & (hm > 0)
        harm_ref, _, _ = L.two_step_reference(chi, Sh, cfg["implausible_voxel_abs_ppm"])
    for name in CHI_ROIS:
        st = L.roi_stats(chi_ref, S, imp, R["rois"][name])
        failed = plaus_row.get(f"fail_{name}", "0") == "1"
        val = st["mean"] if (name in KEEP_ON_FAIL or not failed) else float("nan")
        row[f"chi_{name}"] = val
        row[f"chi_{name}_ref2"] = val - ref2 if (val == val and ref2 == ref2) else float("nan")
        row[f"chi_{name}_raw"] = st["mean"]
        row[f"inmask_{name}"] = st["inmask"]
        row[f"fimp_{name}"] = st["fimp"]
        row[f"fail_{name}"] = int(failed)
        if name in ("SNVTA", "SNonly", "DN_R"):
            row[f"chi_{name}_median"] = st["median"]
            row[f"chi_{name}_impkeep"] = st["mean_impkeep"]
            row[f"chi_{name}_harm"] = val + (ref - harm_ref) * 1000.0 if (val == val and harm_ref == harm_ref) else float("nan")
    r2p = sorted(glob.glob(os.path.join(lx, "qsm", "primary", "derivatives", "qsmxt", f"sub-{code}", f"ses-{ses}",
                                        "anat", f"{lab}*_R2starmap.nii*")))
    if r2p:
        rimg, r2s = L.load_nii(r2p[0])
        if L.same_grid(rimg, cimg):
            r2s = np.clip(r2s, *cfg["r2s_voxel_clip_s-1"])
            vr = (m3 > 0.5) & np.isfinite(r2s)
            for name in R2S_ROIS:
                v = L.roi_stats(r2s, vr, None, R["rois"][name])["mean"]
                row[f"r2s_{name}"] = v if plaus_row.get(f"fail_r2s_{name}", "0") != "1" else float("nan")
    row["chi_CSF_ref"] = plaus_row.get("csf_minus_wb_ppb")
    return row


def roi_columns():
    cols = ["code", "ses_order", "engine", "pipeline", "variant", "status",
            "chi_SNVTA", "chi_SNVTA_ref2", "ref2_type", "r2s_SNVTA", "roi_vol", "Q_phase",
            "chi_STR", "str_vol", "chi_GP", "chi_PUT", "chi_CAU", "chi_THAL", "chi_RN",
            "chi_DN_R", "chi_DGM", "chi_PUT_R", "DN_R_vol", "abs_pitch", "chi_DNIP_R",
            "fimp_SNVTA", "inmask_SNVTA", "fimp_DN_R", "inmask_DN_R",
            "chi_SNonly", "chi_SNonly_ref2", "r2s_SNonly", "roi_vol_SNonly", "fimp_SNonly", "inmask_SNonly",
            "chi_SNVTA_median", "chi_SNVTA_impkeep", "chi_SNonly_median", "chi_SNonly_impkeep", "chi_DN_R_median",
            "chi_DN_R_impkeep", "chi_SNVTA_harm", "chi_SNonly_harm", "chi_DN_R_harm", "chi_VTAPBP", "chi_STN",
            "chi_DN_L", "chi_STR_ref2", "chi_DN_R_ref2", "ref_ppb", "ref2_minus_wb_ppb", "chi_CSF_ref",
            "ref_frac_wm", "ref_frac_gm", "ref_frac_csf", "fov_brain_frac", "motion_gsr", "dice_min",
            "pitch_deg", "roll_deg", "yaw_deg", "b0_tilt_deg", "orient_frame", "slabcov_SNVTA", "slabcov_DN_R",
            "icv_seg_ml", "v4_vol_ml", "g5_pass"]
    cols += [f"r2s_{n}" for n in R2S_ROIS if f"r2s_{n}" not in cols]
    for n in CHI_ROIS:
        for pre in ("fimp_", "inmask_", "fail_"):
            if f"{pre}{n}" not in cols:
                cols.append(f"{pre}{n}")
        if f"chi_{n}_raw" not in cols:
            cols.append(f"chi_{n}_raw")
    return cols


def stage_roi(args, cfg):
    lx = L.lx_root()
    tag = f"{args.engine}_{args.variant}"
    out = os.path.join(lx, "tables_coded", f"roi_{tag}.csv")
    ptab = os.path.join(lx, "tables_coded", f"plaus_indices_{tag}.csv")
    qtab = os.path.join(lx, "tables_coded", "plaus_indices_qsmxt_primary.csv")
    gtab = os.path.join(lx, "tables_coded", "session_geom.csv")
    L.say(f"stage roi: {tag}; inputs {ptab}, {qtab} (Q_phase), {gtab}; output {out} (hidden)", args.dry_run)
    if args.dry_run:
        return
    for p in (ptab, qtab, gtab):
        if not os.path.exists(p):
            sys.exit(f"{p} missing (order: register_rois.sh -> extract_rois.py --stage geom -> "
                     "qsm_plausibility.py indices/flags -> extract_rois.py --stage roi)")
    if os.path.exists(out):
        sys.exit(f"{out} exists")
    choices = L.pipeline_choices()
    plaus = {(r["code"], int(r["ses_order"])): r for r in L.read_csv(ptab)}
    qph = {(r["code"], int(r["ses_order"])): (float(r["q_phase"]) if r["q_phase"] not in ("", "NA") else float("nan"))
           for r in L.read_csv(qtab)}
    geom = {(r["code"], int(r["ses_order"])): {k: (float(v) if k not in ("code", "orient_frame") and v not in ("", "NA")
                                                 else (None if v in ("", "NA") else v)) for k, v in r.items()}
            for r in L.read_csv(gtab)}
    rr = L.RunRecord(__file__, False)
    for p in (ptab, qtab, gtab):
        rr.add_input(p)
    jobs = [(lx, args.engine, args.variant, c, s, cfg, choices, args.chi_pattern, args.mask_pattern, plaus.get((c, s)),
             geom.get((c, s), {}), qph.get((c, s), float("nan"))) for c, s in L.coded_sessions(lx)]
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        rows = list(ex.map(roi_one, jobs))
    rows.sort(key=lambda r: (r["code"], r["ses_order"]))
    L.write_csv(out, rows, roi_columns())
    rr.add_output(out)
    rr.note(f"rows={len(rows)} ok={sum(r.get('status') == 'ok' for r in rows)}")
    rr.finish()
    L.say(f"wrote {out}: {len(rows)} rows, {sum(r.get('status') == 'ok' for r in rows)} ok (values not displayed)")


# ------------------------------------------------------------------------------------------------ multiverse

MV_NODES = [(b, i) for b in ("vsharp", "pdf", "lbv") for i in ("medi", "ilsqr", "rts")]


def node_variant(b, i):
    return "primary" if (b, i) == ("vsharp", "medi") else f"mv_{b}_{i}"


def mv_one(job):
    import numpy as np
    lx, engine, code, ses, cfg, plaus = job
    lab = L.sess_label(code, ses)
    out = []
    rd = reg_dir(lx, code, ses)
    mp = os.path.join(lx, "blind", "bids_all", "derivatives", "mask3", f"sub-{code}", f"ses-{ses}", "anat", f"{lab}_mask.nii.gz")
    if not os.path.exists(os.path.join(rd, "route.json")) or not os.path.exists(mp):
        return out
    _, m3 = L.load_nii(mp)
    maps, supports = {}, []
    for b, i in MV_NODES:
        v = node_variant(b, i)
        pr = plaus.get(v, {}).get((code, ses))
        vo = L.variant_outputs(lx, v, code, ses)
        if pr is None or pr.get("status") != "ok" or len(vo["chi"]) != 1 or len(vo["mask"]) != 1:
            continue
        L.forbid_raw_inputs(vo["chi"] + vo["mask"])
        _, chi = L.load_nii(vo["chi"][0])
        _, qm = L.load_nii(vo["mask"][0])
        maps[(b, i)] = (chi, pr)
        supports.append(L.support_mask(m3, chi, qm))
    if not maps:
        return out
    Sc = np.logical_and.reduce(supports)            # common support of the available nodes
    R = L.build_rois(rd, cfg)
    for (b, i), (chi, pr) in maps.items():
        ref, imp, _ = L.two_step_reference(chi, Sc, cfg["implausible_voxel_abs_ppm"])
        cr = (chi - ref) * 1000.0
        typ, csf, _ = L.csf_reference(cr, Sc, imp, R["aseg"], None, cfg)
        row = {"code": code, "ses_order": ses, "engine": engine, "bf": b, "inv": i, "n_nodes_available": len(maps),
               "csf_type": typ, "common_support_frac_m3": float(Sc.sum() / max(1, (m3 > 0.5).sum()))}
        for roi in ("SNVTA", "SNonly"):
            st = L.roi_stats(cr, Sc, imp, R["rois"][roi])
            bad = pr.get(f"fail_{roi}", "0") == "1"
            for k, val in (("wb_mean", st["mean"]), ("wb_median", st["median"])):
                row[f"chi_{roi}_{k}"] = float("nan") if bad else val
                row[f"chi_{roi}_{k.replace('wb', 'csf')}"] = float("nan") if (bad or typ == "none") else val - csf
        out.append(row)
    return out


MV_COLS = ["code", "ses_order", "engine", "bf", "inv", "n_nodes_available", "csf_type", "common_support_frac_m3"] + \
    [f"chi_{r}_{k}" for r in ("SNVTA", "SNonly") for k in ("wb_mean", "wb_median", "csf_mean", "csf_median")]


def stage_multiverse(args, cfg):
    lx = L.lx_root()
    out = os.path.join(lx, "tables_coded", f"multiverse_{args.engine}.csv")
    L.say(f"stage multiverse: 9 nodes, common support -> {out} (hidden)", args.dry_run)
    if args.dry_run:
        return
    if os.path.exists(out):
        sys.exit(f"{out} exists")
    plaus = {}
    for b, i in MV_NODES:
        v = node_variant(b, i)
        p = os.path.join(lx, "flags_coded", f"plaus_flags_{args.engine}_{v}.csv")
        if os.path.exists(p):
            plaus[v] = {(r["code"], int(r["ses_order"])): r for r in L.read_csv(p)}
    rr = L.RunRecord(__file__, False)
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        rows = [r for lst in ex.map(mv_one, [(lx, args.engine, c, s, cfg, plaus) for c, s in L.coded_sessions(lx)]) for r in lst]
    rows.sort(key=lambda r: (r["code"], r["ses_order"], r["bf"], r["inv"]))
    L.write_csv(out, rows, MV_COLS)
    rr.add_output(out)
    rr.note(f"rows={len(rows)}")
    rr.finish()
    L.say(f"wrote {out}: {len(rows)} session x node rows (values not displayed)")


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--stage", choices=["geom", "g5", "roi", "multiverse"], required=True)
    a.add_argument("--variant", default="primary")
    a.add_argument("--engine", default="qsmxt", choices=["qsmxt", "sepia"])
    a.add_argument("--chi-pattern", default=None, help="glob for non-QSMxT chi maps ({lx} {code} {ses} {lab})")
    a.add_argument("--mask-pattern", default=None, help="glob for a non-QSMxT reconstruction mask (optional)")
    a.add_argument("--jobs", type=int, default=int(os.environ.get("LX_MAXJOBS", "8")))
    a.add_argument("--dry-run", action="store_true")
    args = a.parse_args()
    cfg = L.load_config()
    {"geom": stage_geom, "g5": stage_g5, "roi": stage_roi, "multiverse": stage_multiverse}[args.stage](args, cfg)


if __name__ == "__main__":
    main()
