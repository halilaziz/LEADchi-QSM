#!/usr/bin/env python3
"""qsm_plausibility.py - label-free QSM physical-plausibility rules on coded data (task P5b, 2026-09-25).

Implements v1.2 delta sec. 5 (Paper 0 lessons i and iv), frozen and hashed before any QSM run. Runs on server
(AFQ_env python) after run_qsmxt.sh and register_rois.sh; reads only the coded tree ($LX), never labels.

Stage `indices` (per coded session, one reconstruction variant):
  support S              Mask 3 AND the reconstruction's own final mask AND finite, non-zero chi; the run is invalid
                         for the session ("mask_mismatch") unless >= 99 % of that mask lies in Mask 3 and it covers
                         >= 90 % of Mask 3 eroded by 2 voxels (post-run mask assertion)
  whole-brain reference  two-step: ref0 = mean chi over S; implausible = S & |chi - ref0| > 1.0 ppm; ref = mean over
                         S minus implausible; chi_ref = chi - ref (registered text sec. 12.1)
  ROI rules              ROI failed if implausible voxels carry > 5 % of its in-mask weight, if < 95 % of its
                         weight lies in Mask 3 (< 90 % dentate), if its mean is outside the sec. 5.1 envelope,
                         or if it has no coverage; R2* ROI means outside [5, 100] s^-1 fail the R2* value
  session indices        f_imp (implausible fraction in Mask 3 eroded 3 voxels), chi SD in eroded WM, chi SD in
                         eroded lateral ventricles, residual background median(|local field|) in the 5-mm shell
                         inside S, posterior-fossa streaking SD(cerebellar WM)/SD(cerebral WM), Q_phase
  secondary reference    CSF only: lateral ventricles eroded 2 (fallback 1) voxels, choroid removed by R2* > median +
                         2 SD, >= 200 voxels and |mean CSF - WB| <= 50 ppb; otherwise 'none' (missing; no PLIC)
Stage `flags` (all sessions of the variant):
  auto grade 2   f_imp > 1 %, or >= 2 of the iron-rich ROIs (GP, SN-VTA, RN, DN_R) failed
  forced review  any session index > Q3 + 3 IQR of the pooled coded distribution
  (c') flag      SN-VTA failed (SN-only failure is flagged separately for the G3 fallback)
  robust z       per ROI over plausibility-passing sessions; |z| > 5 flagged and RETAINED; counts only visible

Outputs
  hidden (code-keyed continuous; read only by hashed scripts):  tables_coded/plaus_indices_<engine>_<variant>.csv
                                                                work/plaus/<engine>_<variant>/<sess>_implausible.nii.gz
  visible (pass/fail flags and pooled aggregates):             flags_coded/plaus_flags_<engine>_<variant>.csv
                                                                flags_coded/plaus_summary_<engine>_<variant>.json
--dry-run prints the planned inputs and outputs only.
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

CHECK_ROIS = ["SNVTA", "SNonly", "GP", "PUT", "PUT_R", "CAU", "RN", "DN_R", "DNIP_R", "THAL", "STR", "DGM"]
IRON_RICH = ["GP", "SNVTA", "RN", "DN_R"]
R2S_ROIS = ["SNVTA", "SNonly", "GP", "PUT", "CAU", "THAL", "RN", "DN_R", "STR"]
INDICES = ["f_imp", "wm_sd", "vent_sd", "resid_bg", "pf_streak", "q_phase"]


def paths_for(lx: str, engine: str, variant: str, code: str, ses: int, choices: dict, chi_pattern: str | None,
              mask_pattern: str | None = None):
    lab = L.sess_label(code, ses)
    vo = L.variant_outputs(lx, variant, code, ses, chi_pattern, mask_pattern)
    m3 = [os.path.join(lx, "blind", "bids_all", "derivatives", "mask3", f"sub-{code}", f"ses-{ses}", "anat",
                       f"{lab}_mask.nii.gz")]
    r2s = sorted(glob.glob(os.path.join(lx, "qsm", "primary", "derivatives", "qsmxt", f"sub-{code}", f"ses-{ses}",
                                        "anat", f"{lab}*_R2starmap.nii*")))
    ph = [os.path.join(lx, "blind", "bids_all", f"sub-{code}", f"ses-{ses}", "anat",
                       f"{lab}_echo-{n}_part-phase_MEGRE.nii.gz") for n in (1, 2, 3, 4)]
    return {"lab": lab, "chi": vo["chi"], "mask3": m3, "qmask": vo["mask"], "localfield": vo["localfield"], "r2s": r2s,
            "phase": ph, "reg": os.path.join(lx, "reg", f"sub-{code}", f"ses-{ses}")}


def one_session(job):
    import numpy as np
    from scipy import ndimage
    lx, engine, variant, code, ses, choices, cfg, chi_pattern, mask_pattern, with_q, dry = job
    P = paths_for(lx, engine, variant, code, ses, choices, chi_pattern, mask_pattern)
    L.forbid_raw_inputs(P["chi"] + P["mask3"] + P["phase"] + P["qmask"])
    row = {"code": code, "ses_order": ses, "variant": variant, "engine": engine, "status": "ok"}
    if dry:
        return row, None
    if len(P["chi"]) != 1 or not os.path.exists(P["mask3"][0]):
        row["status"] = "missing_chi_or_mask"
        return row, None
    if engine == "qsmxt" and len(P["qmask"]) != 1:
        row["status"] = "missing_reconstruction_mask"
        return row, None
    if not os.path.exists(os.path.join(P["reg"], "route.json")):
        row["status"] = "not_registered"
        return row, None
    cimg, chi = L.load_nii(P["chi"][0])
    mimg, m3 = L.load_nii(P["mask3"][0])
    if not L.same_grid(cimg, mimg):
        row["status"] = "grid_mismatch"
        return row, None
    qm = None
    if P["qmask"]:
        qimg, qm = L.load_nii(P["qmask"][0])
        if not L.same_grid(qimg, cimg):
            row["status"] = "grid_mismatch_mask"
            return row, None
        row.update(L.mask_assertion(qm, m3))
        if not row["mask_ok"]:
            row["status"] = "mask_mismatch"      # QSMxT did not use Mask 3: invalid for this run (Annex B.2)
            return row, None
    R = L.build_rois(P["reg"], cfg)
    if not L.same_grid(cimg, R["img"]):
        row["status"] = "grid_mismatch_reg"
        return row, None
    S = L.support_mask(m3, chi, qm)
    ref, imp, ref0 = L.two_step_reference(chi, S, cfg["implausible_voxel_abs_ppm"])
    chi_ref_ppb = (chi - ref) * 1000.0
    row["ref0_ppb"], row["ref_ppb"] = ref0 * 1000.0, ref * 1000.0
    row["support_frac_m3"] = float(S.sum() / max(1, (m3 > 0.5).sum()))
    er3 = L.erode(S, 3)
    row["f_imp"] = float(imp[er3].mean()) if er3.any() else float("nan")
    rois = R["rois"]
    sd = lambda sel: float(np.std(chi_ref_ppb[sel], ddof=1)) if sel.sum() > 1 else float("nan")
    row["wm_sd"] = sd(rois["WM"]["sup"] & S & ~imp)
    row["vent_sd"] = sd(rois["LV1"]["sup"] & S & ~imp)
    cwm_sd = sd(rois["CWM"]["sup"] & S & ~imp)
    row["pf_streak"] = cwm_sd / row["wm_sd"] if row["wm_sd"] and row["wm_sd"] == row["wm_sd"] else float("nan")
    row["resid_bg"] = float("nan")
    if P["localfield"]:
        limg, lfd = L.load_nii(P["localfield"][0])
        if L.same_grid(limg, cimg):
            dist = ndimage.distance_transform_edt(S, sampling=cimg.header.get_zooms()[:3])
            shell = S & (dist <= 5.0) & np.isfinite(lfd)
            row["resid_bg"] = float(np.median(np.abs(lfd[shell]))) if shell.any() else float("nan")   # median(|f|)
    row["q_phase"] = float("nan")
    if with_q:
        phases = [L.load_nii(p)[1] for p in P["phase"]]
        row["q_phase"] = L.q_phase(phases, m3 > 0.5)
    r2s = None
    if P["r2s"]:
        rimg, r2s = L.load_nii(P["r2s"][0])
        r2s = np.clip(r2s, *cfg["r2s_voxel_clip_s-1"]) if L.same_grid(rimg, cimg) else None
    env = cfg["envelopes_ppb"]
    for name in CHECK_ROIS:
        st = L.roi_stats(chi_ref_ppb, S, imp, rois[name])
        thr = cfg["roi_fail_inmask"].get(name, cfg["roi_fail_inmask"]["default"])
        reasons = []
        if st["sumw"] <= 0 or st["nvox"] == 0:
            reasons.append("no_coverage")
        else:
            if st["fimp"] > cfg["roi_fail_fimp"]:
                reasons.append("fimp")
            if st["inmask"] < thr:
                reasons.append("inmask")
            lo, hi = env[name]
            if not (lo <= st["mean"] <= hi):
                reasons.append("envelope")
        row[f"mean_{name}"] = st["mean"]
        row[f"inmask_{name}"] = st["inmask"]
        row[f"fimp_{name}"] = st["fimp"]
        row[f"fail_{name}"] = int(bool(reasons))
        row[f"why_{name}"] = ";".join(reasons)
    for name in R2S_ROIS:
        v = float("nan")
        if r2s is not None:
            valid_r = (m3 > 0.5) & np.isfinite(r2s)
            v = L.roi_stats(r2s, valid_r, None, rois[name])["mean"]
        lo, hi = cfg["r2s_envelope_s-1"]
        row[f"r2s_{name}"] = v
        row[f"fail_r2s_{name}"] = int(not (v == v and lo <= v <= hi))
    typ, csf_mean, n = L.csf_reference(chi_ref_ppb, S, imp, R["aseg"], r2s, cfg)
    row["csf_minus_wb_ppb"], row["csf_nvox"] = csf_mean, n
    row["ref2_type"] = typ
    row["ref2_mean_ppb"] = csf_mean if typ != "none" else float("nan")
    row.update(L.tissue_fractions(R["aseg"], S, imp, cfg))
    return row, (imp.astype("uint8"), cimg.affine, cimg.header)


def stage_indices(args, cfg, choices):
    import nibabel as nib
    lx = L.lx_root()
    sess = L.coded_sessions(lx)
    tag = f"{args.engine}_{args.variant}"
    out_tab = os.path.join(lx, "tables_coded", f"plaus_indices_{tag}.csv")
    out_dir = os.path.join(lx, "work", "plaus", tag)
    with_q = args.variant == "primary" and args.engine == "qsmxt"
    L.say(f"stage indices: {len(sess)} coded sessions, variant {tag}, Q_phase {'yes' if with_q else 'no'}", args.dry_run)
    L.say(f"outputs: {out_tab}; implausible masks in {out_dir}", args.dry_run)
    if args.dry_run:
        ex = paths_for(lx, args.engine, args.variant, sess[0][0], sess[0][1], choices, args.chi_pattern) if sess else {}
        L.say(f"example inputs (first coded session): {json.dumps(ex, indent=1)}", True)
        return
    if os.path.exists(out_tab):
        sys.exit(f"{out_tab} exists; plausibility indices are computed once per variant")
    rr = L.RunRecord(__file__, False)
    jobs = [(lx, args.engine, args.variant, c, s, choices, cfg, args.chi_pattern, args.mask_pattern, with_q, False)
            for c, s in sess]
    rows = []
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for row, res in ex.map(one_session, jobs):
            rows.append(row)
            if res is not None:
                imp, aff, hdr0 = res
                p = os.path.join(out_dir, f"{L.sess_label(row['code'], row['ses_order'])}_implausible.nii.gz")
                L.ensure_parent(p)
                hdr = hdr0.copy()
                hdr.set_data_dtype("uint8")
                hdr.set_slope_inter(1, 0)
                nib.save(nib.Nifti1Image(imp, aff, hdr), p)
    cols = ["code", "ses_order", "engine", "variant", "status", "ref0_ppb", "ref_ppb", "support_frac_m3",
            "mask_inside_m3", "mask_covers_m3core", "mask_ok"] + INDICES + \
           ["csf_minus_wb_ppb", "csf_nvox", "ref2_type", "ref2_mean_ppb", "ref_frac_wm", "ref_frac_gm", "ref_frac_csf"]
    for n in CHECK_ROIS:
        cols += [f"mean_{n}", f"inmask_{n}", f"fimp_{n}", f"fail_{n}", f"why_{n}"]
    for n in R2S_ROIS:
        cols += [f"r2s_{n}", f"fail_r2s_{n}"]
    rows.sort(key=lambda r: (r["code"], r["ses_order"]))
    L.write_csv(out_tab, rows, cols)
    rr.add_output(out_tab)
    rr.note(f"sessions={len(rows)} ok={sum(r['status'] == 'ok' for r in rows)}")
    rr.finish()
    L.say(f"wrote {out_tab} ({len(rows)} rows; values hidden by commitment)")


def q3iqr(vals, k):
    import numpy as np
    v = np.asarray([x for x in vals if x == x and x not in (float("inf"), float("-inf"))], float)
    if v.size < 4:
        return float("nan")
    q1, q3 = np.percentile(v, [25, 75])
    return float(q3 + k * (q3 - q1))


def stage_flags(args, cfg):
    import numpy as np
    lx = L.lx_root()
    tag = f"{args.engine}_{args.variant}"
    tab = os.path.join(lx, "tables_coded", f"plaus_indices_{tag}.csv")
    out = os.path.join(lx, "flags_coded", f"plaus_flags_{tag}.csv")
    summ = os.path.join(lx, "flags_coded", f"plaus_summary_{tag}.json")
    L.say(f"stage flags: {tab} -> {out}, {summ}", args.dry_run)
    if args.dry_run:
        return
    rows = L.read_csv(tab)
    num = lambda r, k: float(r[k]) if r.get(k) not in (None, "", "NA") else float("nan")
    ok = [r for r in rows if r["status"] == "ok"]
    k = cfg["session_review_iqr_k"]
    thr = {i: q3iqr([num(r, i) for r in ok], k) for i in INDICES}
    qp_primary = {}
    if not any(num(r, "q_phase") == num(r, "q_phase") for r in ok):
        ptab = os.path.join(lx, "tables_coded", f"plaus_indices_{args.engine}_primary.csv")
        if os.path.exists(ptab):   # Q_phase is session-level: take it from the primary run for the review rule
            qp_primary = {(r["code"], r["ses_order"]): num(r, "q_phase") for r in L.read_csv(ptab)}
            thr["q_phase"] = q3iqr(list(qp_primary.values()), k)
    flags, rz_counts = [], {}
    for r in rows:
        f = {"code": r["code"], "ses_order": int(r["ses_order"]), "status": r["status"]}
        if r["status"] != "ok":
            f.update({"auto_grade2": 1, "auto_reason": r["status"], "forced_review": 1, "review_reasons": "",
                      "snvta_fail": 1, "snonly_fail": 1})
            flags.append(f)
            continue
        auto = []
        if num(r, "f_imp") > cfg["session_auto_fail_fimp"]:
            auto.append("f_imp>1%")
        n_iron = sum(int(r[f"fail_{n}"]) for n in IRON_RICH)
        if n_iron >= 2:
            auto.append(f"{n_iron}_iron_rich_ROIs_failed")
        rev = []
        for i in INDICES:
            v = num(r, i) if i != "q_phase" or not qp_primary else qp_primary.get((r["code"], r["ses_order"]),
                                                                                  float("nan"))
            if thr[i] == thr[i] and v == v and v > thr[i]:
                rev.append(i)
        f.update({"auto_grade2": int(bool(auto)), "auto_reason": ";".join(auto),
                  "forced_review": int(bool(rev) or bool(auto)), "review_reasons": ";".join(rev),
                  "snvta_fail": int(r["fail_SNVTA"]), "snonly_fail": int(r["fail_SNonly"])})
        for n in CHECK_ROIS:
            f[f"fail_{n}"] = int(r[f"fail_{n}"])
        for n in R2S_ROIS:
            f[f"fail_r2s_{n}"] = int(r[f"fail_r2s_{n}"])
        flags.append(f)
    # robust z over plausibility-passing sessions (flag retained; counts only)
    for n in CHECK_ROIS:
        x = np.array([num(r, f"mean_{n}") for r in ok if r[f"fail_{n}"] == "0"], float)
        x = x[np.isfinite(x)]
        if x.size < 5:
            rz_counts[n] = {"n": int(x.size), "n_absz_gt5": None}
            continue
        med = np.median(x)
        mad = 1.4826 * np.median(np.abs(x - med))
        z = (x - med) / mad if mad > 0 else np.zeros_like(x)
        rz_counts[n] = {"n": int(x.size), "n_absz_gt5": int((np.abs(z) > cfg["robust_z_flag"]).sum()),
                        "sd_over_robust_sd": float(np.std(x, ddof=1) / mad) if mad > 0 else None}
    cols = ["code", "ses_order", "status", "auto_grade2", "auto_reason", "forced_review", "review_reasons",
            "snvta_fail", "snonly_fail"] + [f"fail_{n}" for n in CHECK_ROIS] + [f"fail_r2s_{n}" for n in R2S_ROIS]
    flags.sort(key=lambda r: (r["code"], r["ses_order"]))
    L.write_csv(out, flags, cols)
    summary = {
        "variant": tag, "n_sessions": len(rows), "n_ok": len(ok),
        "n_auto_grade2": sum(f["auto_grade2"] for f in flags), "n_forced_review": sum(f["forced_review"] for f in flags),
        "n_snvta_fail": sum(f["snvta_fail"] for f in flags), "n_snonly_fail": sum(f["snonly_fail"] for f in flags),
        "roi_fail_counts": {n: sum(int(r.get(f"fail_{n}", 0)) for r in ok) for n in CHECK_ROIS},
        "review_thresholds_q3_plus_3iqr": thr,
        "ref2_type_counts": {t: sum(r["ref2_type"] == t for r in ok) for t in ("CSF", "CSF_erode1", "none")},
        "status_counts": {st: sum(r["status"] == st for r in rows) for st in sorted({r["status"] for r in rows})},
        "robust_z": rz_counts, "utc": L.utcnow(),
    }
    L.write_json(summ, summary)
    rr = L.RunRecord(__file__, False)
    rr.add_input(tab)
    rr.add_output(out)
    rr.add_output(summ)
    rr.finish()
    print(json.dumps({k: summary[k] for k in ("n_sessions", "n_ok", "n_auto_grade2", "n_forced_review",
                                              "n_snvta_fail")}, indent=1))


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--stage", choices=["indices", "flags"], required=True)
    a.add_argument("--variant", required=True, help="primary | split13 | split24 | mv_<bf>_<inv> | pilot")
    a.add_argument("--engine", default="qsmxt", choices=["qsmxt", "sepia"])
    a.add_argument("--chi-pattern", default=None,
                   help="glob for non-QSMxT chi maps, e.g. {lx}/qsm/sepia/{lab}_Chimap.nii.gz (G2)")
    a.add_argument("--mask-pattern", default=None, help="glob for a non-QSMxT reconstruction mask (optional)")
    a.add_argument("--jobs", type=int, default=int(os.environ.get("LX_MAXJOBS", "8")))
    a.add_argument("--dry-run", action="store_true")
    args = a.parse_args()
    cfg = L.load_config()
    if args.stage == "indices":
        stage_indices(args, cfg, L.pipeline_choices() if not args.dry_run or
                      os.path.exists(os.path.join(L.lx_root(), "prereg", "pipeline_choices.json"))
                      else {"mask_source": "custom_mask3"})
    else:
        stage_flags(args, cfg)


if __name__ == "__main__":
    main()
