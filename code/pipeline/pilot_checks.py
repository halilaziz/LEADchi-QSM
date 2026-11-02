#!/usr/bin/env python3
"""pilot_checks.py - pilot sign check and the L-curve lambda rule (registration Annex B.4; amendment 1).

Label-free; reads coded pilot outputs only; run through audit_run.py. Output JSON is aggregate-only.

Stage `sign` (the 10 sessions of blind/pilot10.txt): pooled median, over the pilot sessions, of the whole-brain-
  referenced ROI means (the *_raw columns of extract_rois.py --variant pilot, i.e. before the envelope rules) of GP,
  SN-VTA, RN and right dentate.
    all four > 0  -> "pass"
    all four < 0  -> "flip": the phase sign is negated for ALL sessions (a new coded phase tree written by
                     blind_setup.py --negate-phase from the raw data is NOT allowed; instead run_qsmxt.sh --stage
                     negate writes phase * -1 copies inside blind/), every run is recomputed and the pilot chi is
                     discarded; logged in amendment 1
    mixed signs   -> "indeterminate": no flip; one registered diagnostic revision is allowed (menu: (1) re-check the
                     phase scaling from the header audit, (2) re-check the Mask 3 coverage of the four ROIs, (3)
                     re-check the reference support share); if the signs are still mixed after it, the study
                     follows the methods-report path (Q24)
Stage `g5cal` (the pilot sessions): G5 ventral-DC Dice threshold calibration (see stage_g5cal).
Stage `lcurve` (the first 2 pilot sessions; blind/lcurve2.txt): for each session, engine and lambda of the grid,
  rho = || M (D * chi - f_local) ||_2 (ppm; D the dipole kernel in k-space with the B0 direction taken from the
  affine; f_local the reconstruction's local field in ppm; M its final mask) and eta = || M grad(chi) ||_1
  (forward differences, ppm per voxel). Corner rule: x = log10 rho, y = log10 eta over the 4 grid points, each
  min-max normalised to [0, 1]; the session's choice is the grid point nearest the origin; if the largest
  perpendicular distance of the two interior points from the chord joining the end points is < 0.10 the session's
  choice is the engine default; the 2 sessions are combined by the geometric mean of their choices, snapped to the
  nearest grid value on the log scale, ties toward the engine default. Writes prereg/amendment1_lambda_<engine>.txt.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

GRIDS = {"qsmxt": [4e-3, 2e-3, 1e-3, 5e-4], "sepia": [250.0, 500.0, 1000.0, 2000.0]}
DEFAULTS = {"qsmxt": 1e-3, "sepia": 1000.0}
SIGN_ROIS = ["chi_GP_raw", "chi_SNVTA_raw", "chi_RN_raw", "chi_DN_R_raw"]


def stage_sign(args):
    import numpy as np
    lx = L.lx_root()
    pilot = [tuple(ln.split("\t")) for ln in open(os.path.join(lx, "blind", "pilot10.txt")).read().split("\n") if ln.strip()]
    pk = {(c, int(s)) for c, s in pilot}
    tab = [r for r in L.read_csv(os.path.join(lx, "tables_coded", f"roi_{args.engine}_pilot.csv"))
           if (r["code"], int(r["ses_order"])) in pk and r.get("status") == "ok"]
    med = {}
    for c in SIGN_ROIS:
        v = [float(r[c]) for r in tab if r.get(c) not in (None, "", "NA")]
        med[c] = float(np.median(v)) if v else None
    signs = [m for m in med.values() if m is not None]
    state = "pass" if signs and len(signs) == 4 and all(m > 0 for m in signs) else \
        ("flip" if signs and len(signs) == 4 and all(m < 0 for m in signs) else "indeterminate")
    out = {"utc": L.utcnow(), "engine": args.engine, "n_pilot_ok": len(tab), "pooled_medians_ppb": med, "sign_check": state}
    L.write_json(args.out, out)
    print(json.dumps(out, indent=1))


def stage_g5cal(args):
    """G5 ventral-DC threshold calibration on the pilot (registration Q24 G5): if more than 20 % of the pilot sessions
    fall below 0.80 on the MEGRE<->T1w ventral-DC Dice alone, the VDC threshold becomes the 10th percentile of the
    pilot's VDC Dice (never below 0.70); otherwise it stays 0.80. Written to prereg/g5_vdc_threshold.txt."""
    import numpy as np
    lx = L.lx_root()
    pilot = {(c, int(s)) for c, s in (ln.split("\t") for ln in open(os.path.join(lx, "blind", "pilot10.txt")).read().split("\n") if ln.strip())}
    v = [float(r["dice_megre_vdc"]) for r in L.read_csv(os.path.join(lx, "tables_coded", "session_geom.csv"))
         if (r["code"], int(r["ses_order"])) in pilot and r.get("dice_megre_vdc") not in ("", "NA")]
    share_low = float(np.mean([x < 0.80 for x in v])) if v else None
    thr = 0.80 if (share_low is None or share_low <= 0.20) else float(max(0.70, np.percentile(v, 10)))
    out = {"utc": L.utcnow(), "n_pilot": len(v), "share_below_0.80": share_low, "vdc_threshold": thr}
    L.write_json(args.out, out)
    tp = os.path.join(lx, "prereg", "g5_vdc_threshold.txt")
    if os.path.exists(tp):
        sys.exit(f"{tp} exists (fixed once, in amendment 1)")
    with open(tp, "w") as f:
        f.write(f"{thr:.4f}\n")
    print(json.dumps(out, indent=1))


def dipole_kernel(shape, zooms, b0_vox):
    import numpy as np
    k = np.meshgrid(*[np.fft.fftfreq(n, d=z) for n, z in zip(shape, zooms)], indexing="ij")
    k2 = k[0] ** 2 + k[1] ** 2 + k[2] ** 2
    kb = k[0] * b0_vox[0] + k[1] * b0_vox[1] + k[2] * b0_vox[2]
    with np.errstate(invalid="ignore", divide="ignore"):
        D = 1.0 / 3.0 - np.where(k2 > 0, kb ** 2 / k2, 0.0)
    D[0, 0, 0] = 0.0
    return D


def b0_in_voxel_frame(aff):
    """Unit B0 (world +z) expressed along the voxel axes (direction cosines of the affine, zooms removed)."""
    import numpy as np
    R = aff[:3, :3] / np.linalg.norm(aff[:3, :3], axis=0)
    b = np.linalg.solve(R, np.array([0.0, 0.0, 1.0]))
    return b / np.linalg.norm(b)


def norms(chi, field, mask, aff, zooms):
    import numpy as np
    D = dipole_kernel(chi.shape, zooms, b0_in_voxel_frame(aff))
    pred = np.real(np.fft.ifftn(D * np.fft.fftn(np.where(mask, chi, 0.0))))
    rho = float(np.sqrt(np.sum(((pred - field) * mask) ** 2)))
    g = [np.diff(chi, axis=a, append=np.take(chi, [-1], axis=a)) for a in range(3)]
    eta = float(sum(np.sum(np.abs(gi) * mask) for gi in g))
    return rho, eta


def corner(rhos, etas, grid, default):
    x = [math.log10(r) for r in rhos]
    y = [math.log10(e) for e in etas]
    nx = [(v - min(x)) / (max(x) - min(x)) if max(x) > min(x) else 0.0 for v in x]
    ny = [(v - min(y)) / (max(y) - min(y)) if max(y) > min(y) else 0.0 for v in y]
    # perpendicular distance of the interior points from the chord (end points 0 and 3)
    x0, y0, x1, y1 = nx[0], ny[0], nx[-1], ny[-1]
    den = math.hypot(x1 - x0, y1 - y0) or 1.0
    dmax = max(abs((y1 - y0) * nx[i] - (x1 - x0) * ny[i] + x1 * y0 - y1 * x0) / den for i in (1, 2))
    if dmax < 0.10:
        return default, {"flat": True, "dmax": dmax, "x": nx, "y": ny}
    d0 = [math.hypot(a, b) for a, b in zip(nx, ny)]
    return grid[int(min(range(len(grid)), key=lambda i: d0[i]))], {"flat": False, "dmax": dmax, "x": nx, "y": ny}


def snap(values, grid, default):
    gm = math.exp(sum(math.log(v) for v in values) / len(values))
    lg = [abs(math.log(g) - math.log(gm)) for g in grid]
    best = min(lg)
    cands = [g for g, d in zip(grid, lg) if abs(d - best) < 1e-12]
    return default if default in cands else cands[0]


def stage_lcurve(args):
    import numpy as np
    lx = L.lx_root()
    grid, default = GRIDS[args.engine], DEFAULTS[args.engine]
    sess = [tuple(ln.split("\t")) for ln in open(os.path.join(lx, "blind", "lcurve2.txt")).read().split("\n") if ln.strip()]
    per = []
    for i, (c, s) in enumerate(sess, start=1):
        rhos, etas = [], []
        for lam in grid:
            v = f"lcurve_{lam:g}" if args.engine == "qsmxt" else f"lcurve_sepia_{lam:g}"
            vo = L.variant_outputs(lx, v, c, int(s), args.chi_pattern and args.chi_pattern.replace("{lam}", f"{lam:g}"),
                                   args.mask_pattern and args.mask_pattern.replace("{lam}", f"{lam:g}"))
            fld = vo["localfield"] if not args.field_pattern else \
                sorted(__import__("glob").glob(args.field_pattern.replace("{lam}", f"{lam:g}").format(lx=lx, code=c, ses=s, lab=L.sess_label(c, int(s)), variant=v)))
            if len(vo["chi"]) != 1 or len(fld) != 1 or len(vo["mask"]) != 1:
                sys.exit(f"L-curve session {i}, lambda {lam:g}: chi / local field / mask not found uniquely")
            L.forbid_raw_inputs(vo["chi"] + fld + vo["mask"])
            cimg, chi = L.load_nii(vo["chi"][0], "float64")
            _, f = L.load_nii(fld[0], "float64")
            _, m = L.load_nii(vo["mask"][0])
            r, e = norms(chi, f, m > 0.5, cimg.affine, cimg.header.get_zooms()[:3])
            rhos.append(r)
            etas.append(e)
        choice, det = corner(rhos, etas, grid, default)
        per.append({"session_index": i, "rho": rhos, "eta": etas, "choice": choice, **det})
    lam = snap([p["choice"] for p in per], grid, default)
    out = {"utc": L.utcnow(), "engine": args.engine, "grid": grid, "default": default, "sessions": per, "lambda": lam,
           "rule": "Annex B.4 corner rule; geometric mean of the 2 session choices snapped to the grid (ties -> default)"}
    L.write_json(args.out, out)
    lf = os.path.join(lx, "prereg", f"amendment1_lambda_{args.engine}.txt")
    L.assert_inside(lf)
    if os.path.exists(lf):
        sys.exit(f"{lf} exists (lambda is fixed once, in amendment 1)")
    with open(lf, "w") as fh:
        fh.write(f"{lam:g}\n")
    print(json.dumps({"lambda": lam, "choices": [p["choice"] for p in per]}, indent=1))


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--stage", choices=["sign", "lcurve", "g5cal"], required=True)
    a.add_argument("--engine", choices=["qsmxt", "sepia"], default="qsmxt")
    a.add_argument("--chi-pattern", default=None, help="SEPIA chi glob ({lam} {lx} {code} {ses} {lab} {variant})")
    a.add_argument("--mask-pattern", default=None)
    a.add_argument("--field-pattern", default=None, help="SEPIA local-field glob in ppm")
    a.add_argument("--out", default=None)
    args = a.parse_args()
    lx = L.lx_root()
    args.out = args.out or os.path.join(lx, "gates", f"pilot_{args.stage}_{args.engine}.json")
    L.assert_inside(args.out)
    if os.path.exists(args.out):
        sys.exit(f"{args.out} exists (computed once)")
    {"sign": stage_sign, "lcurve": stage_lcurve, "g5cal": stage_g5cal}[args.stage](args)


if __name__ == "__main__":
    main()
