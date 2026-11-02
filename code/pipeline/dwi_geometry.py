#!/usr/bin/env python3
"""dwi_geometry.py - real-ID DWI node geometry, pSCP_R segment values and n_min for LEAD-chi S3 (registration 12.3, B.7).

Runs on server (python with nibabel, scipy and dipy, e.g. the qsirecon container's python as for Paper 0) on EXISTING
DWI outputs only (pyAFQ 2.1 bundles, QSIRecon AMICO NODDI maps, QSIPrep anatomy and b0 reference, PennLINC DWI QC).
DWI is non-blind by construction; nothing here reads MEGRE data, chi, or any group label. Outputs are real-ID tables
of node geometry and node/segment DWI values (no volumes, no chi).

Stage `geometry` (before posting; fixes n_min):
  per DWI session with an SCP_R bundle:
    - SCP_R, MCP_R and VTAuncross_R streamlines (RAS mm, ACPC space), each streamline oriented so that its first
      point is nearer the right-dentate centroid, resampled to 100 nodes; node centroids = mean of the resampled points
    - right dentate (SUIT atl-Anatom label 30, p >= 0.5) voxel centres and the pyAFQ SCP_R intermediate-ROI centroid
      mapped from template space to the session's ACPC space with ANTs antsApplyTransformsToPoints (NLin6 -> 2009c
      template-to-template SyN, then the QSIPrep from-ACPC_to-MNI152NLin2009cAsym transform, which maps template
      points to ACPC points); a built-in check requires the mapped template brainstem centroid to fall inside the
      session's SynthSeg brainstem (label 16), otherwise the session is 'transform_check_failed'
    - node rules on SCP_R nodes 0-30: (i) upper limit min(30, n_int - 1), n_int = node nearest the intermediate-ROI
      centroid; (ii) drop nodes whose centroid lies < 4 mm from the nearest mapped dentate point (3 mm plus the
      1-voxel dilation at the 1 mm atlas grid); (iii) drop nodes with ISOVF > 0.30
    - profiles (dipy afq_profile with Gaussian weights, as pyAFQ): ICVF, ISOVF, TF = 1 - ISOVF, ICVF_mod = ICVF x TF;
      node tw-ICVF = profile(ICVF_mod) / profile(TF); N4-corrected b0 (ANTs N4BiasFieldCorrection) divided by the
      median N4 b0 of the right cerebellar WM (SynthSeg 46, eroded 1 voxel, excluding voxels within 2 voxels of the
      mapped dentate)
    - segment values at each candidate n_min: twICVF_pSCP_R = mean node tw-ICVF over retained nodes; b0loc_pSCP_R;
      tract_len (mean SCP_R streamline length, mm); nSL (streamline count); twICVF_MCP_R_matched = mean MCP_R node
      tw-ICVF over MCP nodes 0-79 whose dentate distance is within +/-1 mm of a retained pSCP_R node's; whole-tract
      MCP_R tw-ICVF (nodes with ISOVF <= 0.30; PC5a); VTAuncross_R nodes 0-30 under the same rules;
      twICVF_pSCP_R_peak (profile shifted by s in [-5, 5] nodes to maximise its Pearson r with the pooled mean SCP_R
      tw-ICVF profile over nodes 0-30, then the same retained-node mean)
    - split-half streamline reliability: 20 random halvings (seed 20260924) of the SCP_R streamlines; pooled over
      sessions, the Spearman-Brown r of the two half-bundle tw-ICVF values is reported (aggregate only)
  n_min = the largest of {10, 8, 6} that keeps >= 90 % of the pre-QC S3 candidate sessions (--s3-candidates; 108
  sessions from metadata) after rules (i)-(iii). Written to dwi/n_min.json.
Stage `t1t2` (optional; S3 sensitivity): T2w rigidly registered to the ACPC T1w (antsRegistrationSyN.sh -t r); the
  ratio T1w/T2w calibrated so that the median in the lateral-ventricle CSF (SynthSeg 4/43, eroded 1) is 0 and in the
  cerebral WM (2/41, eroded 2) is 1; profiled along SCP_R; t1t2_pSCP_R = mean over the retained nodes.
Stage `table` (after the DWI bundle-QC freeze): dwi/dwi_s3_table.csv with s3_eligible = PennLINC scalar pass AND
  SCP_R bundle present AND coded bundle-QC pass AND >= n_min retained nodes, the registered S3 DWI columns and the
  round-1 AFQ-QC list flag; this table is frozen, hashed and merged by unlock_and_merge.py --dwi-table.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

N_NODES = 100
SEG_MAX = 30
DIST_MM = 4.0
ISOVF_MAX = 0.30
NMIN_GRID = (10, 8, 6)
KEEP_FRAC = 0.90
BUNDLES = {"SCP_R": "SCPR", "MCP_R": "MCPR", "VTAuncross_R": "VTAuncrossR"}   # pyAFQ desc- labels (VERIFY at install)


def cfg_default(lx):
    return {"pyafq": os.path.join(L.pl_root(), "pyAFQ"), "noddi": os.path.join(L.pl_root(), "QSIrecon_noddi"),
            "qsiprep": os.path.join(L.pl_root(), "QSIprep"), "bids": os.path.join(L.pl_root(), "ds007116"),
            "trk": "{pyafq}/{sub}/{ses}/**/{sub}_{ses}*desc-{label}_tractography.trk",
            "noddi_map": "{noddi}/**/{sub}_{ses}*model-noddi_param-{param}_dwimap.nii.gz",
            "b0": "{qsiprep}/{sub}/{ses}/dwi/{sub}_{ses}*space-ACPC*dwiref.nii.gz",
            "aseg": "{qsiprep}/{sub}/{ses}/anat/{sub}_{ses}_space-ACPC_desc-aseg_dseg.nii.gz",
            "t1w": "{qsiprep}/{sub}/{ses}/anat/{sub}_{ses}_space-ACPC_desc-preproc_T1w.nii.gz",
            "t2w": "{bids}/{sub}/{ses}/anat/{sub}_{ses}*T2w.nii.gz",
            "xfm_acpc2mni": "{qsiprep}/{sub}/{ses}/anat/{sub}_{ses}_from-ACPC_to-MNI152NLin2009cAsym_mode-image_xfm.h5",
            "suit": os.path.join(lx, "rois", "SUIT", "atl-Anatom_space-MNI_probseg.nii"),
            "suit_dn_r_index": 29,
            "nlin6_to_2009c": os.path.join(lx, "rois", "xfm", "nlin6_to_2009c_"),
            "scp_inter": os.path.join(lx, "rois", "afq_templates", "SCP_R_inter_prob.nii.gz"),
            "tpl_brainstem_2009c_mm": [0.0, -26.0, -28.0],   # upper pons, MNI152NLin2009cAsym (VERIFY on one session)
            "pennlinc_qc": os.path.join(lx, "prereg", "inputs", "EF_desc-diffusion_qc.csv")}


def one(pattern: str):
    m = sorted(glob.glob(pattern, recursive=True))
    return m[0] if len(m) == 1 else None


# ------------------------------------------------------------------------------------------ geometry helpers

def orient_and_resample(streamlines, anchor):
    """Orient each streamline so that its first point is nearer `anchor` (RAS mm) and resample to N_NODES points."""
    import numpy as np
    from dipy.tracking.streamline import set_number_of_points
    out = []
    for s in streamlines:
        s = np.asarray(s, float)
        if np.linalg.norm(s[-1] - anchor) < np.linalg.norm(s[0] - anchor):
            s = s[::-1]
        out.append(s)
    return set_number_of_points(out, N_NODES)


def node_centroids(res):
    import numpy as np
    return np.mean(np.stack([np.asarray(s) for s in res]), axis=0)


def nearest_dist(points, cloud):
    import numpy as np
    from scipy.spatial import cKDTree
    return cKDTree(cloud).query(points)[0] if len(cloud) else np.full(len(points), np.inf)


def profile(img_data, affine, res):
    """dipy afq_profile with Gaussian weights on streamlines already in voxel space of the map."""
    import numpy as np
    from dipy.stats.analysis import afq_profile, gaussian_weights
    from dipy.tracking.streamline import Streamlines, transform_streamlines
    vox = Streamlines(transform_streamlines(res, np.linalg.inv(affine)))
    w = gaussian_weights(vox, n_points=N_NODES)
    return np.asarray(afq_profile(img_data, vox, np.eye(4), weights=w, n_points=N_NODES), float)


def retained_nodes(dist, isovf, n_int, upper=SEG_MAX):
    """Registered node rules (i)-(iii) on nodes 0..30; returns a boolean vector over the 100 nodes."""
    import numpy as np
    keep = np.zeros(N_NODES, bool)
    top = min(upper, int(n_int) - 1) if n_int is not None else upper
    keep[: max(0, top) + 1] = True
    keep &= dist >= DIST_MM
    keep &= ~(isovf > ISOVF_MAX)
    return keep


def peak_shift(prof, template, maxs=5):
    import numpy as np
    best, bs = -np.inf, 0
    for s in range(-maxs, maxs + 1):
        idx = np.arange(0, SEG_MAX + 1) + s
        ok = (idx >= 0) & (idx < N_NODES)
        a, b = prof[idx[ok]], template[np.arange(0, SEG_MAX + 1)[ok]]
        if ok.sum() > 5 and np.isfinite(a).all() and np.std(a) > 0:
            r = np.corrcoef(a, b)[0, 1]
            if r > best:
                best, bs = r, s
    return bs


# ------------------------------------------------------------------------------------------ template points -> ACPC

def lps(p):
    return [-p[0], -p[1], p[2]]


def map_points(points_mm, transforms, identity=False):
    """Template RAS points -> session ACPC RAS points with antsApplyTransformsToPoints (points are passed in LPS).

    `transforms` is the ordered list of -t arguments (template-to-template inverse, then the QSIPrep ACPC->MNI h5,
    which maps template points to ACPC points). identity=True is for the synthetic test only."""
    import numpy as np
    pts = np.asarray(points_mm, float).reshape(-1, 3)
    if identity:
        if os.environ.get("LX_TEST") != "1":
            sys.exit("identity point mapping is a test-only mode (LX_TEST=1)")
        return pts
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as td:
        fi, fo = os.path.join(td, "in.csv"), os.path.join(td, "out.csv")
        with open(fi, "w") as f:
            f.write("x,y,z,t\n")
            for p in pts:
                q = lps(p)
                f.write(f"{q[0]},{q[1]},{q[2]},0\n")
        cmd = ["antsApplyTransformsToPoints", "-d", "3", "-i", fi, "-o", fo]
        for t in transforms:
            cmd += ["-t", t]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        with open(fo) as f:
            rows = list(csv.DictReader(f))
    return np.array([lps([float(r["x"]), float(r["y"]), float(r["z"])]) for r in rows])


def template_points(cfg):
    """Right-dentate voxel centres (SUIT NLin6, p >= 0.5) and the SCP_R intermediate-ROI centroid (template RAS mm)."""
    import nibabel as nib
    import numpy as np
    im = nib.load(cfg["suit"])
    p = np.asarray(im.dataobj, np.float32)[..., cfg["suit_dn_r_index"]]
    ijk = np.argwhere(p >= 0.5)
    dn = nib.affines.apply_affine(im.affine, ijk)
    si = nib.load(cfg["scp_inter"])
    q = np.asarray(si.dataobj, np.float32)
    ij = np.argwhere(q >= 0.5 * q.max())
    inter = nib.affines.apply_affine(si.affine, ij).mean(axis=0)
    return dn, inter


# ------------------------------------------------------------------------------------------ one session

def session_job(job):
    import nibabel as nib
    import numpy as np
    from dipy.io.streamline import load_tractogram
    from scipy import ndimage
    sub, ses, cfg, tpl, identity, n4 = job
    fm = dict(cfg, sub=sub, ses=ses)
    row = {"participant_id": sub, "session_id": ses, "status": "ok"}
    trk = {k: one(cfg["trk"].format(**fm, label=v)) for k, v in BUNDLES.items()}
    maps = {k: one(cfg["noddi_map"].format(**fm, param=k)) for k in ("icvf", "isovf")}
    b0p, asegp, xfm = one(cfg["b0"].format(**fm)), one(cfg["aseg"].format(**fm)), one(cfg["xfm_acpc2mni"].format(**fm))
    if not trk["SCP_R"]:
        row["status"] = "no_SCP_R_bundle"
        return row, None
    if not (maps["icvf"] and maps["isovf"] and asegp and (xfm or identity)):
        row["status"] = "missing_inputs"
        return row, None
    t6 = cfg["nlin6_to_2009c"]
    chain = [f"[{t6}0GenericAffine.mat,1]", f"{t6}1InverseWarp.nii.gz", xfm] if not identity else []
    dn_t, inter_t = tpl
    dn = map_points(dn_t, chain, identity)
    inter = map_points(inter_t, chain, identity)[0]
    aimg = nib.load(asegp)
    aseg = np.asarray(aimg.dataobj).astype(int)
    # transform-direction check: 7 template brainstem points (centre and +/-3 mm per axis) must map into the
    # session's SynthSeg brainstem (label 16) in >= 5 of 7 cases
    c = np.asarray(cfg["tpl_brainstem_2009c_mm"], float)
    probe = np.vstack([c] + [c + d for d in (np.eye(3) * 3.0)] + [c - d for d in (np.eye(3) * 3.0)])
    if identity:
        probe = nib.affines.apply_affine(aimg.affine, np.argwhere(aseg == 16).mean(axis=0))[None, :].repeat(7, 0)
        bsp = probe
    else:
        bsp = map_points(probe, [xfm], identity)
    vv = np.rint(nib.affines.apply_affine(np.linalg.inv(aimg.affine), bsp)).astype(int)
    n_in = sum(1 for v in vv if all(0 <= v[i] < aseg.shape[i] for i in range(3)) and aseg[tuple(v)] == 16)
    row["transform_check_inside"] = n_in
    if n_in < 5:
        row["status"] = "transform_check_failed"
        return row, None
    anchor = dn.mean(axis=0)
    icvf_img, isovf_img = nib.load(maps["icvf"]), nib.load(maps["isovf"])
    icvf = np.asarray(icvf_img.get_fdata(), float)
    isovf = np.asarray(isovf_img.get_fdata(), float)
    tf = 1.0 - isovf
    b0n = None
    if b0p:
        b0img = nib.load(b0p)
        b0 = np.asarray(b0img.get_fdata(), float)
        if n4 and not identity:
            with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as td:
                o = os.path.join(td, "b0n4.nii.gz")
                subprocess.run(["N4BiasFieldCorrection", "-d", "3", "-i", b0p, "-o", o], check=True, stdout=subprocess.DEVNULL)
                b0 = np.asarray(nib.load(o).get_fdata(), float)
        # right cerebellar WM on the b0 grid (nearest neighbour through the shared ACPC world), eroded 1, minus the
        # dentate dilated by 2 voxels
        ijk = np.indices(b0.shape).reshape(3, -1).T
        wmm = nib.affines.apply_affine(b0img.affine, ijk)
        av = np.rint(nib.affines.apply_affine(np.linalg.inv(aimg.affine), wmm)).astype(int)
        ok = np.all((av >= 0) & (av < np.array(aseg.shape)), axis=1)
        lab = np.zeros(len(ijk), int)
        lab[ok] = aseg[av[ok, 0], av[ok, 1], av[ok, 2]]
        cwm = ndimage.binary_erosion((lab == 46).reshape(b0.shape), iterations=1)
        dmask = np.zeros(b0.shape, bool)
        dv = np.rint(nib.affines.apply_affine(np.linalg.inv(b0img.affine), dn)).astype(int)
        dv = dv[np.all((dv >= 0) & (dv < np.array(b0.shape)), axis=1)]
        dmask[dv[:, 0], dv[:, 1], dv[:, 2]] = True
        cwm &= ~ndimage.binary_dilation(dmask, iterations=2)
        if cwm.sum() > 20:
            b0n = (b0 / float(np.median(b0[cwm])), b0img.affine)
    out = {}
    nodes = []
    for key, path in trk.items():
        if not path:
            continue
        sft = load_tractogram(path, "same", bbox_valid_check=False)
        sft.to_rasmm()
        sl = list(sft.streamlines)
        if not sl:
            continue
        res = orient_and_resample(sl, anchor)
        cen = node_centroids(res)
        dist = nearest_dist(cen, dn)
        p_icvf = profile(icvf, icvf_img.affine, res)
        p_iso = profile(isovf, isovf_img.affine, res)
        p_tf = profile(tf, isovf_img.affine, res)
        p_mod = profile(icvf * tf, icvf_img.affine, res)
        tw = p_mod / p_tf
        p_b0 = profile(b0n[0], b0n[1], res) if b0n is not None else np.full(N_NODES, np.nan)
        lens = [float(np.sum(np.linalg.norm(np.diff(np.asarray(s), axis=0), axis=1))) for s in sl]
        out[key] = {"res": res, "cen": cen, "dist": dist, "icvf": p_icvf, "isovf": p_iso, "tw": tw, "b0": p_b0,
                    "len": float(np.mean(lens)), "n": len(sl)}
        for i in range(N_NODES):
            nodes.append({"participant_id": sub, "session_id": ses, "bundle": key, "node": i, "cx": cen[i, 0],
                          "cy": cen[i, 1], "cz": cen[i, 2], "dist_dentate_mm": dist[i], "icvf": p_icvf[i],
                          "isovf": p_iso[i], "tf": 1 - p_iso[i], "twicvf": tw[i], "b0loc": p_b0[i]})
    S = out["SCP_R"]
    n_int = int(np.argmin(np.linalg.norm(S["cen"] - inter, axis=1)))
    keep = retained_nodes(S["dist"], S["isovf"], n_int)
    row.update({"n_int": n_int, "n_retained": int(keep.sum()), "tract_len": S["len"], "nSL": S["n"],
                "twICVF_pSCP_R": float(np.nanmean(S["tw"][keep])) if keep.any() else float("nan"),
                "b0loc_pSCP_R": float(np.nanmean(S["b0"][keep])) if keep.any() else float("nan")})
    for nd in nodes:
        if nd["bundle"] == "SCP_R":
            nd["retained"] = int(keep[nd["node"]])
    if "MCP_R" in out:
        M = out["MCP_R"]
        cand = np.arange(N_NODES) < 80
        match = np.zeros(N_NODES, bool)
        for d in S["dist"][keep]:
            match |= cand & (np.abs(M["dist"] - d) <= 1.0)
        row["twICVF_MCP_R_matched"] = float(np.nanmean(M["tw"][match])) if match.any() else float("nan")
        okw = ~(M["isovf"] > ISOVF_MAX)
        row["twICVF_MCP_R_whole"] = float(np.nanmean(M["tw"][okw])) if okw.any() else float("nan")
    if "VTAuncross_R" in out:
        V = out["VTAuncross_R"]
        kv = retained_nodes(V["dist"], V["isovf"], None)
        row["twICVF_VTAuncross_R"] = float(np.nanmean(V["tw"][kv])) if kv.any() else float("nan")
    # split-half streamline reliability (20 halvings), per session values kept in memory for the pooled aggregate
    rng = np.random.default_rng(L.SEED)
    halves = []
    for _ in range(20):
        perm = rng.permutation(len(S["res"]))
        h = len(perm) // 2
        vals = []
        for part in (perm[:h], perm[h:2 * h]):
            sub_res = [S["res"][j] for j in part]
            if len(sub_res) < 5:
                vals.append(float("nan"))
                continue
            try:
                pm = profile(icvf * tf, icvf_img.affine, sub_res)
                pt = profile(tf, isovf_img.affine, sub_res)
                vals.append(float(np.nanmean((pm / pt)[keep])) if keep.any() else float("nan"))
            except Exception:  # noqa: BLE001 (a degenerate half-bundle is recorded as missing)
                vals.append(float("nan"))
        halves.append(vals)
    return row, {"nodes": nodes, "tw_prof": S["tw"], "keep": keep, "halves": halves}


# ------------------------------------------------------------------------------------------ stages

def read_candidates(p):
    with open(p, newline="") as f:
        return {(r["participant_id"], r["session_id"]) for r in csv.DictReader(f)}


def stage_geometry(args, cfg):
    import numpy as np
    lx = L.lx_root()
    outd = os.path.join(lx, "dwi")
    L.forbid_label_inputs([cfg["pennlinc_qc"], args.s3_candidates])
    found = sorted(glob.glob(cfg["trk"].format(**dict(cfg, sub="sub-*", ses="ses-*", label=BUNDLES["SCP_R"])), recursive=True))
    import re
    sess = sorted({m.groups() for m in (re.search(r"(sub-[A-Za-z0-9]+)_(ses-[A-Za-z0-9]+)", os.path.basename(p)) for p in found) if m})
    L.say(f"stage geometry: {len(sess)} DWI sessions with an SCP_R bundle -> {outd}", args.dry_run)
    if args.dry_run:
        return
    for fn in ("dwi_sessions.csv", "dwi_nodes.csv", "n_min.json"):
        if os.path.exists(os.path.join(outd, fn)):
            sys.exit(f"{os.path.join(outd, fn)} exists (computed once)")
    tpl = template_points(cfg)
    rr = L.RunRecord(__file__, False)
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        res = list(ex.map(session_job, [(s, k, cfg, tpl, args.test_identity_points, not args.no_n4) for s, k in sess]))
    rows = [r for r, _ in res]
    extra = {(r["participant_id"], r["session_id"]): e for r, e in res if e is not None}
    # peak alignment against the pooled mean SCP_R tw-ICVF profile (label-free)
    if extra:
        tmpl = np.nanmean(np.stack([e["tw_prof"] for e in extra.values()]), axis=0)
        for r in rows:
            e = extra.get((r["participant_id"], r["session_id"]))
            if e is None:
                continue
            s = peak_shift(e["tw_prof"], tmpl)
            idx = np.arange(N_NODES)[e["keep"]] + s
            idx = idx[(idx >= 0) & (idx < N_NODES)]
            r["peak_shift"] = s
            r["twICVF_pSCP_R_peak"] = float(np.nanmean(e["tw_prof"][idx])) if idx.size else float("nan")
    # PennLINC QC (scalar pass, FD, NDC)
    qc = {}
    for q in L.read_csv(cfg["pennlinc_qc"]):
        k = (f"sub-{q['participant_id'].strip()}", f"ses-{q['session_id'].strip()}")
        qc[k] = q
    for r in rows:
        q = qc.get((r["participant_id"], r["session_id"]), {})
        r["scalar_pass"] = int((q.get("qc_determination_scalar_maps") or "").strip() == "pass")
        r["FD"] = q.get("mean_FD", "NA")
        r["NDC"] = q.get("raw_neighbor_corr", "NA")
    # n_min over the pre-QC S3 candidates
    cand = read_candidates(args.s3_candidates)
    cr = [r for r in rows if (r["participant_id"], r["session_id"]) in cand and r.get("status") == "ok"]
    counts = {n: sum(1 for r in cr if r.get("n_retained", 0) >= n) for n in NMIN_GRID}
    nmin = next((n for n in NMIN_GRID if counts[n] >= KEEP_FRAC * len(cand)), None)
    # split-half reliability pooled over sessions
    sh = []
    H = [e["halves"] for e in extra.values()]
    for k in range(20):
        a = np.array([h[k][0] for h in H], float)
        b = np.array([h[k][1] for h in H], float)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() > 10:
            rr_ = np.corrcoef(a[ok], b[ok])[0, 1]
            sh.append(2 * rr_ / (1 + rr_))
    summ = {"utc": L.utcnow(), "n_sessions": len(rows), "status_counts": {s: sum(r["status"] == s for r in rows) for s in sorted({r["status"] for r in rows})},
            "n_candidates": len(cand), "n_candidates_ok": len(cr), "kept_at_n_min": counts, "n_min": nmin,
            "rule": "largest of {10, 8, 6} keeping >= 90 % of the pre-QC S3 candidate sessions after rules (i)-(iii)",
            "split_half_SB_mean_of_20": float(np.mean(sh)) if sh else None,
            "fallback_required": nmin is None}
    cols = ["participant_id", "session_id", "status", "n_int", "n_retained", "tract_len", "nSL", "twICVF_pSCP_R", "b0loc_pSCP_R",
            "twICVF_MCP_R_matched", "twICVF_MCP_R_whole", "twICVF_VTAuncross_R", "peak_shift", "twICVF_pSCP_R_peak",
            "scalar_pass", "FD", "NDC"]
    L.write_csv(os.path.join(outd, "dwi_sessions.csv"), rows, cols)
    ncols = ["participant_id", "session_id", "bundle", "node", "cx", "cy", "cz", "dist_dentate_mm", "icvf", "isovf", "tf",
             "twicvf", "b0loc", "retained"]
    L.write_csv(os.path.join(outd, "dwi_nodes.csv"), [n for e in extra.values() for n in e["nodes"]], ncols)
    L.write_json(os.path.join(outd, "n_min.json"), summ)
    for fn in ("dwi_sessions.csv", "dwi_nodes.csv", "n_min.json"):
        rr.add_output(os.path.join(outd, fn))
    rr.finish()
    print(json.dumps(summ, indent=1))


def stage_table(args, cfg):
    lx = L.lx_root()
    outd = os.path.join(lx, "dwi")
    nm = json.load(open(os.path.join(outd, "n_min.json")))
    nmin = nm["n_min"]
    if nmin is None:
        sys.exit("n_min undefined: the registered fallback applies (S3 exploratory); no confirmatory S3 table")
    bq = {(r["participant_id"], r["session_id"]): r for r in L.read_csv(args.bundle_qc)}
    afq = set()
    if args.round1_list and os.path.exists(args.round1_list):
        afq = {ln.split("#")[0].strip() for ln in open(args.round1_list) if ln.split("#")[0].strip().startswith("sub-")}
    rows = []
    for r in L.read_csv(os.path.join(outd, "dwi_sessions.csv")):
        k = (r["participant_id"], r["session_id"])
        q = bq.get(k, {})
        ok = r["status"] == "ok" and r["scalar_pass"] == "1" and q.get("scp_r_pass") == "1" and \
            r["n_retained"] not in ("", "NA") and int(float(r["n_retained"])) >= nmin
        r["s3_eligible"] = int(ok)
        r["afq_round1_list"] = int(r["participant_id"] in afq)
        rows.append(r)
    out = os.path.join(outd, "dwi_s3_table.csv")
    if os.path.exists(out):
        sys.exit(f"{out} exists (frozen once)")
    cols = ["participant_id", "session_id", "s3_eligible", "twICVF_pSCP_R", "b0loc_pSCP_R", "tract_len", "nSL", "FD", "NDC",
            "twICVF_MCP_R_matched", "twICVF_VTAuncross_R", "twICVF_pSCP_R_peak", "t1t2_pSCP_R", "afq_round1_list", "n_retained"]
    t1t2 = os.path.join(outd, "dwi_t1t2.csv")
    if os.path.exists(t1t2):
        tt = {(x["participant_id"], x["session_id"]): x for x in L.read_csv(t1t2)}
        for r in rows:
            r["t1t2_pSCP_R"] = tt.get((r["participant_id"], r["session_id"]), {}).get("t1t2_pSCP_R", "NA")
    L.write_csv(out, rows, cols)
    L.say(f"wrote {out}: {len(rows)} sessions, {sum(r['s3_eligible'] for r in rows)} S3-eligible (n_min {nmin})")


def stage_t1t2(args, cfg):
    """Calibrated T1w/T2w along SCP_R (optional S3 sensitivity)."""
    import nibabel as nib
    import numpy as np
    from dipy.io.streamline import load_tractogram
    from scipy import ndimage
    lx = L.lx_root()
    outd = os.path.join(lx, "dwi")
    nodes = {}
    for n in L.read_csv(os.path.join(outd, "dwi_nodes.csv")):
        if n["bundle"] == "SCP_R":
            nodes.setdefault((n["participant_id"], n["session_id"]), []).append(n)
    rows = []
    for (sub, ses), nl in sorted(nodes.items()):
        fm = dict(cfg, sub=sub, ses=ses)
        t1, t2, asegp = one(cfg["t1w"].format(**fm)), one(cfg["t2w"].format(**fm)), one(cfg["aseg"].format(**fm))
        trk = one(cfg["trk"].format(**fm, label=BUNDLES["SCP_R"]))
        row = {"participant_id": sub, "session_id": ses, "t1t2_pSCP_R": float("nan")}
        if not (t1 and t2 and asegp and trk):
            rows.append(row)
            continue
        with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as td:
            pre = os.path.join(td, "t2_to_t1_")
            if not args.test_identity_points:
                subprocess.run(["antsRegistrationSyN.sh", "-d", "3", "-f", t1, "-m", t2, "-o", pre, "-t", "r", "-n", "2"],
                               check=True, stdout=subprocess.DEVNULL)
                t2r = pre + "Warped.nii.gz"
            else:
                t2r = t2
            T1 = np.asarray(nib.load(t1).get_fdata(), float)
            T2 = np.asarray(nib.load(t2r).get_fdata(), float)
            aimg = nib.load(asegp)
            aseg = np.asarray(aimg.dataobj).astype(int)
            ratio = np.where(T2 > 0, T1 / np.maximum(T2, 1e-6), np.nan)
            csf = ndimage.binary_erosion(np.isin(aseg, [4, 43]), iterations=1)
            wm = ndimage.binary_erosion(np.isin(aseg, [2, 41]), iterations=2)
            if csf.sum() < 20 or wm.sum() < 20 or ratio.shape != aseg.shape:
                rows.append(row)
                continue
            c0, w1 = np.nanmedian(ratio[csf]), np.nanmedian(ratio[wm])
            cal = (ratio - c0) / (w1 - c0)
            sft = load_tractogram(trk, "same", bbox_valid_check=False)
            sft.to_rasmm()
            cen = np.array([[float(n["cx"]), float(n["cy"]), float(n["cz"])] for n in sorted(nl, key=lambda x: int(x["node"]))])
            anchor = cen[0]
            res = orient_and_resample(list(sft.streamlines), anchor)
            p = profile(np.nan_to_num(cal), nib.load(t1).affine, res)
            keep = np.array([n.get("retained") == "1" for n in sorted(nl, key=lambda x: int(x["node"]))])
            row["t1t2_pSCP_R"] = float(np.nanmean(p[keep])) if keep.any() else float("nan")
        rows.append(row)
    L.write_csv(os.path.join(outd, "dwi_t1t2.csv"), rows, ["participant_id", "session_id", "t1t2_pSCP_R"])
    L.say(f"wrote {os.path.join(outd, 'dwi_t1t2.csv')} ({len(rows)} sessions)")


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--stage", choices=["geometry", "t1t2", "table"], required=True)
    a.add_argument("--config", default=None, help="JSON overriding the input patterns (VERIFY at install)")
    a.add_argument("--s3-candidates", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "s3_candidates_preqc.csv"))
    a.add_argument("--bundle-qc", default=None, help="frozen DWI bundle-QC table (participant_id, session_id, scp_r_pass)")
    a.add_argument("--round1-list", default=None, help="the 20-subject round-1 AFQ-QC list (conservative S3 sample)")
    a.add_argument("--jobs", type=int, default=int(os.environ.get("LX_MAXJOBS", "8")))
    a.add_argument("--no-n4", action="store_true", help="skip N4 on the b0 (test only)")
    a.add_argument("--test-identity-points", action="store_true", help="synthetic test only (requires LX_TEST=1)")
    a.add_argument("--dry-run", action="store_true")
    args = a.parse_args()
    if args.test_identity_points and os.environ.get("LX_TEST") != "1":
        sys.exit("--test-identity-points requires LX_TEST=1")
    cfg = cfg_default(L.lx_root())
    if args.config:
        cfg.update(json.load(open(args.config)))
    {"geometry": stage_geometry, "t1t2": stage_t1t2, "table": stage_table}[args.stage](args, cfg)


if __name__ == "__main__":
    main()
