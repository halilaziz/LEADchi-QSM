#!/usr/bin/env python3
"""blind_setup.py - coded, header-scrubbed copies of the Penn LEAD MEGRE and anatomy (task P5b, 2026-09-25).

Procedural blinding, v1.2 delta sec. 2.3-2.4 (user decision U1: no escrow holder). Written and hashed before the
OSF registration; run ONCE on server (python from AFQ_env: numpy, nibabel) after approvals A4/A6, before posting.

What it does
  1. Enumerates every ds007116 session with all 8 MEGRE NIfTI + JSON (expected 181; a partial set aborts).
  2. Validates the sidecars: TE = 7.5/15/22.5/30 ms for every echo; one ImagingFrequency per session.
  3. Draws random subject codes (secrets.token_hex(3); >= 1 letter, no 4 consecutive digits). Coded ses-<k> keep
     within-subject chronological order (k = rank of the BIDS session number). Processing order is randomised.
  4. Writes, per session, into $LX/blind/bids_all (MEGRE only; QSMxT input) and $LX/blind/anat_all:
       MEGRE magnitude  float32, header scrubbed, geometry unchanged
       MEGRE phase      float32 radians = raw * pi / 4096 (registered --phase-mode signed12, verified per echo),
                        so neither QSMxT nor SEPIA rescales it (E3 finding A3)
       JSON sidecars    reduced to EchoTime, RepetitionTime, FlipAngle and
                        MagneticFieldStrength = round(ImagingFrequency / 42.576, 4) T  (QSMxT's gamma, so that
                        gamma_engine * B0 = f0; E3 finding A1; rounding error negligible)
       raw T1w / T2w    header scrubbed (anatomy for the ROI chain and T1w/T2w)
       QSIPrep anat     ACPC T1w and SynthSeg aseg / dseg (images only; no QSIPrep transform is copied: the
                        registered route is R, own SyN, for every session)
     Every image of a session receives the SAME random world translation (secrets; integer mm in [-40, 40] per
     axis) in qform and sform: relative geometry, rotations and the B0 direction are unchanged, but a coded header no
     longer matches its raw header (registration Annex A.3 item 5). The offsets are stored only in the encrypted key.
     NIfTI descrip (holds acquisition time of day), aux_file, intent_name, db_name, data_type and all header
     extensions are blanked; every written file is scanned for the real subject ID.
  5. Builds the key [code, ses_order, participant_id, session_id, epoch, f0_cluster, offset_x/y/z] in memory only
     (epoch = 'early' for ProtocolName QSM_SWI_1.5mm; f0 cluster = rank of the ImagingFrequency cluster, clusters
     split at gaps > 0.01 MHz) and writes it encrypted to the
     user's lock_pub.pem (lx_common.hybrid_encrypt): blind/key.csv.enc + blind/key.sessionkey.enc.
  6. Writes the visible, label-free per-code flags (blind/session_flags.csv: T1w presence and PennLINC T1 QC
     status, QSIPrep anatomy availability, SoftwareVersions equality class), the processing order, the pilot
     (10), L-curve (2 of the pilot) and G2 (20, distinct subjects, disjoint from the pilot) draws (seed 20260924),
     coded_files.sha256 and setup_manifest.json (counts and hashes only).

It never opens participants.tsv, *_sessions.tsv or phenotype files, prints no code-ID pair and no per-session
value. Modes: --dry-run (metadata only, prints the plan, writes nothing); --header-audit (NIfTI headers + JSON
only, no voxel data; aggregate summary, used in week 1 to settle scl_slope/scl_inter and SoftwareVersions, the
slab obliquity, the sequence parameters and the uniqueness of the per-code visible channels).
"""
from __future__ import annotations

import argparse
import collections
import glob
import io
import json
import os
import random
import re
import secrets
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

PARTS = [(n, p) for n in (1, 2, 3, 4) for p in ("mag", "phase")]
SEQ_KEYS = ("CoilCombinationMethod", "ProtocolName", "ReceiveCoilName", "SequenceName", "SequenceVariant", "ScanOptions",
            "PixelBandwidth", "PhaseEncodingDirection", "InPlanePhaseEncodingDirectionDICOM", "ParallelReductionFactorInPlane",
            "ParallelAcquisitionTechnique", "PartialFourier", "PercentPhaseFOV", "PercentSampling", "MRAcquisitionType",
            "AcquisitionDuration", "BandwidthPerPixelPhaseEncode", "MTState", "SpoilingState")
EARLY_PROTOCOL = "QSM_SWI_1.5mm"


def f0_clusters(values, gap=0.01):
    """Rank of the ImagingFrequency cluster (distinct values sorted, split at gaps > `gap` MHz), 1-based."""
    v = sorted(set(values))
    lab, k = {}, 1
    for i, x in enumerate(v):
        if i and x - v[i - 1] > gap:
            k += 1
        lab[x] = k
    return lab


def obliquity_deg(aff):
    """Largest angle between a voxel axis and its nearest world axis, and the angle of the slice normal to world z."""
    import numpy as np
    R = aff[:3, :3] / np.linalg.norm(aff[:3, :3], axis=0)
    worst = max(float(np.degrees(np.arccos(np.clip(np.max(np.abs(R[:, j])), -1, 1)))) for j in range(3))
    return worst, float(np.degrees(np.arccos(np.clip(abs(R[2, 2]), -1, 1))))


def parse_args():
    lx, pl = L.lx_root(), L.pl_root()
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--bids", default=os.path.join(pl, "ds007116"))
    a.add_argument("--qsiprep", default=os.path.join(pl, "QSIprep"))
    a.add_argument("--smriprep", default=None,
                   help="optional dir with sMRIPrep/fMRIPrep anat transforms (route S; only if downloaded)")
    a.add_argument("--t1qc", default=os.path.join(lx, "prereg", "inputs", "EF_desc-T1_qc.csv"))
    a.add_argument("--out", default=os.path.join(lx, "blind"))
    a.add_argument("--lock-pubkey", default=os.path.join(lx, "prereg", "lock_pub.pem"))
    a.add_argument("--seed", type=int, default=L.SEED)
    a.add_argument("--pilot", type=int, default=10)
    a.add_argument("--lcurve", type=int, default=2)
    a.add_argument("--gate", type=int, default=20)
    a.add_argument("--expect-sessions", type=int, default=181)
    a.add_argument("--phase-mode", choices=["signed12", "unsigned12", "radians", "auto12"], default="signed12")
    a.add_argument("--xfm-route", choices=["R"], default="R", help="registered: R for every session")
    a.add_argument("--offset-mm", type=int, default=40, help="max |translation| per axis of the per-session offset")
    a.add_argument("--dry-run", action="store_true")
    a.add_argument("--header-audit", action="store_true")
    return a.parse_args()


# ----------------------------------------------------------------------------------------------- inventory

def one(pattern: str):
    m = sorted(glob.glob(pattern))
    return m


def enumerate_sessions(bids: str) -> tuple[list[dict], list[str]]:
    sessions, partial = [], []
    for sd in sorted(glob.glob(os.path.join(bids, "sub-*"))):
        sub = os.path.basename(sd)
        if not os.path.isdir(sd):
            continue
        for ed in sorted(glob.glob(os.path.join(sd, "ses-*"))):
            ses = os.path.basename(ed)
            anat = os.path.join(ed, "anat")
            files, jsons = {}, {}
            for n, p in PARTS:
                m = one(os.path.join(anat, f"{sub}_{ses}*echo-{n}_part-{p}_MEGRE.nii.gz"))
                if len(m) == 1:
                    files[(n, p)] = m[0]
                    j = m[0][:-7] + ".json"
                    if os.path.exists(j):
                        jsons[(n, p)] = j
                elif len(m) > 1:
                    sys.exit(f"{sub} {ses}: {len(m)} files match echo-{n} part-{p}; ambiguous")
            if not files:
                continue
            if len(files) != 8 or len(jsons) != 8:
                partial.append(f"{sub}_{ses}")
                continue
            mnum = re.fullmatch(r"ses-(\d+)", ses)
            if not mnum:
                sys.exit(f"{sub} {ses}: non-numeric session label; chronological order undefined")
            sessions.append({"sub": sub, "ses": ses, "ses_num": int(mnum.group(1)), "anat": anat,
                             "files": files, "jsons": jsons})
    return sessions, partial


def read_meta(s: dict) -> dict:
    te, tr, fa, f0, fs, sw = {}, set(), set(), set(), set(), set()
    proto = None
    extra = collections.Counter()
    for (n, p), j in s["jsons"].items():
        with open(j, encoding="utf-8") as f:
            d = json.load(f)
        te[(n, p)] = float(d["EchoTime"])
        tr.add(round(float(d["RepetitionTime"]), 6))
        fa.add(round(float(d["FlipAngle"]), 3))
        if "ImagingFrequency" not in d:
            sys.exit(f"{s['sub']} {s['ses']}: ImagingFrequency missing in {os.path.basename(j)}")
        f0.add(round(float(d["ImagingFrequency"]), 6))
        fs.add(d.get("MagneticFieldStrength"))
        sw.add(str(d.get("SoftwareVersions")))
        for k in SEQ_KEYS:
            if k in d and (n, p) == (1, "phase"):
                extra[(k, str(d[k]))] += 1
        if (n, p) == (1, "phase"):
            proto = str(d.get("ProtocolName"))
    for (n, p), v in te.items():
        if abs(v - L.TE_REGISTERED_S[n - 1]) > 1e-6:
            sys.exit(f"{s['sub']} {s['ses']}: echo-{n} {p} EchoTime {v} s differs from the registered TE set")
    if len(tr) != 1 or len(fa) != 1 or len(f0) != 1 or len(sw) != 1:
        sys.exit(f"{s['sub']} {s['ses']}: TR/FA/ImagingFrequency/SoftwareVersions differ between echoes")
    f0v = f0.pop()
    return {"protocol": proto, "TE": [round(L.TE_REGISTERED_S[n - 1], 6) for n in (1, 2, 3, 4)], "TR": tr.pop(), "FA": fa.pop(),
            "f0_MHz": f0v, "B0_T": round(f0v / L.GAMMA_QSMXT_MHZ_T, 4), "field_json": fs.pop(),
            "sw": sw.pop(), "extra": extra}


def anat_inventory(s: dict, qsiprep: str, smriprep: str | None) -> dict:
    sub, ses = s["sub"], s["ses"]

    def pick(pattern):
        m = sorted(glob.glob(pattern))
        if not m:
            return None
        r1 = [x for x in m if "run-01" in os.path.basename(x)]
        return (r1 or m)[0], len(m)

    inv = {}
    t1 = pick(os.path.join(s["anat"], f"{sub}_{ses}*_T1w.nii.gz"))
    t2 = pick(os.path.join(s["anat"], f"{sub}_{ses}*_T2w.nii.gz"))
    inv["t1w_raw"], inv["t1w_n"] = (t1 if t1 else (None, 0))
    inv["t2w_raw"], inv["t2w_n"] = (t2 if t2 else (None, 0))
    qa = os.path.join(qsiprep, sub, ses, "anat")
    q = {
        "t1w_acpc": f"{sub}_{ses}_space-ACPC_desc-preproc_T1w.nii.gz",
        "aseg_acpc": f"{sub}_{ses}_space-ACPC_desc-aseg_dseg.nii.gz",
        "dseg_acpc": f"{sub}_{ses}_space-ACPC_dseg.nii.gz",
        "xfm_acpc2anat": f"{sub}_{ses}_from-ACPC_to-anat_mode-image_xfm.mat",
        "xfm_mni2acpc": f"{sub}_{ses}_from-MNI152NLin2009cAsym_to-ACPC_mode-image_xfm.h5",
    }
    for k, fn in q.items():
        p = os.path.join(qa, fn)
        inv[k] = p if os.path.exists(p) else None
    inv["qsiprep_anat"] = all(inv[k] for k in ("t1w_acpc", "aseg_acpc"))   # images only (route R needs no transform)
    inv["xfm_mni2t1w_smriprep"] = None
    if smriprep:
        m = sorted(glob.glob(os.path.join(smriprep, sub, ses, "anat",
                                          f"{sub}_{ses}*_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5")))
        inv["xfm_mni2t1w_smriprep"] = m[0] if m else None
    return inv


def read_t1qc(path: str) -> dict:
    out = {}
    if not os.path.exists(path):
        sys.exit(f"T1 QC table {path} not found (copy PennLINC EF_desc-T1_qc.csv there; public QC file)")
    for r in L.read_csv(path):
        out[(r["participant_id"].strip(), r["session_id"].strip())] = r["qc_determination"].strip().lower()
    return out


def draw_code(taken: set) -> str:
    while True:
        c = secrets.token_hex(3)
        if re.search(r"[a-f]", c) and not re.search(r"\d{4}", c) and c not in taken:
            taken.add(c)
            return c


# ----------------------------------------------------------------------------------------------- image IO

def phase_to_rad(arr, mode: str):
    import numpy as np
    mn, mx = float(np.nanmin(arr)), float(np.nanmax(arr))
    if mode == "auto12":   # amendment-1 technical deviation: per-file Siemens 12-bit encoding (175 signed, 6 unsigned sessions)
        mode = "signed12" if mn < 0 else "unsigned12"
    if mode == "signed12":
        ok = mn >= -4096.5 and mx <= 4095.5 and mn < -2048 and mx > 2048
        out = arr * (np.pi / 4096.0)
    elif mode == "unsigned12":
        ok = mn >= -0.5 and mx <= 4095.5 and mx > 2048
        out = (arr - 2048.0) * (np.pi / 2048.0)
    else:
        ok = mn >= -np.pi - 1e-3 and mx <= np.pi + 1e-3
        out = arr
    return np.asarray(out, dtype=np.float32), ok, mn, mx


def intensity_data(img):
    import numpy as np
    dt = img.get_data_dtype()
    slope, inter = img.header.get_slope_inter()
    if np.issubdtype(dt, np.integer) and slope in (None, 1.0) and inter in (None, 0.0):
        return np.asarray(img.dataobj).astype(dt)
    return np.asarray(img.get_fdata(dtype=np.float32), dtype=np.float32)


def scan_ids(path: str, tokens: list[bytes]) -> list[str]:
    if path.endswith((".nii", ".nii.gz")):
        blob = L.header_text_blob(path)
    else:
        with open(path, "rb") as f:
            blob = f.read()
    return [t.decode() for t in tokens if t in blob]


# ----------------------------------------------------------------------------------------------- header audit

def header_audit(sessions: list[dict], args) -> None:
    import nibabel as nib
    import numpy as np
    agg = collections.Counter()
    b0 = collections.Counter()
    f0 = collections.Counter()
    sw = collections.Counter()
    extra = collections.Counter()
    obl, obl_z = [], []
    meta = {}
    for s in sessions:
        m = read_meta(s)
        meta[(s["sub"], s["ses"])] = m
        f0[m["f0_MHz"]] += 1
        b0[m["B0_T"]] += 1
        sw[m["sw"]] += 1
        extra.update(m["extra"])
        for (n, p), fpath in s["files"].items():
            h = nib.load(fpath).header          # header only; no voxel data are read
            sl, it = h.get_slope_inter()
            agg[("part", p, "dtype", str(h.get_data_dtype()))] += 1
            agg[("part", p, "scl_slope", str(sl), "scl_inter", str(it))] += 1
            agg[("part", p, "shape", str(h.get_data_shape()), "zooms", str(tuple(round(float(z), 3) for z in h.get_zooms())))] += 1
            agg[("part", p, "qform_code", int(h["qform_code"]), "sform_code", int(h["sform_code"]))] += 1
            agg[("part", p, "descrip_has_Time", b"Time=" in bytes(h["descrip"]))] += 1
            agg[("part", p, "n_extensions", len(h.extensions))] += 1
            if (n, p) == (1, "phase"):
                w, z = obliquity_deg(h.get_best_affine())
                obl.append(w)
                obl_z.append(z)
    # uniqueness of the visible per-code channels (sessions per subject, f0 cluster, software class, epoch)
    fcl = f0_clusters([m["f0_MHz"] for m in meta.values()])
    nses = collections.Counter(s["sub"] for s in sessions)
    tup = collections.Counter()
    keyed = {}
    for s in sessions:
        m = meta[(s["sub"], s["ses"])]
        t = (nses[s["sub"]], fcl[m["f0_MHz"]], m["sw"], m["protocol"] == EARLY_PROTOCOL)
        keyed[(s["sub"], s["ses"])] = t
        tup[t] += 1
    uniq = sum(1 for k in keyed if tup[keyed[k]] == 1)
    q = lambda a: [float(np.min(a)), float(np.median(a)), float(np.max(a))] if a else None
    out = {"n_sessions": len(sessions), "f0_MHz_counts": dict(sorted(f0.items())),
           "f0_cluster_counts": dict(collections.Counter(fcl[m["f0_MHz"]] for m in meta.values())),
           "B0_T_qsmxt_counts": dict(sorted(b0.items())), "SoftwareVersions_counts": dict(sw),
           "obliquity_worst_axis_deg_min_median_max": q(obl), "slice_normal_vs_world_z_deg_min_median_max": q(obl_z),
           "n_sessions_obliquity_gt_0p5deg": int(sum(x > 0.5 for x in obl)),
           "visible_channel_uniqueness": {"channels": "sessions per subject x f0 cluster x software class x epoch",
                                          "n_sessions_unique": uniq, "n_distinct_combinations": len(tup)},
           "sidecar_field_counts": {f"{k}={v}": c for (k, v), c in sorted(extra.items())},
           "nifti_header_counts": {" | ".join(map(str, k)): c for k, c in sorted(agg.items(), key=str)}}
    print(json.dumps(out, indent=2, default=str))


# ----------------------------------------------------------------------------------------------- main

def main():
    args = parse_args()
    dry = args.dry_run
    L.forbid_label_inputs([args.t1qc])
    sessions, partial = enumerate_sessions(args.bids)
    L.say(f"MEGRE sessions complete: {len(sessions)}; partial: {len(partial)}", dry)
    if partial:
        sys.exit("partial MEGRE sessions found (the inventory expected none): " + ", ".join(partial))
    if args.header_audit:
        header_audit(sessions, args)
        return
    if len(sessions) != args.expect_sessions:
        sys.exit(f"expected {args.expect_sessions} MEGRE sessions, found {len(sessions)}")
    out = os.path.realpath(args.out)
    L.assert_inside(out)
    if os.path.exists(os.path.join(out, "key.csv.enc")):
        sys.exit(f"REFUSED: {out}/key.csv.enc exists; blind setup runs once")
    rr = L.RunRecord(__file__, dry)
    if not dry:
        if not os.path.exists(args.lock_pubkey):
            sys.exit(f"{args.lock_pubkey} missing: the user generates the lock key pair first (README step 1)")
        bits = L.rsa_bits(args.lock_pubkey)
        if bits < 3072:
            sys.exit(f"lock public key has {bits} bits; >= 3072 required")
    t1qc = read_t1qc(args.t1qc)
    rr.add_input(args.t1qc)

    # ---- metadata and anatomy inventory (no voxel data)
    meta, inv = {}, {}
    for s in sessions:
        k = (s["sub"], s["ses"])
        meta[k] = read_meta(s)
        inv[k] = anat_inventory(s, args.qsiprep, args.smriprep)
    with_t1 = [k for k in inv if inv[k]["qsiprep_anat"] or inv[k]["t1w_raw"]]
    n_q = sum(1 for k in with_t1 if inv[k]["qsiprep_anat"] and inv[k]["xfm_mni2acpc"])
    n_s = 0
    route = "R"
    L.say(f"registration route: R (registered; own ANTs SyN for all {len(with_t1)} sessions with a T1w)", dry)
    fcl = f0_clusters([m["f0_MHz"] for m in meta.values()])
    offsets = {k: [secrets.randbelow(2 * args.offset_mm + 1) - args.offset_mm for _ in range(3)] for k in meta}

    # ---- codes, chronological coded session index, randomised processing order
    taken: set = set()
    subs = sorted({s["sub"] for s in sessions})
    code = {sub: draw_code(taken) for sub in subs}
    by_sub = collections.defaultdict(list)
    for s in sessions:
        by_sub[s["sub"]].append(s)
    ses_order = {}
    for sub, lst in by_sub.items():
        for k, s in enumerate(sorted(lst, key=lambda x: x["ses_num"]), start=1):
            ses_order[(sub, s["ses"])] = k
    order = list(sessions)
    secrets.SystemRandom().shuffle(order)

    bids_all = os.path.join(out, "bids_all")
    anat_all = os.path.join(out, "anat_all")
    written = []
    phase_ok_all, phase_min, phase_max = True, float("inf"), float("-inf")
    flags = []
    sw_classes = {v: chr(ord("A") + i) for i, v in enumerate(sorted({m["sw"] for m in meta.values()}))}
    t1qc_counts = collections.Counter()

    import nibabel as nib  # noqa: E402  (needed only past this point)

    for i, s in enumerate(order, start=1):
        sub, ses = s["sub"], s["ses"]
        c, k = code[sub], ses_order[(sub, ses)]
        lab = L.sess_label(c, k)
        m, a = meta[(sub, ses)], inv[(sub, ses)]
        tokens = [sub.encode(), sub.replace("sub-", "").encode()]
        dmeg = os.path.join(bids_all, f"sub-{c}", f"ses-{k}", "anat")
        dana = os.path.join(anat_all, f"sub-{c}", f"ses-{k}", "anat")
        ref_aff, ref_shape = None, None
        sess_modes = set()
        for n, p in PARTS:
            src = s["files"][(n, p)]
            dst = os.path.join(dmeg, f"{lab}_echo-{n}_part-{p}_MEGRE.nii.gz")
            if dry:
                if i <= 1:
                    L.say(f"would write {os.path.join(bids_all, 'sub-<code>', 'ses-<k>', 'anat', os.path.basename(dst).replace(lab, 'sub-<code>_ses-<k>'))}", dry)
                continue
            img = nib.load(src)
            if ref_aff is None:
                ref_aff, ref_shape = img.affine, img.shape
            elif img.shape != ref_shape or not (abs(img.affine - ref_aff) < 1e-5).all():
                sys.exit(f"{sub} {ses}: MEGRE files do not share one grid")
            if p == "phase":
                data, ok, mn, mx = phase_to_rad(img.get_fdata(dtype="float64"), args.phase_mode)
                if args.phase_mode == "auto12":
                    sess_modes.add("signed12" if mn < 0 else "unsigned12")
                    if len(sess_modes) > 1:
                        sys.exit(f"{sub} {ses}: echoes use different 12-bit phase encodings; stop")
                phase_ok_all &= ok
                phase_min, phase_max = min(phase_min, mn), max(phase_max, mx)
                if not ok:
                    sys.exit(f"{sub} {ses} echo-{n}: phase range does not match --phase-mode "
                             f"{args.phase_mode}; stop, delete {out} yourself and re-run the header audit")
            else:
                data = img.get_fdata(dtype="float32").astype("float32")
            L.write_scrubbed(img, data, dst, offset_mm=offsets[(sub, ses)])
            js = dst[:-7] + ".json"
            with open(js, "w", encoding="utf-8") as f:
                json.dump({"EchoTime": m["TE"][n - 1], "RepetitionTime": m["TR"], "FlipAngle": m["FA"],
                           "MagneticFieldStrength": m["B0_T"]}, f, indent=2)
            written += [dst, js]
        # anatomy
        copies = []
        if a["t1w_raw"]:
            copies.append(("img", a["t1w_raw"], os.path.join(dana, f"{lab}_T1w.nii.gz")))
        if a["t2w_raw"]:
            copies.append(("img", a["t2w_raw"], os.path.join(dana, f"{lab}_T2w.nii.gz")))
        if a["qsiprep_anat"]:
            copies += [("img", a["t1w_acpc"], os.path.join(dana, f"{lab}_space-ACPC_desc-preproc_T1w.nii.gz")),
                       ("lab", a["aseg_acpc"], os.path.join(dana, f"{lab}_space-ACPC_desc-aseg_dseg.nii.gz"))]
            if a["dseg_acpc"]:
                copies.append(("lab", a["dseg_acpc"], os.path.join(dana, f"{lab}_space-ACPC_dseg.nii.gz")))
        for kind, src, dst in copies:
            if dry:
                continue
            L.ensure_parent(dst)
            if kind == "raw":
                shutil.copyfile(src, dst)
            else:
                img = nib.load(src)
                data = intensity_data(img) if kind == "img" else L.load_labels(src)[1].astype("int16")
                L.write_scrubbed(img, data, dst, offset_mm=offsets[(sub, ses)])
            written.append(dst)
        if not dry:
            for f in [x for x in written if os.path.basename(x).startswith(lab + "_")]:
                hit = scan_ids(f, tokens)
                if hit:
                    sys.exit(f"ID string found in a coded file ({os.path.basename(f)}); aborting")
        st = t1qc.get((sub, ses), "missing") if a["t1w_raw"] or a["qsiprep_anat"] else "no_t1w"
        t1qc_counts[st] += 1
        flags.append({"code": c, "ses_order": k,
                      "t1w_present": int(bool(a["t1w_raw"] or a["qsiprep_anat"])),
                      "t1qc": st, "t1_ok": int(st in ("pass", "artifact", "missing")),
                      "t1_artifact": int(st == "artifact"), "t1qc_missing": int(st == "missing"),
                      "t2w_present": int(bool(a["t2w_raw"])), "qsiprep_anat": int(bool(a["qsiprep_anat"])),
                      "anat_src": "qsiprep" if a["qsiprep_anat"] else ("raw" if a["t1w_raw"] else "none"),
                      "sw_class": sw_classes[m["sw"]]})
        if i % 20 == 0 and not dry:
            L.say(f"{i}/{len(order)} sessions written")

    # ---- key (encrypted), flags, draws, manifest
    key_rows = sorted(((code[s["sub"]], ses_order[(s["sub"], s["ses"])], s["sub"], s["ses"],
                        "early" if meta[(s["sub"], s["ses"])]["protocol"] == EARLY_PROTOCOL else "late",
                        fcl[meta[(s["sub"], s["ses"])]["f0_MHz"]], *offsets[(s["sub"], s["ses"])]) for s in sessions))
    buf = io.StringIO()
    buf.write("code,ses_order,participant_id,session_id,epoch,f0_cluster,offset_x,offset_y,offset_z\n")
    for r in key_rows:
        buf.write(",".join(map(str, r)) + "\n")
    enc = L.hybrid_encrypt(buf.getvalue().encode(), args.lock_pubkey, os.path.join(out, "key.csv.enc"),
                           os.path.join(out, "key.sessionkey.enc"), dry=dry)
    del buf, key_rows

    flags.sort(key=lambda r: (r["code"], r["ses_order"]))
    fcols = ["code", "ses_order", "t1w_present", "t1qc", "t1_ok", "t1_artifact", "t1qc_missing", "t2w_present",
             "qsiprep_anat", "anat_src", "sw_class"]
    L.write_csv(os.path.join(out, "session_flags.csv"), flags, fcols, dry)
    po = [f"{code[s['sub']]}\t{ses_order[(s['sub'], s['ses'])]}" for s in order]
    # draws only among sessions that can be registered (T1w present, PennLINC T1 QC not "fail"; label-free)
    eligible = {f"{r['code']}\t{r['ses_order']}" for r in flags if r["t1_ok"] == 1}
    all_sorted = sorted(x for x in po if x in eligible)
    rng = random.Random(args.seed)
    if len(all_sorted) < args.pilot + args.gate:
        sys.exit(f"only {len(all_sorted)} registrable sessions for a pilot of {args.pilot} and a gate set of {args.gate}")
    pilot = rng.sample(all_sorted, args.pilot)
    rest = [x for x in all_sorted if x not in pilot]
    rng.shuffle(rest)
    gate, seen = [], set()
    for x in rest:
        cc = x.split("\t")[0]
        if cc not in seen:
            gate.append(x)
            seen.add(cc)
        if len(gate) == args.gate:
            break
    lists = {"processing_order.txt": po, "pilot10.txt": pilot, "lcurve2.txt": pilot[:args.lcurve],
             "gate20.txt": gate}
    for fn, lst in lists.items():
        p = os.path.join(out, fn)
        L.ensure_parent(p, dry=dry)
        if dry:
            L.say(f"would write {p} ({len(lst)} coded sessions)", dry)
            continue
        with open(p, "w", encoding="utf-8") as f:
            f.write("\n".join(lst) + "\n")
    if not dry:
        with open(os.path.join(bids_all, "dataset_description.json"), "w", encoding="utf-8") as f:
            json.dump({"Name": "LEAD-chi coded MEGRE (procedurally blinded copy of ds007116 v1.0.6)",
                       "BIDSVersion": "1.9.0", "DatasetType": "raw",
                       "GeneratedBy": [{"Name": "blind_setup.py", "Version": rr.rec["script_sha256"][:12]}]},
                      f, indent=2)
        with open(os.path.join(out, "coded_files.sha256"), "w", encoding="utf-8") as f:
            for p in sorted(written):
                f.write(f"{L.sha256_file(p)}  {os.path.relpath(p, out)}\n")
    b0c = collections.Counter(m["B0_T"] for m in meta.values())
    manifest = {
        "utc": L.utcnow(), "script_sha256": rr.rec["script_sha256"], "lx_common_sha256": rr.rec["lx_common_sha256"],
        "n_subjects": len(subs), "n_sessions": len(sessions), "n_files_written": len(written),
        "te_s": list(L.TE_REGISTERED_S), "tr_s_counts": dict(collections.Counter(m["TR"] for m in meta.values())),
        "fa_counts": dict(collections.Counter(m["FA"] for m in meta.values())),
        "b0_qsmxt_T_counts": {str(k): v for k, v in sorted(b0c.items())},
        "gamma_MHz_per_T": L.GAMMA_QSMXT_MHZ_T, "phase_mode": args.phase_mode,
        "phase_raw_range_all_sessions": [phase_min, phase_max] if not dry else None,
        "xfm_route": route, "route_availability": {"with_t1w": len(with_t1), "qsiprep_inverse": n_q,
                                                   "smriprep": n_s},
        "t1qc_counts": dict(t1qc_counts), "software_classes": {v: sw_classes[v] for v in sw_classes},
        "software_class_counts": dict(collections.Counter(sw_classes[m["sw"]] for m in meta.values())),
        "anat_src_counts": dict(collections.Counter(r["anat_src"] for r in flags)),
        "key_csv_enc_sha256": enc.get("data_sha256"), "key_sessionkey_enc_sha256": enc.get("key_sha256"),
        "lock_pub_sha256": L.sha256_file(args.lock_pubkey) if os.path.exists(args.lock_pubkey) else None,
        "seed": args.seed, "draws": {fn: len(lst) for fn, lst in lists.items()},
        "coded_files_sha256_of_list": (L.sha256_file(os.path.join(out, "coded_files.sha256")) if not dry else None),
    }
    L.write_json(os.path.join(out, "setup_manifest.json"), manifest, dry)
    for fn in lists:
        if not dry:
            rr.add_output(os.path.join(out, fn))
    if not dry:
        for fn in ("key.csv.enc", "key.sessionkey.enc", "session_flags.csv", "setup_manifest.json",
                   "coded_files.sha256"):
            rr.add_output(os.path.join(out, fn))
    rr.note(f"route={route}; sessions={len(sessions)}; subjects={len(subs)}")
    rr.finish()
    L.say(json.dumps({k: manifest[k] for k in ("n_subjects", "n_sessions", "n_files_written", "xfm_route",
                                                "b0_qsmxt_T_counts", "t1qc_counts", "key_csv_enc_sha256",
                                                "key_sessionkey_enc_sha256", "lock_pub_sha256")}, indent=2), dry)


if __name__ == "__main__":
    main()
