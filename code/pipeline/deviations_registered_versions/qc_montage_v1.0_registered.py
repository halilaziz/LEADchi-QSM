#!/usr/bin/env python3
"""qc_montage.py - group-blind visual QC material, rating sheets and kappa (task P5b, 2026-09-25).

Implements the procedure in qc_rubric.md (protocol v1.1 sec. 11.4; v1.2 sec. 2.6; gate G4). Runs on server (AFQ_env
python with matplotlib, Agg backend) on the coded tree only. The rater never sees a group label, an age, a code-ID
link, a chi or R2* number, or the within-subject session order: montages are named by random QC IDs
(qc/qc_ids.csv maps them to codes; the rater does not open it) and carry no colour-bar values.

Subcommands
  montage   --kind gross     axial slices at the dentate, SN-VTA and GP levels: echo-1 and echo-4 magnitude, total
                             field, local field, chi (fixed window -0.2..+0.2 ppm, grey) with the Mask-3 outline
            --kind placement T1w (in the MEGRE grid) and echo-1 magnitude with the SN-VTA (w >= 0.5) outline in three
                             planes through its centroid, and the dentate outlines at the dentate level
  sheet     --pass 1         rating sheet for every session, random order (seed 20260924)
            --pass 2         re-rating sheet for a random 20 % (>= 14 days after pass 1 was closed)
            --pass 3         full re-rating after ONE blind rubric refinement (only if G4 failed on passes 1/2)
            --pass 4         new random 20 % re-rating of pass 3 (>= 14 days after pass 3 was closed)
  kappa                      linear-weighted kappa of the 0/1/2 grade (decision, G4 >= 0.60): pass 1 vs 2, or pass
                             3 vs 4 after the refinement; quadratic-weighted and placement kappa co-reported
  finalize                   qc/grades_final.csv (code, ses_order, qc_grade, placement, failure_codes, dental_D) for
                             build_usable.py. Stop rule: at most ONE refinement; if G4 fails again, the pass-3 grades
                             are final, G4 is recorded as failed and every H1 statement carries the caveat.
--dry-run prints the plan only.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import json
import os
import random
import secrets
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

QC = lambda lx, *p: os.path.join(lx, "qc", *p)


def qc_ids(lx, dry):
    p = QC(lx, "qc_ids.csv")
    if os.path.exists(p):
        return {(r["code"], int(r["ses_order"])): r["qc_id"] for r in L.read_csv(p)}
    ids, taken = {}, set()
    for c, s in L.coded_sessions(lx):
        while True:
            q = "q" + secrets.token_hex(4)
            if q not in taken:
                break
        taken.add(q)
        ids[(c, s)] = q
    L.write_csv(p, [{"code": c, "ses_order": s, "qc_id": q} for (c, s), q in sorted(ids.items())],
                ["code", "ses_order", "qc_id"], dry)
    return ids


def _find(lx, variant, code, ses, token):
    lab = L.sess_label(code, ses)
    hits = sorted(h for h in glob.glob(os.path.join(lx, "qsm", variant, "**", f"*{lab}*{token}*.nii*"), recursive=True)
                  if "desc-" not in os.path.basename(h))      # never a desc-singlepass / desc-reliable variant
    return hits[0] if hits else None


def montage_one(job):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import nibabel as nib
    import numpy as np
    lx, kind, variant, code, ses, qid, cfg = job
    lab = L.sess_label(code, ses)
    meg = os.path.join(lx, "blind", "bids_all", f"sub-{code}", f"ses-{ses}", "anat")
    rd = os.path.join(lx, "reg", f"sub-{code}", f"ses-{ses}")
    out = QC(lx, "montages", kind, f"{qid}.png")
    if os.path.exists(out) or not os.path.exists(os.path.join(rd, "route.json")):
        return out, False
    can = lambda p: np.asarray(nib.as_closest_canonical(nib.load(p)).get_fdata(dtype=np.float32))
    e1 = can(os.path.join(meg, f"{lab}_echo-1_part-mag_MEGRE.nii.gz"))
    R = L.build_rois(rd, cfg)
    ref = R["img"]
    tocan = lambda a: np.asarray(nib.as_closest_canonical(nib.Nifti1Image(a.astype(np.float32), ref.affine)).dataobj)
    sn = tocan(R["rois"]["SNVTA"]["w"]) >= 0.5
    dn = tocan(np.maximum(R["rois"]["DN_R"]["w"], R["rois"]["DN_L"]["w"])) >= 0.5
    gp = tocan(R["rois"]["GP"]["w"]) > 0.5

    def zc(m, frac):
        return int(round(np.argwhere(m)[:, 2].mean())) if m.any() else int(frac * (e1.shape[2] - 1))

    def show(ax, img, sl, cmap="gray", vmin=None, vmax=None, outline=None, color="r"):
        ax.imshow(img[sl].T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
        if outline is not None and outline[sl].any():
            ax.contour(outline[sl].T.astype(float), levels=[0.5], colors=color, linewidths=0.6)
        ax.set_xticks([])
        ax.set_yticks([])

    pr = lambda a, m=None: (np.percentile(a[m] if m is not None and m.any() else a, 1),
                            np.percentile(a[m] if m is not None and m.any() else a, 99))
    if kind == "gross":
        e4 = can(os.path.join(meg, f"{lab}_echo-4_part-mag_MEGRE.nii.gz"))
        chi_p = _find(lx, variant, code, ses, "Chimap")
        panels = [("mag e1", e1, "gray"), ("mag e4", e4, "gray")]
        for name, tok in (("total field", "field-ppm"), ("local field", "localfield")):
            p = _find(lx, variant, code, ses, tok)
            if p:
                panels.append((name, can(p), "gray"))
        chi = can(chi_p) if chi_p else None
        m3p = os.path.join(lx, "blind", "bids_all", "derivatives", "mask3", f"sub-{code}", f"ses-{ses}", "anat",
                           f"{lab}_mask.nii.gz")
        m3 = can(m3p) > 0.5 if os.path.exists(m3p) else None
        levels = [("dentate", zc(dn, 0.2)), ("midbrain", zc(sn, 0.45)), ("basal ganglia", zc(gp, 0.6))]
        ncol = len(panels) + (1 if chi is not None else 0)
        fig, axes = plt.subplots(3, ncol, figsize=(2.2 * ncol, 6.6), squeeze=False)
        for i, (lvl, z) in enumerate(levels):
            sl = (slice(None), slice(None), int(np.clip(z, 0, e1.shape[2] - 1)))
            for j, (name, img, cm) in enumerate(panels):
                if "field" in name:
                    v = np.nanpercentile(np.abs(img[m3]) if m3 is not None else np.abs(img), 99) or 1.0
                    show(axes[i, j], img, sl, cm, -v, v)
                else:
                    show(axes[i, j], img, sl, cm, *pr(img))
                if i == 0:
                    axes[i, j].set_title(name, fontsize=8)
            if chi is not None:
                show(axes[i, -1], np.nan_to_num(chi), sl, "gray", -0.2, 0.2, outline=m3, color="y")
                if i == 0:
                    axes[i, -1].set_title("chi (+/-0.2 ppm)", fontsize=8)
            axes[i, 0].set_ylabel(lvl, fontsize=8)
    else:
        t1 = can(os.path.join(rd, "t1_megre.nii.gz"))
        c = np.argwhere(sn).mean(0).round().astype(int) if sn.any() else (np.array(e1.shape) // 2)
        zd = zc(dn, 0.2)
        cuts = [("axial SN-VTA", (slice(None), slice(None), c[2])), ("coronal SN-VTA", (slice(None), c[1], slice(None))),
                ("sagittal SN-VTA", (c[0], slice(None), slice(None))), ("axial dentate", (slice(None), slice(None), zd))]
        fig, axes = plt.subplots(2, 4, figsize=(9, 4.8), squeeze=False)
        for j, (name, sl) in enumerate(cuts):
            ol = dn if "dentate" in name else sn
            show(axes[0, j], t1, sl, "gray", *pr(t1), outline=ol)
            show(axes[1, j], e1, sl, "gray", *pr(e1), outline=ol)
            axes[0, j].set_title(name, fontsize=8)
        axes[0, 0].set_ylabel("T1w", fontsize=8)
        axes[1, 0].set_ylabel("echo-1 mag", fontsize=8)
    fig.suptitle(f"{qid}  ({kind})", fontsize=9)
    fig.tight_layout()
    L.ensure_parent(out)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out, True


def cmd_montage(args, cfg):
    from concurrent.futures import ProcessPoolExecutor
    lx = L.lx_root()
    ids = qc_ids(lx, args.dry_run)
    L.say(f"montage {args.kind}: {len(ids)} sessions -> {QC(lx, 'montages', args.kind)}/<qc_id>.png", args.dry_run)
    if args.dry_run:
        return
    jobs = [(lx, args.kind, args.variant, c, s, q, cfg) for (c, s), q in ids.items()]
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        n = sum(made for _, made in ex.map(montage_one, jobs))
    L.say(f"{n} montages written")


def cmd_sheet(args):
    lx = L.lx_root()
    ids = qc_ids(lx, args.dry_run)
    out = QC(lx, f"rating_pass{args.pass_}.csv")
    rev = {}
    for f in glob.glob(os.path.join(lx, "flags_coded", "plaus_flags_qsmxt_primary.csv")) + \
            glob.glob(os.path.join(lx, "flags_coded", "g5_flags.csv")):
        for r in L.read_csv(f):
            k = (r["code"], int(r["ses_order"]))
            rev[k] = rev.get(k, 0) | int(r.get("forced_review", "0") == "1" or r.get("slab_review", "0") == "1")
    items = sorted(ids.items(), key=lambda kv: kv[1])
    rng = random.Random(L.SEED + args.pass_)
    if args.pass_ in (2, 4):
        ref = 1 if args.pass_ == 2 else 3
        closed = QC(lx, f"pass{ref}_closed_utc.txt")
        if not os.path.exists(closed) and not args.dry_run:
            sys.exit(f"pass {ref} is not closed (qc_montage.py finalize --close-pass{ref} writes the time stamp)")
        if os.path.exists(closed):
            t1 = datetime.datetime.strptime(open(closed).read().strip(), "%Y-%m-%dT%H:%M:%SZ")
            if (datetime.datetime.utcnow() - t1).days < 14 and not args.dry_run:
                sys.exit(f"re-rating must start >= 14 days after pass {ref} was closed (registration A.5)")
        items = rng.sample(items, max(1, round(args.fraction * len(items))))
    rng.shuffle(items)
    rows = [{"order": i + 1, "qc_id": q, "gross_png": f"montages/gross/{q}.png",
             "placement_png": f"montages/placement/{q}.png", "forced_review": rev.get(k, 0),
             "grade": "", "placement": "", "failure_codes": "", "notes": ""} for i, (k, q) in enumerate(items)]
    if os.path.exists(out):
        sys.exit(f"{out} exists")
    L.write_csv(out, rows, ["order", "qc_id", "gross_png", "placement_png", "forced_review", "grade", "placement",
                            "failure_codes", "notes"], args.dry_run)
    L.say(f"rating sheet pass {args.pass_}: {len(rows)} items -> {out}", args.dry_run)


def _read_sheet(p):
    rows = L.read_csv(p)
    bad = [r["qc_id"] for r in rows if r["grade"] not in ("0", "1", "2") or r["placement"] not in ("pass", "fail")]
    if bad:
        sys.exit(f"{p}: {len(bad)} items without a valid grade (0/1/2) and placement (pass/fail)")
    return {r["qc_id"]: r for r in rows}


def wkappa(pairs, k, weights):
    import numpy as np
    O = np.zeros((k, k))
    for x, y in pairs:
        O[x, y] += 1
    n = O.sum()
    i, j = np.indices((k, k))
    W = np.abs(i - j) / (k - 1) if weights == "linear" else ((i - j) / (k - 1)) ** 2 if weights == "quadratic" \
        else (i != j).astype(float)
    E = np.outer(O.sum(1), O.sum(0)) / n
    den = (W * E).sum()
    return (float(1 - (W * O).sum() / den) if den > 0 else float("nan")), O.astype(int).tolist()


def cmd_kappa(args):
    lx = L.lx_root()
    refined = os.path.exists(QC(lx, "rating_pass3.csv"))
    first, second = ("rating_pass3.csv", "rating_pass4.csv") if refined else ("rating_pass1.csv", "rating_pass2.csv")
    L.say(f"kappa: {first} vs {second} -> qc/kappa.json", args.dry_run)
    if args.dry_run:
        return
    a, b = _read_sheet(QC(lx, first)), _read_sheet(QC(lx, second))
    common = sorted(set(a) & set(b))
    g = [(int(a[q]["grade"]), int(b[q]["grade"])) for q in common]
    p = [(int(a[q]["placement"] == "fail"), int(b[q]["placement"] == "fail")) for q in common]
    kl, tab = wkappa(g, 3, "linear")
    kq, _ = wkappa(g, 3, "quadratic")
    kp, ptab = wkappa(p, 2, "nominal")
    res = {"n": len(common), "reference_pass": first, "kappa_linear": kl, "kappa_quadratic": kq,
           "kappa_placement": kp, "grade_table_rows_pass_ref": tab, "placement_table": ptab,
           "G4_pass": bool(kl == kl and kl >= 0.60), "decision_weight": "linear", "refined": refined,
           "rerating_pass": second, "utc": L.utcnow()}
    L.write_json(QC(lx, "kappa.json"), res)
    print(json.dumps(res, indent=1))


def cmd_finalize(args):
    lx = L.lx_root()
    for k in (1, 3):
        if getattr(args, f"close_pass{k}", False):
            p = QC(lx, f"pass{k}_closed_utc.txt")
            _read_sheet(QC(lx, f"rating_pass{k}.csv"))
            if not args.dry_run:
                with open(p, "w") as f:
                    f.write(L.utcnow() + "\n")
            L.say(f"pass {k} closed -> {p}", args.dry_run)
            return
    L.say("finalize -> qc/grades_final.csv", args.dry_run)
    if args.dry_run:
        return
    kj = QC(lx, "kappa.json")
    if not os.path.exists(kj):
        sys.exit("run kappa first")
    kap = json.load(open(kj))
    if not kap["G4_pass"] and not kap.get("refined"):
        sys.exit("G4 failed (kappa < 0.60): refine the rubric blind ONCE (log the new rubric hash), re-rate every "
                 "session (sheet --pass 3, finalize --close-pass3), >= 14 days later sheet --pass 4, then kappa")
    base = _read_sheet(QC(lx, kap["reference_pass"]))
    rer = _read_sheet(QC(lx, kap.get("rerating_pass", "rating_pass2.csv")))
    inv = {q: k for k, q in qc_ids(lx, False).items()}
    rows = []
    for q, r in base.items():
        g, pl, again = int(r["grade"]), r["placement"], 0
        codes = {x.strip().upper() for x in (r.get("failure_codes") or "").split(";") if x.strip()}
        if q in rer:     # conservative rule for the re-rated subset: worse grade, fail if either rating fails
            g = max(g, int(rer[q]["grade"]))
            pl = "fail" if "fail" in (pl, rer[q]["placement"]) else "pass"
            codes |= {x.strip().upper() for x in (rer[q].get("failure_codes") or "").split(";") if x.strip()}
            again = 1
        c, s = inv[q]
        rows.append({"code": c, "ses_order": s, "qc_grade": g, "placement": pl, "re_rated": again,
                     "failure_codes": ";".join(sorted(codes)), "dental_D": int("D" in codes)})
    rows.sort(key=lambda r: (r["code"], r["ses_order"]))
    out = QC(lx, "grades_final.csv")
    if os.path.exists(out):
        sys.exit(f"{out} exists")
    L.write_csv(out, rows, ["code", "ses_order", "qc_grade", "placement", "re_rated", "failure_codes", "dental_D"])
    L.write_json(QC(lx, "g4_final.json"), {"G4_pass": bool(kap["G4_pass"]), "refined": bool(kap.get("refined")),
                                           "kappa_linear": kap["kappa_linear"], "utc": L.utcnow()})
    L.say(f"wrote {out} ({len(rows)} sessions); G4 pass: {kap['G4_pass']}")


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = a.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("montage")
    m.add_argument("--kind", choices=["gross", "placement"], required=True)
    m.add_argument("--variant", default="primary")
    m.add_argument("--jobs", type=int, default=int(os.environ.get("LX_MAXJOBS", "8")))
    s = sub.add_parser("sheet")
    s.add_argument("--pass", dest="pass_", type=int, choices=[1, 2, 3, 4], required=True)
    s.add_argument("--fraction", type=float, default=0.20)
    sub.add_parser("kappa")
    f = sub.add_parser("finalize")
    f.add_argument("--close-pass1", action="store_true")
    f.add_argument("--close-pass3", action="store_true")
    for p in (m, s, f):
        p.add_argument("--dry-run", action="store_true")
    args, rest = a.parse_known_args()
    if not hasattr(args, "dry_run"):
        args.dry_run = "--dry-run" in rest
    cfg = L.load_config()
    {"montage": lambda: cmd_montage(args, cfg), "sheet": lambda: cmd_sheet(args), "kappa": lambda: cmd_kappa(args),
     "finalize": lambda: cmd_finalize(args)}[args.cmd]()


if __name__ == "__main__":
    main()
