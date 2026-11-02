#!/usr/bin/env python3
"""build_kit.py - LEAD-chi blinded HUMAN visual-QC re-rating kit (POST HOC; built 2026-10-07).

Purpose: a QSM-experienced co-author, blind to group and to the AI grade, re-grades a selection of MEGRE sessions
with the registered rubric (qc_rubric.md sections 1-3) so that the AI rater of amendment 1c can be validated.
This is a post hoc validation; it changes no frozen file and no registered result.

Sample (seed 20261007, numpy.random.default_rng):
  1. every session with AI final grade 2 (grades_final.csv),
  2. a simple random sample of 40 sessions with AI final grade 0 or 1 (pool sorted by code, ses_order),
  3. 10 hidden duplicates drawn at random from the 1+2 set (intra-rater reliability).
Each item gets a new random kit ID ('h' + 6 hex, same rng); the order is shuffled with the same rng, and a shuffle is
redrawn until every duplicate pair is >= MIN_GAP positions apart (number of draws recorded).
Montages are copied from the coded montage folder; the matplotlib super-title that prints the original QC ID is
blanked and replaced by the kit ID (no other pixel is changed); PNG text metadata is dropped.

Outputs (all inside human_qc_kit/):
  rater_kit/images/{gross,placement}/<kit_id>.png, rater_kit/rating_sheet_blank.csv, analysis/items_for_html.json
  (embedded in rate.html by make_html.py), _key_DO_NOT_SEND/kit_key.csv, analysis/build_log.json
Reads only: mirror/qc/grades_final.csv, mirror/qc/qc_ids.csv, mirror/qc/rating_pass1_blank_template.csv (forced-review
flag; the blank template has no AI grades or notes) and the montage PNGs. Never reads the AI notes.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

SEED = 20261007
N_PASS_SAMPLE = 40
N_DUP = 10
MIN_GAP = 20

LX = "/path/to/Penn_LEAD/LEADchi"
KIT = os.path.join(LX, "human_qc_kit")
GRADES = os.path.join(LX, "mirror/qc/grades_final.csv")
QCIDS = os.path.join(LX, "mirror/qc/qc_ids.csv")
BLANK1 = os.path.join(LX, "mirror/qc/rating_pass1_blank_template.csv")
# Source montages: the permanent copy of the 180 + 180 PNGs the AI rater graded (SHA256SUMS in that folder).
# POST HOC CHANGE 2026-10-07 (review item H1): the kit was built on 2026-10-07 from the session-scratchpad copy of the
# same files (a temporary folder, no longer referenced); all 202 source files used by the kit are byte-identical to the
# archive (build_log.json src_sha256 == archive SHA256SUMS; checked by verify_kit.py).
MONT = os.path.join(LX, "qc_montages_archive", "montages")
FONT = None
try:
    import matplotlib
    FONT = os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans.ttf")
except Exception:  # pragma: no cover
    pass


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def title_band(arr):
    """Rows of the figure super-title: the first block of non-white rows, which must be narrow and centred."""
    rgb = arr[:, :, :3].astype(int)
    nonwhite = (rgb < 245).any(axis=2)
    rows = np.where(nonwhite.any(axis=1))[0]
    if not len(rows):
        raise ValueError("blank image")
    y0 = rows[0]
    y1 = y0
    while y1 + 1 < arr.shape[0] and nonwhite[y1 + 1].any():
        y1 += 1
    cols = np.where(nonwhite[y0:y1 + 1].any(axis=0))[0]
    W = arr.shape[1]
    ok = (y1 - y0 + 1) <= 25 and (cols[-1] - cols[0]) < 0.5 * W and abs((cols[0] + cols[-1]) / 2 - W / 2) < 0.15 * W
    nxt = rows[rows > y1]
    gap = int(nxt[0] - y1 - 1) if len(nxt) else 999
    if not ok or gap < 3:
        raise ValueError(f"title band not recognised (y {y0}-{y1}, x {cols[0]}-{cols[-1]}, gap {gap})")
    return int(y0), int(y1), gap


def relabel(src, dst, kit_id, kind):
    im = Image.open(src).convert("RGB")
    arr = np.asarray(im)
    y0, y1, gap = title_band(arr)
    a0, a1 = max(0, y0 - 2), y1 + min(2, gap - 1)
    arr = arr.copy()
    arr[a0:a1 + 1, :, :] = 255
    out = Image.fromarray(arr)
    d = ImageDraw.Draw(out)
    font = ImageFont.truetype(FONT, 14) if FONT and os.path.exists(FONT) else ImageFont.load_default()
    text = f"{kit_id}  ({kind})"
    bb = d.textbbox((0, 0), text, font=font)
    d.text(((out.width - (bb[2] - bb[0])) / 2, y0 - 1), text, fill=(0, 0, 0), font=font)
    out.save(dst, format="PNG", optimize=True)          # no pnginfo -> no text chunks
    return {"title_rows": [y0, y1], "blanked_rows": [a0, a1]}


def main():
    grades = read_csv(GRADES)
    qmap = {(r["code"], int(r["ses_order"])): r["qc_id"] for r in read_csv(QCIDS)}
    frev = {r["qc_id"]: int(r["forced_review"]) for r in read_csv(BLANK1)}
    for r in grades:
        r["ses_order"] = int(r["ses_order"])
        r["qc_grade"] = int(r["qc_grade"])
        r["qc_id"] = qmap[(r["code"], r["ses_order"])]
        for k in ("gross", "placement"):
            if not os.path.exists(os.path.join(MONT, k, r["qc_id"] + ".png")):
                sys.exit(f"missing montage {k}/{r['qc_id']}.png")
    grades.sort(key=lambda r: (r["code"], r["ses_order"]))
    fail = [r for r in grades if r["qc_grade"] == 2]
    pool = [r for r in grades if r["qc_grade"] in (0, 1)]
    rng = np.random.default_rng(SEED)
    pick = sorted(rng.choice(len(pool), size=N_PASS_SAMPLE, replace=False).tolist())
    passed = [pool[i] for i in pick]
    selected = [dict(r, stratum="AI_fail_all") for r in fail] + [dict(r, stratum="AI_pass_sample") for r in passed]
    dup_idx = sorted(rng.choice(len(selected), size=N_DUP, replace=False).tolist())
    items = [dict(r, dup_of_index=None) for r in selected] + [dict(selected[i], dup_of_index=i) for i in dup_idx]
    ids = set()
    for it in items:
        while True:
            k = "h" + format(int(rng.integers(0, 16 ** 6)), "06x")
            if k not in ids:
                ids.add(k)
                it["kit_id"] = k
                break
    for it in items:
        it["is_duplicate_of"] = items[it["dup_of_index"]]["kit_id"] if it["dup_of_index"] is not None else ""
    draws = 0
    while True:
        draws += 1
        perm = rng.permutation(len(items))
        pos = {int(p): i for i, p in enumerate(perm)}           # item index -> position
        gaps = [abs(pos[j] - pos[it["dup_of_index"]]) for j, it in enumerate(items) if it["dup_of_index"] is not None]
        if min(gaps) >= MIN_GAP:
            break
    ordered = [items[int(p)] for p in perm]
    for o, it in enumerate(ordered, 1):
        it["order"] = o

    rk = os.path.join(KIT, "rater_kit")
    log = {"seed": SEED, "rng": "numpy.random.default_rng", "numpy": np.__version__, "n_ai_fail": len(fail),
           "n_ai_pass_pool": len(pool), "n_ai_pass_sampled": len(passed), "n_duplicates": N_DUP,
           "n_items": len(ordered), "min_duplicate_gap": MIN_GAP, "shuffle_draws": draws,
           "duplicate_gaps": sorted(gaps), "images": {}, "inputs": {p: sha256(p) for p in (GRADES, QCIDS, BLANK1)}}
    for it in ordered:
        for kind in ("gross", "placement"):
            src = os.path.join(MONT, kind, it["qc_id"] + ".png")
            dst = os.path.join(rk, "images", kind, it["kit_id"] + ".png")
            if os.path.exists(dst):
                sys.exit(f"{dst} exists; the kit is built once (delete nothing, start a new folder)")
            info = relabel(src, dst, it["kit_id"], kind)
            log["images"][f"{kind}/{it['kit_id']}"] = dict(info, src_sha256=sha256(src), dst_sha256=sha256(dst))

    # rater-facing item list (no AI information)
    with open(os.path.join(rk, "rating_sheet_blank.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["order", "kit_id", "forced_review", "grade", "failure_codes", "placement", "note", "timestamp"])
        for it in ordered:
            w.writerow([it["order"], it["kit_id"], frev[it["qc_id"]], "", "", "", "", ""])
    with open(os.path.join(KIT, "analysis", "items_for_html.json"), "w") as f:
        json.dump([{"order": it["order"], "kit_id": it["kit_id"], "forced_review": frev[it["qc_id"]]} for it in ordered],
                  f, indent=0)

    # sealed key
    w_pass = len(pool) / len(passed)
    keyp = os.path.join(KIT, "_key_DO_NOT_SEND", "kit_key.csv")
    if os.path.exists(keyp):
        sys.exit(f"{keyp} exists")
    with open(keyp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["kit_id", "order", "code", "ses_order", "qc_id", "ai_grade", "ai_placement", "ai_failure_codes",
                    "ai_dental_D", "forced_review", "stratum", "sampling_weight", "is_duplicate_of"])
        for it in ordered:
            w.writerow([it["kit_id"], it["order"], it["code"], it["ses_order"], it["qc_id"], it["qc_grade"],
                        it["placement"], it["failure_codes"], it["dental_D"], frev[it["qc_id"]], it["stratum"],
                        1.0 if it["stratum"] == "AI_fail_all" else round(w_pass, 6), it["is_duplicate_of"]])
    log["kit_key_sha256"] = sha256(keyp)
    log["sampling_weight_ai_pass"] = w_pass
    with open(os.path.join(KIT, "analysis", "build_log.json"), "w") as f:
        json.dump(log, f, indent=1)
    print(json.dumps({k: v for k, v in log.items() if k != "images"}, indent=1))


if __name__ == "__main__":
    main()
