#!/usr/bin/env python3
"""build_kit2.py - LEAD-chi blinded human visual-QC kit 2 (POST HOC; built 2026-10-28; plan: ../PLAN_kit2.md).

Kit 2 gives the kit-1 rater the AI-pass sessions that kit 1 did not hold, so that every one of the 180 AI-graded
sessions gets a human grade, plus hidden repeats of kit-1 items (test-retest; they also hide the composition).

Sample (seed 20261028, numpy.random.default_rng; one stream, used in this order):
  0. new items  = every AI grade-0/1 session (grades_final.csv) that is not a kit-1 primary item (79; no sampling)
  1. 11 repeats drawn from the 61 kit-1 primary AI-fail items  (sorted by code, ses_order)
  2. 10 repeats drawn from the 40 kit-1 primary AI-pass items  (sorted by code, ses_order)
  3. a new kit ID per item: 'h' + 6 hex, unique, and different from every kit-1 kit ID
  4. shuffle; redrawn until the repeats are spread: each block of 20 consecutive positions (1-20, ..., 81-100) holds
     between SPREAD_MIN and SPREAD_MAX repeats (number of draws recorded).
Image treatment is kit 1's, called unchanged: human_qc_kit/analysis/build_kit.py relabel() (the matplotlib title that
prints the QC ID is blanked and replaced by the kit ID; no other pixel changes; no PNG text metadata). The forced-review
flag comes from the same source as kit 1 (mirror/qc/rating_pass1_blank_template.csv; it holds no AI grade).

Reads only (all read-only): mirror/qc/grades_final.csv, mirror/qc/qc_ids.csv, mirror/qc/rating_pass1_blank_template.csv,
human_qc_kit/_key_DO_NOT_SEND/kit_key.csv, qc_montages_archive/montages/*. Writes only inside kit2_2026-10-28/ and
refuses to overwrite anything.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True                       # importing kit-1 build_kit must not write a __pycache__ there
import numpy as np  # noqa: E402

K1 = "/path/to/Penn_LEAD/LEADchi/human_qc_kit"
K2 = os.path.join(K1, "kit2_2026-10-28")
sys.path.insert(0, os.path.join(K1, "analysis"))
import build_kit as B1  # noqa: E402  (relabel(), paths; main() is not run)

SEED = 20261028
N_REP_FAIL = 11
N_REP_PASS = 10
BLOCK = 20
SPREAD_MIN, SPREAD_MAX = 3, 6
KEY1 = os.path.join(K1, "_key_DO_NOT_SEND", "kit_key.csv")
RK = os.path.join(K2, "rater_kit2")
KEY2 = os.path.join(K2, "_key_DO_NOT_SEND", "kit2_key.csv")
LOG = os.path.join(K2, "analysis", "build_log_kit2.json")
ITEMS = os.path.join(K2, "analysis", "items_for_html.json")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def select(grades, key1, qmap):
    """Deterministic selection (shared with verify_kit2.py). Returns (items in final order, info dict)."""
    for r in grades:
        r["ses_order"] = int(r["ses_order"])
        r["qc_grade"] = int(r["qc_grade"])
        r["qc_id"] = qmap[(r["code"], r["ses_order"])]
    grades = sorted(grades, key=lambda r: (r["code"], r["ses_order"]))
    byq = {r["qc_id"]: r for r in grades}
    prim1 = [r for r in key1 if not r["is_duplicate_of"]]
    k1q = {r["qc_id"] for r in prim1}
    new = [dict(r, stratum="AI_pass_new", kit1_id="") for r in grades if r["qc_grade"] in (0, 1) and r["qc_id"] not in k1q]
    srt = lambda rows: sorted(rows, key=lambda r: (r["code"], int(r["ses_order"])))
    p_fail = srt([r for r in prim1 if r["stratum"] == "AI_fail_all"])
    p_pass = srt([r for r in prim1 if r["stratum"] == "AI_pass_sample"])
    rng = np.random.default_rng(SEED)
    i_f = sorted(rng.choice(len(p_fail), size=N_REP_FAIL, replace=False).tolist())
    i_p = sorted(rng.choice(len(p_pass), size=N_REP_PASS, replace=False).tolist())
    reps = ([dict(byq[p_fail[i]["qc_id"]], stratum="repeat_of_kit1_AI_fail", kit1_id=p_fail[i]["kit_id"]) for i in i_f] +
            [dict(byq[p_pass[i]["qc_id"]], stratum="repeat_of_kit1_AI_pass", kit1_id=p_pass[i]["kit_id"]) for i in i_p])
    items = new + reps
    used = {r["kit_id"] for r in key1}
    ids = set()
    for it in items:
        while True:
            k = "h" + format(int(rng.integers(0, 16 ** 6)), "06x")
            if k not in ids and k not in used:
                ids.add(k)
                it["kit2_id"] = k
                break
    draws = 0
    while True:
        draws += 1
        perm = rng.permutation(len(items))
        isrep = np.array([bool(items[int(p)]["kit1_id"]) for p in perm])
        per_block = [int(isrep[b:b + BLOCK].sum()) for b in range(0, len(items), BLOCK)]
        if all(SPREAD_MIN <= c <= SPREAD_MAX for c in per_block):
            break
    ordered = [items[int(p)] for p in perm]
    for o, it in enumerate(ordered, 1):
        it["order"] = o
    info = {"n_new": len(new), "n_rep_fail": len(i_f), "n_rep_pass": len(i_p), "shuffle_draws": draws,
            "repeats_per_block_of_20": per_block, "repeat_positions": [o for o, it in enumerate(ordered, 1) if it["kit1_id"]]}
    return ordered, info


def main():
    for p in (KEY2, LOG, ITEMS, os.path.join(RK, "rating_sheet_blank.csv")):
        if os.path.exists(p):
            sys.exit(f"{p} exists; kit 2 is built once (nothing is overwritten)")
    grades = read_csv(B1.GRADES)
    qmap = {(r["code"], int(r["ses_order"])): r["qc_id"] for r in read_csv(B1.QCIDS)}
    frev = {r["qc_id"]: int(r["forced_review"]) for r in read_csv(B1.BLANK1)}
    key1 = read_csv(KEY1)
    ordered, info = select(grades, key1, qmap)
    assert info["n_new"] == 79 and len(ordered) == 100, info
    log = {"seed": SEED, "rng": "numpy.random.default_rng", "numpy": np.__version__, **info,
           "n_items": len(ordered), "spread_rule": f"each block of {BLOCK} positions holds {SPREAD_MIN}-{SPREAD_MAX} repeats",
           "image_treatment": "human_qc_kit/analysis/build_kit.py relabel(), unchanged",
           "inputs": {p: sha256(p) for p in (B1.GRADES, B1.QCIDS, B1.BLANK1, KEY1,
                                              os.path.join(os.path.dirname(B1.MONT), "SHA256SUMS"),
                                              os.path.join(K1, "analysis", "build_kit.py"))},
           "images": {}}
    for kind in ("gross", "placement"):
        os.makedirs(os.path.join(RK, "images", kind), exist_ok=True)
    for it in ordered:
        for kind in ("gross", "placement"):
            src = os.path.join(B1.MONT, kind, it["qc_id"] + ".png")
            dst = os.path.join(RK, "images", kind, it["kit2_id"] + ".png")
            if os.path.exists(dst):
                sys.exit(f"{dst} exists")
            inf = B1.relabel(src, dst, it["kit2_id"], kind)
            log["images"][f"{kind}/{it['kit2_id']}"] = dict(inf, src_sha256=sha256(src), dst_sha256=sha256(dst))
    with open(os.path.join(RK, "rating_sheet_blank.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["order", "kit_id", "forced_review", "grade", "failure_codes", "placement", "note", "timestamp"])
        for it in ordered:
            w.writerow([it["order"], it["kit2_id"], frev[it["qc_id"]], "", "", "", "", ""])
    with open(ITEMS, "w") as f:
        json.dump([{"order": it["order"], "kit_id": it["kit2_id"], "forced_review": frev[it["qc_id"]]} for it in ordered],
                  f, indent=0)
    with open(KEY2, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["kit2_id", "order", "code", "ses_order", "qc_id", "ai_grade", "ai_placement", "ai_failure_codes",
                    "ai_dental_D", "forced_review", "stratum", "is_kit1_repeat", "kit1_id"])
        for it in ordered:
            w.writerow([it["kit2_id"], it["order"], it["code"], it["ses_order"], it["qc_id"], it["qc_grade"],
                        it["placement"], it["failure_codes"], it["dental_D"], frev[it["qc_id"]], it["stratum"],
                        int(bool(it["kit1_id"])), it["kit1_id"]])
    log["kit2_key_sha256"] = sha256(KEY2)
    with open(LOG, "w") as f:
        json.dump(log, f, indent=1)
    print(json.dumps({k: v for k, v in log.items() if k != "images"}, indent=1))


if __name__ == "__main__":
    main()
