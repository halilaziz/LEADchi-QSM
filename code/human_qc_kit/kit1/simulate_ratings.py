#!/usr/bin/env python3
"""simulate_ratings.py - SIMULATED human ratings for an end-to-end test of analyze_human_qc.py.

THE OUTPUT IS SIMULATED. It is generated from the AI grades in kit_key.csv plus random disagreement and is used only
to test the pipeline (analysis/test_run/). It must never be analysed or reported as real ratings.
Format = the CSV that rate.html exports (kit_id, grade, failure_codes, placement, note, timestamp; CRLF).
One item is left unrated on purpose to exercise the missing-rating path.
"""
import csv
import datetime
import os

import numpy as np

KIT = "/path/to/Penn_LEAD/LEADchi/human_qc_kit"
OUT = os.path.join(KIT, "analysis", "test_run", "SIMULATED_ratings.csv")
rng = np.random.default_rng(777)
P = {2: [0.03, 0.17, 0.80], 1: [0.12, 0.78, 0.10], 0: [0.70, 0.30, 0.0]}

key = list(csv.DictReader(open(os.path.join(KIT, "_key_DO_NOT_SEND", "kit_key.csv"))))
by_id = {r["kit_id"]: r for r in key}


def draw(r):
    g = int(rng.choice(3, p=P[int(r["ai_grade"])]))
    ai = [c for c in r["ai_failure_codes"].split(";") if c]
    codes = set()
    if g == 2:
        codes = set(c for c in ai if rng.random() < 0.85) or {str(rng.choice(list("MWB")))}
        if rng.random() < 0.15:
            codes.add(str(rng.choice(list("MWSBKA"))))
    if "D" in ai and rng.random() < 0.8:
        codes.add("D")
    pl = "fail" if rng.random() < 0.03 else "pass"
    note = "SIMULATED" + ("; outline displaced (simulated)" if pl == "fail" else "")
    return {"grade": str(g), "failure_codes": ";".join(sorted(codes)), "placement": pl, "note": note}


t0 = datetime.datetime(2026, 10, 7, 12, 0, tzinfo=datetime.timezone.utc)
rows, done = [], {}
for k, r in enumerate(sorted(key, key=lambda r: int(r["order"]))):
    if r["is_duplicate_of"] and r["is_duplicate_of"] in done and rng.random() < 0.8:
        x = dict(done[r["is_duplicate_of"]])
    else:
        x = draw(r)
    done[r["kit_id"]] = x
    ts = (t0 + datetime.timedelta(seconds=45 * k)).strftime("%Y-%m-%dT%H:%M:%S.000Z")
    rows.append([r["kit_id"], x["grade"], x["failure_codes"], x["placement"], x["note"], ts])
rows[57][1:] = ["", "", "", "", ""]          # one unrated item (tests the missing-rating path)
with open(OUT, "w", newline="") as f:
    w = csv.writer(f, lineterminator="\r\n")
    w.writerow(["kit_id", "grade", "failure_codes", "placement", "note", "timestamp"])
    w.writerows(rows)
print(f"wrote {OUT} (SIMULATED, {len(rows)} rows)")
