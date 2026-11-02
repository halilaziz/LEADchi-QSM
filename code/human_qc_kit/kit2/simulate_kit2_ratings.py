#!/usr/bin/env python3
"""simulate_kit2_ratings.py - SIMULATED kit-1 and kit-2 rating files for testing analyze_kit2.py ONLY.

Never real ratings. The grades are the AI grade with random disagreement (seed fixed), in the rate.html export format
(kit_id, grade, failure_codes, placement, note, timestamp, rater). Every file name starts with SIMULATED_.

  python3 simulate_kit2_ratings.py --outdir test_run_SIMULATED/inputs [--seed 777] [--blank 1] [--p-flip 0.25]
"""
import argparse
import csv
import os

import numpy as np

K1 = "/path/to/Penn_LEAD/LEADchi/human_qc_kit"
K2 = os.path.join(K1, "kit2_2026-10-28")


def sim(key, idcol, rng, p_flip, p_place_fail, blank, out):
    rows = list(csv.DictReader(open(key, newline="")))
    blank_ids = set(rng.choice([r[idcol] for r in rows], size=blank, replace=False).tolist()) if blank else set()
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["kit_id", "grade", "failure_codes", "placement", "note", "timestamp", "rater"])
        for r in rows:
            if r[idcol] in blank_ids:
                w.writerow([r[idcol], "", "", "", "SIMULATED blank", "", "SIM"])
                continue
            g = int(r["ai_grade"])
            if rng.random() < p_flip:
                g = int(np.clip(g + rng.choice([-1, 1]), 0, 2))
            pl = "fail" if rng.random() < p_place_fail else "pass"
            w.writerow([r[idcol], g, "M" if g == 2 else "", pl, "SIMULATED", "2026-10-28T00:00:00.000Z", "SIM"])
    return sorted(blank_ids)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=777)
    ap.add_argument("--blank", type=int, default=1, help="kit-2 items left blank")
    ap.add_argument("--p-flip", dest="p_flip", type=float, default=0.25)
    ap.add_argument("--p-place-fail", dest="p_place_fail", type=float, default=0.02)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    o1 = os.path.join(a.outdir, f"SIMULATED_kit1_ratings_seed{a.seed}.csv")
    o2 = os.path.join(a.outdir, f"SIMULATED_kit2_ratings_seed{a.seed}_blank{a.blank}.csv")
    for o in (o1, o2):
        if os.path.exists(o):
            raise SystemExit(f"{o} exists")
    sim(os.path.join(K1, "_key_DO_NOT_SEND", "kit_key.csv"), "kit_id", rng, a.p_flip, a.p_place_fail, 0, o1)
    b = sim(os.path.join(K2, "_key_DO_NOT_SEND", "kit2_key.csv"), "kit2_id", rng, a.p_flip, a.p_place_fail, a.blank, o2)
    print(f"wrote {o1}\nwrote {o2} (blank kit-2 items: {b})")


if __name__ == "__main__":
    main()
