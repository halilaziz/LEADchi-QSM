#!/usr/bin/env python3
"""build_input.py - POST HOC (2026-10-28): build the merged ROI table of the complete-human-QC sample (H_all180).

Reads (never writes) the frozen merged table mirror/unlocked/merged_qsmxt_primary_SNVTA.csv and the kit-2 scenario
table human_qc_kit/kit2_2026-10-28/analysis/real_run_20261028/kit2_scenario_usable.csv. Every cell is read and
written as TEXT, so the output differs from the frozen table only in the replaced columns:

  usable        <- usable__<scenario>
  first_usable  <- first_usable__<scenario>
  fail_reason   <- fail_reason__<scenario>
  qc_grade      <- the scenario's visual grade (H_all180: human_grade_all; AI_frozen: the frozen AI grade)
  (the frozen AI grade is kept in a new column qc_grade_ai for H_all180)

Check: rebuilding with scenario AI_frozen must reproduce the frozen table byte for byte (it is written to a scratch
path given by --check-out, never into the frozen tree).

  python3 build_input.py --scenario H_all180 --out ../input/merged_H_all180.csv
"""
import argparse, csv, hashlib, io, sys

LX = "/path/to/Penn_LEAD/LEADchi"
MERGED = f"{LX}/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv"
SCEN = f"{LX}/human_qc_kit/kit2_2026-10-28/analysis/real_run_20261028/kit2_scenario_usable.csv"


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def build(scenario):
    raw = open(MERGED, newline="").read()
    rd = list(csv.reader(io.StringIO(raw)))
    hdr, rows = rd[0], rd[1:]
    sc = {r["sess_key"]: r for r in csv.DictReader(open(SCEN, newline=""))}
    ip, iss = hdr.index("participant_id"), hdr.index("session_id")
    iu, ifu, ifr, iq = (hdr.index(c) for c in ("usable", "first_usable", "fail_reason", "qc_grade"))
    if len(sc) != len(rows):
        sys.exit("scenario table does not cover the merged table one-to-one")
    out_hdr = hdr + (["qc_grade_ai"] if scenario != "AI_frozen" else [])
    out = []
    for r in rows:
        k = f"{r[ip]}_{r[iss]}"
        s = sc[k]
        r2 = list(r)
        r2[iu] = str(int(float(s[f"usable__{scenario}"])))
        r2[ifu] = str(int(float(s[f"first_usable__{scenario}"])))
        r2[ifr] = s[f"fail_reason__{scenario}"] or "NA"         # frozen table writes a missing reason as NA
        if scenario == "H_all180":
            g = s["human_grade_all"]
            r2[iq] = "NA" if g == "" else repr(float(g))   # frozen style "1.0"; missing as NA
            r2.append(r[iq])
        elif scenario == "K1_H_replace":
            g = s["human_grade_kit1"] or s["ai_grade"]
            r2[iq] = "NA" if g == "" else repr(float(g)); r2.append(r[iq])
        out.append(r2)
    # keep the frozen table's text conventions (frozen first_usable/usable are written as 0/1; check below)
    buf = io.StringIO(newline="")
    w = csv.writer(buf, lineterminator="\n" if "\r\n" not in raw else "\r\n")
    w.writerow(out_hdr); w.writerows(out)
    txt = buf.getvalue()
    if not raw.endswith("\n"):
        txt = txt.rstrip("\n")                         # the frozen table has no final newline
    return txt, raw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="H_all180")
    ap.add_argument("--out", required=True)
    ap.add_argument("--check-out", default=None, help="scratch path for the AI_frozen byte-identity check")
    a = ap.parse_args()
    if a.check_out:
        txt, raw = build("AI_frozen")
        open(a.check_out, "w", newline="").write(txt)
        same = txt == raw
        print(f"AI_frozen rebuild byte-identical to the frozen merged table: {same} "
              f"(frozen sha256 {sha(MERGED)}, rebuild sha256 {sha(a.check_out)})")
        if not same:
            sys.exit("rebuild check failed")
    txt, _ = build(a.scenario)
    open(a.out, "w", newline="").write(txt)
    print(f"wrote {a.out} sha256 {sha(a.out)}")


if __name__ == "__main__":
    main()
