#!/usr/bin/env python3
"""assemble_sheet.py <workflow_output.json> <sheet_template.csv> <out.csv> : fill grade/placement/failure_codes/notes from
the AI rater output into the server's rating sheet (same columns/order). Rows without a montage stay blank."""
import csv, json, sys
d = json.load(open(sys.argv[1])); res = d["result"] if isinstance(d, dict) and "result" in d else d
G = {}
for b in res:
    if not b or not b.get("result"): continue
    for g in b["result"]["grades"]:
        if g["qc_id"] in G: sys.exit(f"duplicate grade for {g['qc_id']}")
        G[g["qc_id"]] = g
rows = [r for r in csv.DictReader(l for l in open(sys.argv[2]) if not l.startswith("rating sheet"))]
n = 0
for r in rows:
    g = G.get(r["qc_id"])
    if g:
        r["grade"] = str(int(g["grade"])); r["placement"] = g["placement"]
        r["failure_codes"] = ";".join(sorted({c.strip().upper() for c in g["failure_codes"].split(";") if c.strip()}))
        r["notes"] = g["notes"].replace("\n", " ").strip(); n += 1
with open(sys.argv[3], "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["order", "qc_id", "gross_png", "placement_png", "forced_review", "grade", "placement", "failure_codes", "notes"])
    w.writeheader(); w.writerows(rows)
missing = [r["qc_id"] for r in rows if not r["grade"]]
print(f"filled {n}/{len(rows)}; ungraded: {missing}")
from collections import Counter
print("grade counts", Counter(r["grade"] for r in rows if r["grade"]), "placement", Counter(r["placement"] for r in rows if r["placement"]))
print("codes", Counter(c for r in rows for c in r["failure_codes"].split(";") if c))
