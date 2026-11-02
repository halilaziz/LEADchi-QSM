#!/usr/bin/env python3
"""dwiqc.py - coded DWI bundle QC for the S3 sample (registration Q13, A.3 item 8; qc_rubric.md section 6).

Stage `setup`: copies the pyAFQ per-session bundle figures of SCP_R and MCP_R (--pattern, VERIFY at install) under a
  separate random DWI code set into $LX/dwiqc/figures/, writes the plaintext code<->ID file $LX/dwiqc/dwiqc_key.csv
  (its SHA-256 is registered; it is not opened until the freeze) and the rating sheet $LX/dwiqc/dwiqc_sheet.csv
  (dqc_id, scp_r_pass, mcp_r_pass, notes) in random order (seed 20260924).
Stage `freeze`: maps the completed sheet back to real IDs -> $LX/dwi/bundle_qc.csv (participant_id, session_id,
  scp_r_pass, mcp_r_pass) and prints its SHA-256 (frozen; read by dwi_geometry.py --stage table).
The rater is the first author and is not blind in practice (DWI group results were seen; declared). This code set
links nothing to chi.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import random
import re
import secrets
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

DEFAULT_PATTERN = "{pyafq}/sub-*/ses-*/**/*desc-{label}*.png"


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--stage", choices=["setup", "freeze"], required=True)
    a.add_argument("--pyafq", default=os.path.join(L.pl_root(), "pyAFQ"))
    a.add_argument("--pattern", default=DEFAULT_PATTERN)
    a.add_argument("--dry-run", action="store_true")
    args = a.parse_args()
    lx = L.lx_root()
    d = os.path.join(lx, "dwiqc")
    key_p, sheet_p = os.path.join(d, "dwiqc_key.csv"), os.path.join(d, "dwiqc_sheet.csv")
    if args.stage == "setup":
        figs = {}
        for label, name in (("SCPR", "scp"), ("MCPR", "mcp")):
            for f in sorted(glob.glob(args.pattern.format(pyafq=args.pyafq, label=label), recursive=True)):
                m = re.search(r"(sub-[A-Za-z0-9]+)_(ses-[A-Za-z0-9]+)", os.path.basename(f))
                if m:
                    figs.setdefault(m.groups(), {})[name] = f
        L.say(f"dwiqc setup: {len(figs)} DWI sessions with bundle figures -> {d}", args.dry_run)
        if args.dry_run:
            return
        if os.path.exists(key_p):
            sys.exit(f"{key_p} exists (set up once)")
        taken, rows, sheet = set(), [], []
        for (sub, ses), fs in sorted(figs.items()):
            while True:
                c = "d" + secrets.token_hex(4)
                if c not in taken:
                    taken.add(c)
                    break
            for name, f in fs.items():
                dst = os.path.join(d, "figures", f"{c}_{name}{os.path.splitext(f)[1]}")
                L.ensure_parent(dst)
                shutil.copyfile(f, dst)
            rows.append({"dqc_id": c, "participant_id": sub, "session_id": ses})
            sheet.append({"dqc_id": c, "scp_r_pass": "", "mcp_r_pass": "", "notes": ""})
        random.Random(L.SEED).shuffle(sheet)
        L.write_csv(key_p, rows, ["dqc_id", "participant_id", "session_id"])
        L.write_csv(sheet_p, sheet, ["dqc_id", "scp_r_pass", "mcp_r_pass", "notes"])
        print(f"key {key_p} sha256 {L.sha256_file(key_p)} (register; do not open before the freeze)")
        return
    out = os.path.join(lx, "dwi", "bundle_qc.csv")
    if os.path.exists(out):
        sys.exit(f"{out} exists (frozen once)")
    sheet = L.read_csv(sheet_p)
    bad = [r["dqc_id"] for r in sheet if r["scp_r_pass"] not in ("0", "1") or r["mcp_r_pass"] not in ("0", "1")]
    if bad:
        sys.exit(f"{len(bad)} sheet rows without 0/1 ratings")
    key = {r["dqc_id"]: r for r in L.read_csv(key_p)}
    rows = [{"participant_id": key[r["dqc_id"]]["participant_id"], "session_id": key[r["dqc_id"]]["session_id"],
             "scp_r_pass": r["scp_r_pass"], "mcp_r_pass": r["mcp_r_pass"]} for r in sheet]
    L.write_csv(out, sorted(rows, key=lambda r: (r["participant_id"], r["session_id"])),
                ["participant_id", "session_id", "scp_r_pass", "mcp_r_pass"])
    print(f"frozen {out} sha256 {L.sha256_file(out)}")


if __name__ == "__main__":
    main()
