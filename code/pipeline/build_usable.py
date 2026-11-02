#!/usr/bin/env python3
"""build_usable.py - frozen, code-keyed first-usable list and the freeze manifest (task P5b, 2026-09-25).

Runs on server before unlock (v1.2 sec. 2.8 step 1). Reads ONLY visible, label-free pass/fail flags:
  blind/session_flags.csv                      (b) same-session T1w present and PennLINC T1 QC not "fail"
  qc/grades_final.csv                          (c) blinded MEGRE QC grade <= 1; (d) SN-VTA placement pass
  flags_coded/plaus_flags_<engine>_<variant>   (c) automatic grade 2; (c') QSM plausibility (SN-VTA or SN-only)
  flags_coded/g5_flags.csv                     (e) both chi-independent registration checks (G5)
(a) complete 4-echo MEGRE holds for every coded session by construction (blind_setup.py).

Stage `usable` writes, for every available branch (engine x pipeline x roi_def):
  freeze/usable_<engine>_<pipeline>_<roi_def>.csv   code, ses_order, usable, fail_reason, first_usable
  freeze/usable_summary.json                        counts only
  pipeline 'primary' = variant primary; 'fallback_PDF_iLSQR' = variant mv_pdf_ilsqr (sec. 14 PC fallback).
  roi_def 'SNVTA' uses the SN-VTA plausibility flag for (c'); 'SN_only' (G3 fallback) uses the SN-only flag.
  The visual grades (c, d) are the same for every branch (graded on the primary montages; declared).
Stage `manifest` writes freeze/freeze_manifest.json: SHA-256 of every registered pre-unlock artefact (coded tables,
  flags, QC sheets and kappa, usable lists, setup/lock manifests and ciphertexts, frozen configs, the pipeline and
  analysis code given with --include, the gate JSON) and the head hash of the run log. Posted in amendment 2.
--dry-run prints the plan only.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

BRANCH_VARIANT = {"primary": "primary", "fallback_PDF_iLSQR": "mv_pdf_ilsqr"}


def key(r):
    return (r["code"], int(r["ses_order"]))


def stage_usable(args):
    lx = L.lx_root()
    inputs = {"flags": os.path.join(lx, "blind", "session_flags.csv"),
              "grades": os.path.join(lx, "qc", "grades_final.csv"),
              "g5": os.path.join(lx, "flags_coded", "g5_flags.csv")}
    L.forbid_label_inputs(inputs.values())
    L.say(f"stage usable: inputs {json.dumps(inputs)}", args.dry_run)
    branches = []
    for eng in ("qsmxt", "sepia"):
        for pipe, var in BRANCH_VARIANT.items():
            p = os.path.join(lx, "flags_coded", f"plaus_flags_{eng}_{var}.csv")
            if os.path.exists(p) or args.dry_run:
                for rd in ("SNVTA", "SN_only"):
                    branches.append((eng, pipe, rd, p))
    L.say(f"branches: {[b[:3] for b in branches]}", args.dry_run)
    if args.dry_run:
        return
    for p in inputs.values():
        if not os.path.exists(p):
            sys.exit(f"{p} missing")
    sflags = {key(r): r for r in L.read_csv(inputs["flags"])}
    grades = {key(r): r for r in L.read_csv(inputs["grades"])}
    g5 = {key(r): r for r in L.read_csv(inputs["g5"])}
    missing = [k for k in sflags if k not in grades and sflags[k]["t1_ok"] == "1"]
    if missing:
        sys.exit(f"{len(missing)} registrable coded sessions have no final QC grade; grade them before freezing")
    rr = L.RunRecord(__file__, False)
    for p in inputs.values():
        rr.add_input(p)
    summary = {"utc": L.utcnow(), "branches": {}}
    for eng, pipe, rd, pf in branches:
        pl = {key(r): r for r in L.read_csv(pf)}
        rr.add_input(pf)
        rows = []
        for k in sorted(sflags):
            s, g, q, e = sflags[k], grades.get(k, {}), pl.get(k), g5.get(k, {})
            why = []
            if s["t1_ok"] != "1":
                why.append("b_t1w")
            if e.get("registered") != "1":
                why.append("e_not_registered")
            if q is None or q.get("status") != "ok":
                why.append("no_qsm")
            if g.get("qc_grade") not in ("0", "1"):
                why.append("c_qc_grade")
            if q is not None and q.get("auto_grade2") == "1":
                why.append("c_auto_grade2")
            if q is not None and q.get("snvta_fail" if rd == "SNVTA" else "snonly_fail") == "1":
                why.append("cprime_plausibility")
            if g.get("placement") != "pass":
                why.append("d_placement")
            if e.get("g5_pass") != "1":
                why.append("e_g5")
            rows.append({"code": k[0], "ses_order": k[1], "usable": int(not why), "fail_reason": ";".join(why)})
        first = {}
        for r in rows:
            if r["usable"] and (r["code"] not in first or r["ses_order"] < first[r["code"]]):
                first[r["code"]] = r["ses_order"]
        for r in rows:
            r["first_usable"] = int(first.get(r["code"]) == r["ses_order"])
        out = os.path.join(lx, "freeze", f"usable_{eng}_{pipe}_{rd}.csv")
        if os.path.exists(out):
            sys.exit(f"{out} exists (the list is frozen once)")
        L.write_csv(out, rows, ["code", "ses_order", "usable", "fail_reason", "first_usable"])
        rr.add_output(out)
        summary["branches"][f"{eng}_{pipe}_{rd}"] = {
            "sessions": len(rows), "usable": sum(r["usable"] for r in rows), "subjects_first_usable": len(first),
            "returners_with_2_usable": sum(1 for c in first if sum(1 for r in rows if r["code"] == c and r["usable"]) >= 2),
            "fail_reason_counts": {w: sum(w in r["fail_reason"].split(";") for r in rows) for w in
                                   ("b_t1w", "e_not_registered", "no_qsm", "c_qc_grade", "c_auto_grade2",
                                    "cprime_plausibility", "d_placement", "e_g5")}}
    L.write_json(os.path.join(lx, "freeze", "usable_summary.json"), summary)
    rr.finish()
    print(json.dumps(summary, indent=1))


def stage_manifest(args):
    lx = L.lx_root()
    pats = ["blind/setup_manifest.json", "blind/session_flags.csv", "blind/key.csv.enc", "blind/key.sessionkey.enc",
            "blind/coded_files.sha256", "blind/processing_order.txt", "blind/pilot10.txt", "blind/lcurve2.txt",
            "blind/gate20.txt", "locked/labels_locked.tar.enc", "locked/labels_locked.sessionkey.enc",
            "locked/labels_locked_manifest.json", "prereg/*.toml", "prereg/*.json", "prereg/*.txt",
            "prereg/*.sha256", "prereg/lock_pub.pem", "tables_coded/*.csv", "flags_coded/*", "freeze/usable_*",
            "qc/rating_pass*.csv", "qc/grades_final.csv", "qc/kappa.json", "qc/qc_ids.csv", "scripts/pipeline/*"]
    files = []
    for p in pats:
        files += sorted(glob.glob(os.path.join(lx, p)))
    for inc in args.include or []:
        files += sorted(glob.glob(inc))
    out = os.path.join(lx, "freeze", "freeze_manifest.json")
    L.say(f"stage manifest: {len(files)} files -> {out}", args.dry_run)
    if args.dry_run:
        for f in files[:20]:
            L.say(f"  {os.path.relpath(f, lx) if L.is_inside(f, lx) else f}", True)
        return
    if os.path.exists(out):
        sys.exit(f"{out} exists")
    runlog = os.path.join(lx, "logs", "pipeline_runs.jsonl")
    head = None
    if os.path.exists(runlog):
        with open(runlog, "rb") as f:
            lines = f.read().splitlines()
        head = L.sha256_bytes(lines[-1]) if lines else None
    audit = os.path.join(lx, "logs", "audit.jsonl")
    audit_head = None
    if os.path.exists(audit):
        with open(audit, "rb") as f:
            al = f.read().splitlines()
        audit_head = L.sha256_bytes(al[-1]) if al else None
    man = {"utc": L.utcnow(), "script_sha256": L.sha256_file(os.path.abspath(__file__)),
           "files": {(os.path.relpath(f, lx) if L.is_inside(f, lx) else f): L.sha256_file(f)
                     for f in dict.fromkeys(files) if os.path.isfile(f)},
           "pipeline_runs_head_sha256": head, "audit_log_head_sha256": audit_head,
           "gates_json": args.gates, "gates_json_sha256": L.sha256_file(args.gates) if args.gates else None,
           "note": "Posted (embargoed) in OSF amendment 2 before unlock_and_merge.py is run (v1.2 sec. 2.8)."}
    L.write_json(out, man)
    print(f"freeze manifest: {len(man['files'])} files; sha256 {L.sha256_file(out)} (quote it in amendment 2)")


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--stage", choices=["usable", "manifest"], required=True)
    a.add_argument("--gates", default=None, help="gates_preunlock.json (aggregate gate results)")
    a.add_argument("--include", nargs="*", help="extra files to hash (e.g. the analysis code and DWI S3 tables)")
    a.add_argument("--dry-run", action="store_true")
    args = a.parse_args()
    (stage_usable if args.stage == "usable" else stage_manifest)(args)


if __name__ == "__main__":
    main()
