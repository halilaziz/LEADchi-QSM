#!/usr/bin/env python3
"""unlock_and_merge.py - the single unlock event of LEAD-chi (task P5b, 2026-09-25).

Procedural blinding, v1.2 delta sec. 2.8 (P5 and the delta call it unlock_merge.py). Run by the USER, once, on the
Mac (standard library + OpenSSL >= 1.1.1), with a local mirror of the frozen server files. The lock private key
never leaves the Mac and its pass phrase is typed only into OpenSSL's own terminal prompt (never into an AI
session, never on a command line).

Confirmatory-stage trigger (all must hold, else nothing is decrypted):
  1. the ciphertexts (key, labels) and the label and freeze manifests match the SHA-256 values registered on OSF
     (--registered JSON typed from the registration and amendment 2)
  2. every file listed in the freeze manifest that is present in the mirror matches its frozen SHA-256, and every
     table this script merges is listed there
  3. --amendment2-utc (the OSF time stamp of amendment 2) is later than the freeze manifest and not in the future
  4. no earlier unlock record exists (a re-merge after unlock needs --rerun and is logged as such)
  5. the user types the confirmation phrase "UNLOCK LEAD-chi"
Then it decrypts the key and the locked label bundle in memory, verifies the bundle against its lock manifest,
joins labels, ages and phenotypes through the key, and writes, per available branch (engine x pipeline x roi_def),
unlocked/merged_<engine>_<pipeline>_<roi_def>.csv in the analysis schema (analysis/schema/roi_table_schema.csv; all
coded sessions, with the frozen usable / first_usable flags, QC grades and code D, the split-echo values, the
visible session flags and the key's epoch / f0 cluster), unlocked/multiverse_<engine>_<roi_def>.csv in
analysis/schema/multiverse_schema.csv (primary / other ROI mapped by the G3 branch), plus unlocked/key_plain.csv and
unlocked/UNLOCK_RECORD.json (UTC time, host, user, script and input/output hashes). The record is also appended to
<mirror>/logs/unlock_log.jsonl. --dry-run performs checks 1-4 and prints the plan without decrypting.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import getpass
import glob
import io
import json
import os
import socket
import sys
import tarfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

PHRASE = "UNLOCK LEAD-chi"
PHENO = {"prime": ("prime_total_score", "PRIME_total"),
         "bdi_1a_ms_child": ("bdi_1a_ms_child_total_score", "BDI_total"),
         "ari": ("ari_total_score", "ARI_total"),
         "short_penn_continuous_performance_test": ("cnb_cptnl_total_fp", "CPT_fp"),
         "rsas": ("rsas_total_score", "RSAS_total")}
TANNER = {"tanner_boy": [f"tanner_boy_{i}" for i in range(1, 6)], "tanner_girl": [f"tanner_girl_{i}" for i in (1, 2, 3, 5, 6)]}
DX = ["dx_psychosis", "dx_mdd", "dx_moodnos"]
SNONLY_MAP = {"chi_SNVTA": "chi_SNonly", "chi_SNVTA_ref2": "chi_SNonly_ref2", "r2s_SNVTA": "r2s_SNonly",
              "roi_vol": "roi_vol_SNonly", "fimp_SNVTA": "fimp_SNonly", "inmask_SNVTA": "inmask_SNonly",
              "chi_SNVTA_median": "chi_SNonly_median", "chi_SNVTA_impkeep": "chi_SNonly_impkeep"}
BRANCH_VARIANT = {"primary": "primary", "fallback_PDF_iLSQR": "mv_pdf_ilsqr"}
DWI_COLS = ["twICVF_pSCP_R", "b0loc_pSCP_R", "tract_len", "nSL", "FD", "NDC", "twICVF_MCP_R_matched",
            "r2s_pSCP_R", "twICVF_VTAuncross_R", "twICVF_pSCP_R_peak", "t1t2_pSCP_R", "afq_round1_list"]


def die(msg):
    sys.exit(f"UNLOCK REFUSED: {msg}")


def parse_utc(s):
    return datetime.datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=datetime.timezone.utc)


def tsv(b: bytes):
    return list(csv.DictReader(io.StringIO(b.decode("utf-8-sig")), delimiter="\t"))


def na(v):
    return "NA" if v in (None, "", "n/a", "N/A", "NaN", "nan") else v


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--mirror", required=True, help="local copy of the frozen $LX subset (blind/, locked/, "
                   "tables_coded/, flags_coded/, freeze/, qc/, logs/)")
    a.add_argument("--privkey", required=True, help="lock_priv.pem (Mac only)")
    a.add_argument("--registered", required=True, help="JSON with the SHA-256 values registered on OSF")
    a.add_argument("--amendment2", required=True, help="OSF identifier (URL/DOI) of amendment 2")
    a.add_argument("--amendment2-utc", required=True, help="OSF time stamp of amendment 2, YYYY-MM-DDTHH:MM:SSZ")
    a.add_argument("--schema", default=os.path.join(here, "..", "analysis", "schema", "roi_table_schema.csv"))
    a.add_argument("--dwi-table", default=None, help="frozen real-ID DWI S3 table (participant_id, session_id, "
                   "s3_eligible, " + ", ".join(DWI_COLS) + ")")
    a.add_argument("--out", default=None)
    a.add_argument("--confirm", default=None, help=f'the phrase "{PHRASE}" (otherwise asked interactively)')
    a.add_argument("--rerun", action="store_true", help="re-merge after the first unlock (logged as a rerun)")
    a.add_argument("--dry-run", action="store_true")
    args = a.parse_args()
    M = os.path.realpath(args.mirror)
    out = os.path.realpath(args.out or os.path.join(M, "unlocked"))
    L.assert_inside(out, M)
    now = datetime.datetime.now(datetime.timezone.utc)
    with open(args.registered, encoding="utf-8") as f:
        reg = json.load(f)
    P = {"key": os.path.join(M, "blind", "key.csv.enc"), "keyk": os.path.join(M, "blind", "key.sessionkey.enc"),
         "lab": os.path.join(M, "locked", "labels_locked.tar.enc"),
         "labk": os.path.join(M, "locked", "labels_locked.sessionkey.enc"),
         "labm": os.path.join(M, "locked", "labels_locked_manifest.json"),
         "freeze": os.path.join(M, "freeze", "freeze_manifest.json")}
    need = {"key": "key_csv_enc_sha256", "keyk": "key_sessionkey_enc_sha256", "lab": "labels_tar_enc_sha256",
            "labk": "labels_sessionkey_enc_sha256", "labm": "labels_manifest_sha256", "freeze": "freeze_manifest_sha256"}
    # ---- check 1: registered hashes
    got = {}
    for k, rk in need.items():
        if not os.path.exists(P[k]):
            die(f"{P[k]} missing from the mirror")
        got[k] = L.sha256_file(P[k])
        if rk not in reg:
            die(f"'{rk}' not in {args.registered}")
        if got[k] != reg[rk].strip().lower():
            die(f"{os.path.basename(P[k])}: SHA-256 {got[k]} != registered {reg[rk]}")
    L.say("check 1: ciphertexts and manifests match the registered SHA-256 values", args.dry_run)
    # ---- check 2: freeze manifest
    with open(P["freeze"], encoding="utf-8") as f:
        fz = json.load(f)
    bad = []
    for rel, h in fz["files"].items():
        p = rel if os.path.isabs(rel) else os.path.join(M, rel)
        if os.path.exists(p) and L.sha256_file(p) != h:
            bad.append(rel)
    if bad:
        die(f"{len(bad)} frozen files changed after the freeze: {bad[:5]}")
    tables = {}
    for eng in ("qsmxt", "sepia"):
        for pipe, var in BRANCH_VARIANT.items():
            t = os.path.join("tables_coded", f"roi_{eng}_{var}.csv")
            for rd in ("SNVTA", "SN_only"):
                u = os.path.join("freeze", f"usable_{eng}_{pipe}_{rd}.csv")
                if os.path.exists(os.path.join(M, t)) and os.path.exists(os.path.join(M, u)):
                    for rel in (t, u):
                        if rel not in fz["files"]:
                            die(f"{rel} is not in the freeze manifest")
                    tables[(eng, pipe, rd)] = (os.path.join(M, t), os.path.join(M, u))
    if not tables:
        die("no (ROI table, usable list) pair found in the mirror")
    extra_tables = {}
    for eng in ("qsmxt", "sepia"):
        for rel in (os.path.join("tables_coded", f"roi_{eng}_split13.csv"), os.path.join("tables_coded", f"roi_{eng}_split24.csv"),
                    os.path.join("tables_coded", f"multiverse_{eng}.csv")):
            if os.path.exists(os.path.join(M, rel)):
                if rel not in fz["files"]:
                    die(f"{rel} is not in the freeze manifest")
                extra_tables[rel] = os.path.join(M, rel)
    for rel in (os.path.join("qc", "grades_final.csv"), os.path.join("blind", "session_flags.csv")):
        if not os.path.exists(os.path.join(M, rel)) or rel not in fz["files"]:
            die(f"{rel} missing from the mirror or the freeze manifest")
    if args.dwi_table:
        dh = L.sha256_file(args.dwi_table)
        if dh not in fz["files"].values():
            die("the DWI table is not in the freeze manifest (hash not found)")
    L.say(f"check 2: freeze manifest verified; branches to merge: {sorted(tables)}", args.dry_run)
    # ---- check 3: order of events
    t_freeze, t_am2 = parse_utc(fz["utc"]), parse_utc(args.amendment2_utc)
    if not (t_freeze <= t_am2 <= now):
        die(f"amendment 2 time {args.amendment2_utc} must lie between the freeze ({fz['utc']}) and now")
    L.say(f"check 3: freeze {fz['utc']} <= amendment 2 {args.amendment2_utc} <= now", args.dry_run)
    # ---- check 4: single unlock
    rec_p = os.path.join(out, "UNLOCK_RECORD.json")
    if os.path.exists(rec_p) and not args.rerun:
        die(f"{rec_p} exists: the unlock already happened (use --rerun to re-merge; it is logged)")
    L.say("check 4: no earlier unlock record" if not os.path.exists(rec_p) else "check 4: RERUN after unlock",
          args.dry_run)
    if args.dry_run:
        L.say(f"would decrypt {P['key']} and {P['lab']} (OpenSSL prompts for the pass phrase) and write "
              f"{len(tables)} merged tables to {out}", True)
        return
    # ---- check 5: confirmation
    phrase = args.confirm if args.confirm is not None else input(f'Type "{PHRASE}" to unlock: ')
    if phrase.strip() != PHRASE:
        die("confirmation phrase not given")
    L.openssl_bin()
    t_unlock = L.utcnow()
    print("OpenSSL will now ask twice for the lock-key pass phrase (key file, then label bundle).")
    key_b = L.hybrid_decrypt(P["key"], P["keyk"], args.privkey)
    lab_b = L.hybrid_decrypt(P["lab"], P["labk"], args.privkey)
    with open(P["labm"], encoding="utf-8") as f:
        labm = json.load(f)
    files = {}
    with tarfile.open(fileobj=io.BytesIO(lab_b), mode="r") as tar:
        for m in tar.getmembers():
            if m.isfile():
                files[m.name] = tar.extractfile(m).read()
    for rel, meta in labm["files"].items():
        if rel not in files or L.sha256_bytes(files[rel]) != meta["sha256"]:
            die(f"label bundle member {rel} does not match its lock manifest")
    key = {(r["code"], int(r["ses_order"])): r for r in csv.DictReader(io.StringIO(key_b.decode()))}
    part = {r["participant_id"]: r for r in tsv(files["participants.tsv"])}
    ages, age1 = {}, {}
    for rel, b in files.items():
        if rel.endswith("_sessions.tsv"):
            pid = os.path.basename(rel).split("_sessions.tsv")[0]
            for r in tsv(b):
                try:
                    ages[(pid, r["session_id"])] = float(r["age"])
                except (KeyError, ValueError):
                    pass
                if r.get("session_id") == "ses-1":
                    age1[pid] = r.get("age")
    pheno = {}
    for name, (col, new) in PHENO.items():
        for r in tsv(files[f"phenotype/{name}.tsv"]):
            pheno[(r["participant_id"], r["session_id"], new)] = na(r.get(col))
    for r in tsv(files["phenotype/prime.tsv"]):     # PRIME items 1-12 without the duration items (S2 sensitivity)
        try:
            pheno[(r["participant_id"], r["session_id"], "PRIME_items12")] = sum(float(r[f"prime_{i}"]) for i in range(1, 13))
        except (KeyError, ValueError, TypeError):
            pheno[(r["participant_id"], r["session_id"], "PRIME_items12")] = "NA"
    for name, cols in TANNER.items():
        for r in tsv(files[f"phenotype/{name}.tsv"]):
            for c in cols:
                pheno[(r["participant_id"], r["session_id"], c)] = na(r.get(c))
    for r in tsv(files["phenotype/substance.tsv"]):  # prescribed-stimulant endorser at the same session
        pheno[(r["participant_id"], r["session_id"], "stim_prescribed")] = int(
            na(r.get("substance_othr_010")) in ("1", "2") and na(r.get("substance_othr_011")) == "1")
    grades = {(r["code"], int(r["ses_order"])): r for r in L.read_csv(os.path.join(M, "qc", "grades_final.csv"))}
    sflags = {(r["code"], int(r["ses_order"])): r for r in L.read_csv(os.path.join(M, "blind", "session_flags.csv"))}
    split = {}
    for eng in ("qsmxt", "sepia"):
        for sp in ("split13", "split24"):
            rel = os.path.join("tables_coded", f"roi_{eng}_{sp}.csv")
            if rel in extra_tables:
                split[(eng, sp)] = {(r["code"], int(r["ses_order"])): r for r in L.read_csv(extra_tables[rel])}
    dwi = {}
    if args.dwi_table:
        with open(args.dwi_table, newline="", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                dwi[(r["participant_id"], r["session_id"])] = r
    schema_cols = []
    if os.path.exists(args.schema):
        with open(args.schema, newline="", encoding="utf-8") as f:
            schema_cols = [r["column"] for r in csv.DictReader(f)]
    os.makedirs(out, exist_ok=True)
    outputs, warnings = {}, []
    for (eng, pipe, rd), (tp, up) in sorted(tables.items()):
        roi = {(r["code"], int(r["ses_order"])): r for r in L.read_csv(tp)}
        use = {(r["code"], int(r["ses_order"])): r for r in L.read_csv(up)}
        rows = []
        for k, r in sorted(roi.items()):
            if k not in key:
                die(f"coded session {k} not in the key")
            kk = key[k]
            pid, sid = kk["participant_id"], kk["session_id"]
            p = part.get(pid)
            if p is None:
                die(f"{pid} not in participants.tsv")
            row = dict(r)
            if rd == "SN_only":
                for dst, src in SNONLY_MAP.items():
                    row[f"{dst}_SNVTAcomplex"] = r.get(dst)
                    row[dst] = r.get(src)
            roi_col = "chi_SNonly" if rd == "SN_only" else "chi_SNVTA"
            for sp in ("split13", "split24"):
                row[f"chi_SNVTA_{sp}"] = split.get((eng, sp), {}).get(k, {}).get(roi_col, "NA")
            row["chi_SNVTA_harm"] = r.get(f"{roi_col}_harm", "NA")
            g_ = grades.get(k, {})
            row["qc_grade"] = g_.get("qc_grade", "NA"); row["dental_D"] = g_.get("dental_D", "NA")
            sf = sflags.get(k, {})
            row["t1_artifact"] = sf.get("t1_artifact", "NA"); row["sw_class"] = sf.get("sw_class", "NA")
            row["epoch"] = kk.get("epoch", "NA"); row["f0_cluster"] = kk.get("f0_cluster", "NA")
            u = use.get(k, {})
            row.update({"data_origin": "real", "engine": eng, "pipeline": pipe, "roi_def": rd,
                        "participant_id": pid, "session_id": sid, "ses_order": k[1], "code": k[0],
                        "study_group": p["study_group"], "dx_adhd": p["dx_adhd"], "sex": p["sex"],
                        "race": p["race"], "ethnicity": na(p.get("ethnicity")),
                        "age": ages.get((pid, sid), "NA"), "age_ses1": na(age1.get(pid)),
                        "usable": u.get("usable", "0"), "fail_reason": u.get("fail_reason", "not_in_usable_list"),
                        "first_usable": u.get("first_usable", "0")})
            for new in [v[1] for v in PHENO.values()] + ["PRIME_items12", "stim_prescribed"] + TANNER["tanner_boy"] + TANNER["tanner_girl"]:
                row[new] = pheno.get((pid, sid, new), "NA")
            for c in DX:
                row[c] = na(p.get(c))
            d = dwi.get((pid, sid), {})
            row["s3_eligible"] = d.get("s3_eligible", "0") if d else "0"
            for c in DWI_COLS:
                row[c] = na(d.get(c))
            row["dwi_FD"] = na(d.get("FD"))
            row["dx_other"] = ";".join(f"{c}={p[c]}" for c in p if c.startswith("dx_") and c != "dx_adhd")
            rows.append(row)
        # chronological check (reported, never corrected)
        by = {}
        for r in rows:
            by.setdefault(r["participant_id"], []).append((int(r["ses_order"]), r["age"]))
        for pid, lst in by.items():
            ag = [x[1] for x in sorted(lst)]
            if any(isinstance(x, str) for x in ag) or any(b <= a_ for a_, b in zip(ag, ag[1:])):
                warnings.append(f"{pid}: ses_order not strictly increasing in age or age missing")
        extra = sorted({c for r in rows for c in r} - set(schema_cols))
        cols = [c for c in schema_cols if any(c in r for r in rows)] + extra
        missing = [c for c in schema_cols if c not in cols]
        if missing:
            warnings.append(f"{eng}/{pipe}/{rd}: schema columns absent: {missing}")
        op = os.path.join(out, f"merged_{eng}_{pipe}_{rd}.csv")
        with open(op, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            for r in rows:
                w.writerow({c: na(r.get(c)) if not isinstance(r.get(c), (int, float)) else r.get(c) for c in cols})
        outputs[op] = L.sha256_file(op)
    for eng in ("qsmxt", "sepia"):
        rel = os.path.join("tables_coded", f"multiverse_{eng}.csv")
        if rel not in extra_tables:
            continue
        mrows = L.read_csv(extra_tables[rel])
        for rd, (prim, oth) in (("SNVTA", ("SNVTA", "SNonly")), ("SN_only", ("SNonly", "SNVTA"))):
            op = os.path.join(out, f"multiverse_{eng}_{rd}.csv")
            with open(op, "w", newline="", encoding="utf-8") as f:
                cols = ["participant_id", "session_id", "bf", "inv", "chi_wb_mean", "chi_wb_median", "chi_other_wb_mean",
                        "chi_csf_mean", "chi_csf_median", "chi_other_csf_mean"]
                w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
                w.writeheader()
                for r in mrows:
                    k = (r["code"], int(r["ses_order"]))
                    if k not in key:
                        die(f"multiverse: coded session {k} not in the key")
                    w.writerow({"participant_id": key[k]["participant_id"], "session_id": key[k]["session_id"],
                                "bf": r["bf"], "inv": r["inv"],
                                "chi_wb_mean": na(r[f"chi_{prim}_wb_mean"]), "chi_wb_median": na(r[f"chi_{prim}_wb_median"]),
                                "chi_other_wb_mean": na(r[f"chi_{oth}_wb_mean"]), "chi_csf_mean": na(r[f"chi_{prim}_csf_mean"]),
                                "chi_csf_median": na(r[f"chi_{prim}_csf_median"]), "chi_other_csf_mean": na(r[f"chi_{oth}_csf_mean"])})
            outputs[op] = L.sha256_file(op)
    kp = os.path.join(out, "key_plain.csv")
    with open(kp, "wb") as f:
        f.write(key_b)
    outputs[kp] = L.sha256_file(kp)
    rec = {"event": "rerun" if os.path.exists(rec_p) else "unlock", "utc_unlock": t_unlock,
           "utc_end": L.utcnow(), "host": socket.gethostname(), "user": getpass.getuser(),
           "script_sha256": L.sha256_file(os.path.abspath(__file__)),
           "lx_common_sha256": L.sha256_file(os.path.abspath(L.__file__)),
           "amendment2": args.amendment2, "amendment2_utc": args.amendment2_utc, "freeze_utc": fz["utc"],
           "inputs_sha256": got, "dwi_table_sha256": L.sha256_file(args.dwi_table) if args.dwi_table else None,
           "branches": [list(b) for b in sorted(tables)], "outputs_sha256": outputs, "warnings": warnings,
           "confirmation": "interactive" if args.confirm is None else "argument"}
    recs = []
    if os.path.exists(rec_p):
        with open(rec_p, encoding="utf-8") as f:
            recs = json.load(f)
            recs = recs if isinstance(recs, list) else [recs]
    with open(rec_p, "w", encoding="utf-8") as f:
        json.dump(recs + [rec], f, indent=2)
    lg = os.path.join(M, "logs", "unlock_log.jsonl")
    os.makedirs(os.path.dirname(lg), exist_ok=True)
    with open(lg, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec, sort_keys=True) + "\n")
    print(f"UNLOCKED at {t_unlock}. Wrote {len(outputs)} files to {out}; record {rec_p} "
          f"(sha256 {L.sha256_file(rec_p)}). Warnings: {len(warnings)}")
    for w_ in warnings:
        print("  WARNING:", w_)


if __name__ == "__main__":
    main()
