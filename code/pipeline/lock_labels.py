#!/usr/bin/env python3
"""lock_labels.py - lock the group labels, ages and phenotypes in one encrypted file (task P5b, 2026-09-25).

User decision U1 (2026-09-25): group labels are locked in a separate file and unlocked only at the confirmatory
stage. RUN ONCE BY THE USER (standard library + openssl >= 1.1.1; server /bin/python3 or the Mac), after the lock
key pair exists and before the OSF registration is posted.

Inputs (read-only; the only pipeline script besides unlock_and_merge.py allowed to open them):
  <bids>/participants.tsv, participants.json          study_group, dx_*, sex, race, ethnicity, ages
  <bids>/sub-*/sub-*_sessions.tsv (+ .json)           age at each session, acq_time, in-scanner task scores
  <bids>/phenotype/<name>.tsv (+ .json)               --phenotypes (default: the instruments the protocol uses)

Outputs ($LX/locked/ by default; nothing else is written):
  labels_locked.tar.enc          AES-256-CBC/PBKDF2 ciphertext of an in-memory tar of the files above
  labels_locked.sessionkey.enc   RSA-OAEP-wrapped session secret (lock_pub.pem)
  labels_locked_manifest.json    per-file SHA-256 and size of the plaintext files, counts, ciphertext SHA-256,
                                 lock-key fingerprint, UTC time, script SHA-256 (no content)

The plaintext never touches the disk. The manifest SHA-256 values let unlock_and_merge.py verify, after
decryption, that the unlocked files are byte-identical to what was locked, and let anyone check the public
dataset against them. --dry-run lists the files and planned outputs only.
"""
from __future__ import annotations

import argparse
import glob
import io
import os
import sys
import tarfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

DEFAULT_PHENO = ["prime", "bdi_1a_ms_child", "ari", "tanner_boy", "tanner_girl", "substance",
                 "short_penn_continuous_performance_test", "rsas"]


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--bids", default=os.path.join(L.pl_root(), "ds007116"))
    a.add_argument("--lock-pubkey", default=os.path.join(L.lx_root(), "prereg", "lock_pub.pem"))
    a.add_argument("--out", default=os.path.join(L.lx_root(), "locked"))
    a.add_argument("--root", default=None, help="write guard root (default $LX; use the Mac mirror when run locally)")
    a.add_argument("--phenotypes", default=",".join(DEFAULT_PHENO))
    a.add_argument("--dry-run", action="store_true")
    args = a.parse_args()
    dry = args.dry_run
    root = args.root or L.lx_root()
    out = os.path.realpath(args.out)
    L.assert_inside(out, root)
    names = ["labels_locked.tar.enc", "labels_locked.sessionkey.enc", "labels_locked_manifest.json"]
    if any(os.path.exists(os.path.join(out, n)) for n in names):
        sys.exit(f"REFUSED: locked label files already exist in {out} (lock_labels.py runs once)")

    files = []
    for fn in ("participants.tsv", "participants.json"):
        p = os.path.join(args.bids, fn)
        if os.path.exists(p):
            files.append(p)
    if not files or not files[0].endswith("participants.tsv"):
        sys.exit(f"{args.bids}/participants.tsv not found")
    files += sorted(glob.glob(os.path.join(args.bids, "sub-*", "sub-*_sessions.tsv")))
    files += sorted(glob.glob(os.path.join(args.bids, "sub-*", "sub-*_sessions.json")))
    missing = []
    for ph in [x.strip() for x in args.phenotypes.split(",") if x.strip()]:
        for ext in (".tsv", ".json"):
            p = os.path.join(args.bids, "phenotype", ph + ext)
            if os.path.exists(p):
                files.append(p)
            elif ext == ".tsv":
                missing.append(ph)
    if missing:
        sys.exit("phenotype files not found: " + ", ".join(missing))

    n_sess = sum(1 for f in files if f.endswith("_sessions.tsv"))
    L.say(f"files to lock: {len(files)} (participants 1, sessions.tsv {n_sess}, phenotype "
          f"{sum(1 for f in files if '/phenotype/' in f)})", dry)
    manifest = {"utc": L.utcnow(), "script_sha256": L.sha256_file(os.path.abspath(__file__)),
                "lx_common_sha256": L.sha256_file(os.path.abspath(L.__file__)),
                "bids_root": os.path.realpath(args.bids), "files": {}}
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for p in files:
            rel = os.path.relpath(p, args.bids)
            with open(p, "rb") as f:
                b = f.read()
            manifest["files"][rel] = {"sha256": L.sha256_bytes(b), "bytes": len(b)}
            ti = tarfile.TarInfo(rel)
            ti.size, ti.mtime, ti.mode, ti.uname, ti.gname = len(b), 0, 0o600, "", ""
            tar.addfile(ti, io.BytesIO(b))
    with open(os.path.join(args.bids, "participants.tsv"), encoding="utf-8") as f:
        manifest["n_participant_rows"] = sum(1 for _ in f) - 1
    manifest["n_sessions_files"] = n_sess
    if dry:
        L.say(f"would lock {len(files)} files ({buf.tell()} bytes in memory) into {out}/{names[0]}", True)
        return
    if not os.path.exists(args.lock_pubkey):
        sys.exit(f"{args.lock_pubkey} missing (README step 1)")
    bits = L.rsa_bits(args.lock_pubkey)
    if bits < 3072:
        sys.exit(f"lock public key has {bits} bits; >= 3072 required")
    enc = L.hybrid_encrypt(buf.getvalue(), args.lock_pubkey, os.path.join(out, names[0]),
                           os.path.join(out, names[1]))
    del buf
    manifest.update({"labels_tar_enc_sha256": enc["data_sha256"], "labels_sessionkey_enc_sha256": enc["key_sha256"],
                     "lock_pub_sha256": L.sha256_file(args.lock_pubkey), "lock_pub_rsa_bits": bits})
    L.write_json(os.path.join(out, names[2]), manifest)
    rr = L.RunRecord(__file__, False, log_dir=os.path.join(root, "logs"))
    for n in names:
        rr.add_output(os.path.join(out, n))
    rr.finish()
    print(f"locked {len(files)} files\n  {names[0]}  sha256 {enc['data_sha256']}\n"
          f"  {names[1]}  sha256 {enc['key_sha256']}\n  manifest     sha256 "
          f"{L.sha256_file(os.path.join(out, names[2]))}\nRecord these three hashes in the OSF registration.")


if __name__ == "__main__":
    main()
