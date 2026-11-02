#!/usr/bin/env python3
"""make_zip.py - build LEADchi_visualQC_rater_kit.zip reproducibly from rater_kit/ (POST HOC tooling, 2026-10-07).

Added after the adversarial review (item M5). The zip holds rater_kit/ ONLY (never the key, analysis/ or the
coordinator README). It is byte-reproducible: entries are sorted by path, every entry carries the same fixed timestamp
(2026-10-07 00:00:00), fixed Unix permissions (dirs 0755, files 0644), no extra fields and no archive comment, and all
entries are STORED (no compression; the PNGs are already compressed), so the bytes do not depend on the zlib version.
.DS_Store and other dot-files are left out.

Usage
  python3 make_zip.py                 # writes ../LEADchi_visualQC_rater_kit.zip (refuses to overwrite)
  python3 make_zip.py --out <path>    # another target (refuses to overwrite)
  python3 make_zip.py --check         # rebuild in memory and compare with the existing zip; exit 1 if different
"""
import argparse
import hashlib
import io
import os
import sys
import zipfile

KIT = "/path/to/Penn_LEAD/LEADchi/human_qc_kit"
SRC = os.path.join(KIT, "rater_kit")
ZIP = os.path.join(KIT, "LEADchi_visualQC_rater_kit.zip")
STAMP = (2026, 10, 7, 0, 0, 0)


def entries():
    """Sorted list of (archive name, absolute path or None for a directory)."""
    out = [("rater_kit/", None)]
    for dp, dns, fns in os.walk(SRC):
        dns[:] = sorted(d for d in dns if not d.startswith("."))
        rel = os.path.relpath(dp, KIT).replace(os.sep, "/")
        for d in dns:
            out.append((f"{rel}/{d}/", None))
        for f in sorted(fns):
            if not f.startswith("."):
                out.append((f"{rel}/{f}", os.path.join(dp, f)))
    return sorted(out, key=lambda e: e[0])


def build_bytes():
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_STORED) as z:
        for name, path in entries():
            zi = zipfile.ZipInfo(name, date_time=STAMP)
            zi.create_system = 3                                  # Unix
            zi.compress_type = zipfile.ZIP_STORED
            if path is None:
                zi.external_attr = (0o40755 << 16) | 0x10         # directory
                z.writestr(zi, b"")
            else:
                zi.external_attr = 0o100644 << 16
                with open(path, "rb") as f:
                    z.writestr(zi, f.read())
    return buf.getvalue()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=ZIP)
    ap.add_argument("--check", action="store_true", help="compare a fresh in-memory build with the existing zip")
    a = ap.parse_args()
    data = build_bytes()
    sha = hashlib.sha256(data).hexdigest()
    n = sum(1 for nm, p in entries() if p is not None)
    if a.check:
        if not os.path.exists(a.out):
            sys.exit(f"{a.out} does not exist")
        old = open(a.out, "rb").read()
        same = old == data
        print(f"{'IDENTICAL' if same else 'DIFFERENT'}: rebuild {sha} vs existing {hashlib.sha256(old).hexdigest()}")
        sys.exit(0 if same else 1)
    if os.path.exists(a.out):
        sys.exit(f"{a.out} exists: move it away first (nothing is overwritten)")
    with open(a.out, "wb") as f:
        f.write(data)
    print(f"wrote {a.out}: {n} files, {len(data)} bytes, sha256 {sha}")


if __name__ == "__main__":
    main()
