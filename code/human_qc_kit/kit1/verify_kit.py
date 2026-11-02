#!/usr/bin/env python3
"""verify_kit.py - integrity and blinding checks of the human visual-QC kit (read-only; prints a report).

Checks
  1. rate.html embedded item list == rating_sheet_blank.csv == image folders (gross and placement) == kit_key.csv.
  2. Sample: every AI grade-2 session once as a primary item; 40 AI grade-0/1 primary items; 10 duplicates of
     primary items, each >= 20 positions from its original; kit IDs unique, 'h' + 6 hex.
  3. Reproducibility: the selection, kit IDs and order are re-derived from seed 20261007 and must equal the key.
  4. Blinding: no QC ID, subject code, participant/session label or AI field appears in any rater-facing file; PNGs carry
     no text metadata; every kit PNG equals its source montage pixel for pixel outside the re-labelled title band.
  5. kit_key.csv SHA-256 equals the value recorded in README_coordinator.md (when that file exists).
Added 2026-10-07 after the adversarial review (POST HOC tooling; see ../CHANGES_2026-10-07.md):
  6. (H1) the source montages are read from the permanent archive LEADchi/qc_montages_archive (not from the temporary
     session scratchpad): the archive's own SHA256SUMS verifies (360 PNGs), and every kit source and kit image equals
     the src_sha256 / dst_sha256 recorded in build_log.json at build time.
  7. (M5) analysis/make_zip.py rebuilds the zip byte-identically from rater_kit/.
  (L1) No bytecode is written: sys.dont_write_bytecode is set before build_kit and make_zip are imported.
"""
import csv
import hashlib
import json
import os
import re
import sys

import numpy as np
from PIL import Image

KIT = "/path/to/Penn_LEAD/LEADchi/human_qc_kit"
LX = "/path/to/Penn_LEAD/LEADchi"
RK = os.path.join(KIT, "rater_kit")
sys.dont_write_bytecode = True            # read-only: importing build_kit must not create analysis/__pycache__ (L1)
sys.path.insert(0, os.path.join(KIT, "analysis"))
import build_kit as B  # noqa: E402  (constants and paths only; main() is not run)
import make_zip as MZ  # noqa: E402  (build_bytes() only; nothing is written)

ok_all = True


def check(cond, msg):
    global ok_all
    ok_all &= bool(cond)
    print(("PASS  " if cond else "FAIL  ") + msg)


def sha256(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


key = list(csv.DictReader(open(os.path.join(KIT, "_key_DO_NOT_SEND", "kit_key.csv"))))
key.sort(key=lambda r: int(r["order"]))
html = open(os.path.join(RK, "rate.html"), encoding="utf-8").read()
m = re.search(r"var ITEMS = (\[.*?\]);", html, re.S)
items = json.loads(m.group(1))
sheet = list(csv.DictReader(open(os.path.join(RK, "rating_sheet_blank.csv"))))
g_imgs = sorted(f[:-4] for f in os.listdir(os.path.join(RK, "images", "gross")) if f.endswith(".png"))
p_imgs = sorted(f[:-4] for f in os.listdir(os.path.join(RK, "images", "placement")) if f.endswith(".png"))

# 1 completeness
ids_key = [r["kit_id"] for r in key]
check([i["kit_id"] for i in items] == ids_key, f"rate.html item list ({len(items)}) == kit_key order ({len(key)})")
check([r["kit_id"] for r in sheet] == ids_key, f"rating_sheet_blank.csv ({len(sheet)}) == kit_key order")
check([i["order"] for i in items] == list(range(1, len(items) + 1)), "orders are 1..N in rate.html")
check(g_imgs == sorted(ids_key) and p_imgs == sorted(ids_key),
      f"image folders hold exactly the kit IDs (gross {len(g_imgs)}, placement {len(p_imgs)})")
check(all(int(i["forced_review"]) == int(r["forced_review"]) == int(s["forced_review"])
          for i, r, s in zip(items, key, sheet)), "forced-review flag identical in rate.html, sheet and key")
check(len(set(ids_key)) == len(ids_key) and all(re.fullmatch(r"h[0-9a-f]{6}", k) for k in ids_key),
      "kit IDs unique, format h + 6 hex")

# 2 sample
grades = list(csv.DictReader(open(B.GRADES)))
qmap = {(r["code"], r["ses_order"]): r["qc_id"] for r in csv.DictReader(open(B.QCIDS))}
ai2 = {qmap[(r["code"], r["ses_order"])] for r in grades if r["qc_grade"] == "2"}
prim = [r for r in key if not r["is_duplicate_of"]]
dups = [r for r in key if r["is_duplicate_of"]]
check({r["qc_id"] for r in prim if r["stratum"] == "AI_fail_all"} == ai2 and len(ai2) == 61,
      f"all {len(ai2)} AI grade-2 sessions are primary items")
pas = [r for r in prim if r["stratum"] == "AI_pass_sample"]
check(len(pas) == 40 and all(r["ai_grade"] in ("0", "1") for r in pas), "40 AI grade-0/1 primary items")
check(len({r["qc_id"] for r in prim}) == len(prim) == 101, "101 distinct primary sessions")
pos = {r["kit_id"]: int(r["order"]) for r in key}
qid = {r["kit_id"]: r["qc_id"] for r in key}
check(len(dups) == 10 and all(r["is_duplicate_of"] in {p["kit_id"] for p in prim} and qid[r["is_duplicate_of"]] == r["qc_id"]
                              for r in dups), "10 duplicates, each of a primary item with the same QC ID")
check(min(abs(pos[r["kit_id"]] - pos[r["is_duplicate_of"]]) for r in dups) >= B.MIN_GAP,
      f"duplicate pairs >= {B.MIN_GAP} positions apart")

# 3 reproducibility (same rng sequence as build_kit.main)
g = sorted(({**r, "ses_order": int(r["ses_order"]), "qc_grade": int(r["qc_grade"])} for r in grades),
           key=lambda r: (r["code"], r["ses_order"]))
fail = [r for r in g if r["qc_grade"] == 2]
pool = [r for r in g if r["qc_grade"] in (0, 1)]
rng = np.random.default_rng(B.SEED)
pick = sorted(rng.choice(len(pool), size=B.N_PASS_SAMPLE, replace=False).tolist())
sel = fail + [pool[i] for i in pick]
didx = sorted(rng.choice(len(sel), size=B.N_DUP, replace=False).tolist())
its = [(r, None) for r in sel] + [(sel[i], i) for i in didx]
ids, kid = set(), []
for _ in its:
    while True:
        k = "h" + format(int(rng.integers(0, 16 ** 6)), "06x")
        if k not in ids:
            ids.add(k)
            kid.append(k)
            break
while True:
    perm = rng.permutation(len(its))
    ps = {int(p): i for i, p in enumerate(perm)}
    if min(abs(ps[j] - ps[d]) for j, (_, d) in enumerate(its) if d is not None) >= B.MIN_GAP:
        break
re_ids = [kid[int(p)] for p in perm]
re_q = [qmap[(its[int(p)][0]["code"], str(its[int(p)][0]["ses_order"]))] for p in perm]
check(re_ids == ids_key and re_q == [r["qc_id"] for r in key], f"selection, kit IDs and order reproduce from seed {B.SEED}")

# 4 blinding
rater_files = [os.path.join(RK, f) for f in os.listdir(RK) if os.path.isfile(os.path.join(RK, f))]
text = "".join(open(f, encoding="utf-8", errors="ignore").read() for f in rater_files if not f.endswith(".png"))
qcids = set(qmap.values())
codes = {r["code"] for r in key} | {c for c, _ in qmap}
leak_q = [q for q in qcids if q in text]
leak_c = [c for c in codes if re.search(r"(?<![0-9a-f])" + re.escape(c) + r"(?![0-9a-f])", text)]
check(not leak_q, f"no QC ID in rater-facing text files ({len(rater_files)} files)")
check(not leak_c, "no subject code in rater-facing text files" + (f" (found {leak_c[:5]})" if leak_c else ""))
check(not re.search(r"sub-\d{4,}|ses-\d", text), "no participant/session label in rater-facing text files")
check(not re.search(r"ai_grade|qc_grade|AI_fail|AI_pass|stratum|is_duplicate", text),
      "no AI-grade or stratum field in rater-facing text files")
meta_bad, pix_bad = [], []
for r in key:
    for kind in ("gross", "placement"):
        dst = os.path.join(RK, "images", kind, r["kit_id"] + ".png")
        im = Image.open(dst)
        if {k for k in im.info if k not in ("dpi", "gamma", "srgb")} or im.text:
            meta_bad.append(dst)
        a = np.asarray(im.convert("RGB")).astype(int)
        s = np.asarray(Image.open(os.path.join(B.MONT, kind, r["qc_id"] + ".png")).convert("RGB")).astype(int)
        rows = np.where((a != s).any(axis=(1, 2)))[0]
        if a.shape != s.shape or (len(rows) and rows.max() > 35):
            pix_bad.append(dst)
check(not meta_bad, "kit PNGs carry no text metadata")
check(not pix_bad, f"kit PNGs equal their source montages outside the title band (rows 0-35); {2 * len(key)} checked")
# 4b zip == rater_kit (only), no key or analysis file inside
import zipfile  # noqa: E402
zp = os.path.join(KIT, "LEADchi_visualQC_rater_kit.zip")
if os.path.exists(zp):
    z = zipfile.ZipFile(zp)
    names = [n for n in z.namelist() if not n.endswith("/")]
    local = sorted(os.path.relpath(os.path.join(dp, f), KIT) for dp, _, fs in os.walk(RK) for f in fs if f != ".DS_Store")
    same = sorted(names) == local and all(hashlib.sha256(z.read(n)).hexdigest() == sha256(os.path.join(KIT, n)) for n in names)
    check(same, f"zip holds exactly rater_kit/ ({len(names)} files), byte-identical")
    check(not any(("key" in n.lower() and "kit_key" in n.lower()) or n.startswith(("analysis", "_key")) or
                  "coordinator" in n.lower() for n in names), "zip holds no key, analysis or coordinator file")
    rc0 = os.path.join(KIT, "README_coordinator.md")
    if os.path.exists(rc0):
        check(sha256(zp) in open(rc0, encoding="utf-8").read(), "zip SHA-256 recorded in README_coordinator.md")
    check(MZ.build_bytes() == open(zp, "rb").read(), "make_zip.py rebuilds the zip byte-identically (sorted entries, fixed "
          "timestamps, stored)")
else:
    check(False, f"{os.path.basename(zp)} exists (build it with analysis/make_zip.py)")
# 6 source archive (H1)
arch = os.path.dirname(B.MONT)
check(os.path.realpath(B.MONT).startswith(os.path.realpath(os.path.join(LX, "qc_montages_archive"))) and
      "/tmp/" not in os.path.realpath(B.MONT), f"build_kit.MONT points at the permanent archive ({B.MONT})")
sums = {}
for ln in open(os.path.join(arch, "SHA256SUMS")):
    h, f = ln.rstrip("\n").split(None, 1)
    sums[f.lstrip("*")] = h
bad_arch = [f for f, h in sums.items() if sha256(os.path.join(arch, f)) != h]
check(len(sums) == 360 and not bad_arch, f"archive SHA256SUMS verifies ({len(sums)} PNGs, {len(bad_arch)} mismatches)")
blog = json.load(open(os.path.join(KIT, "analysis", "build_log.json")))
bad_src, bad_dst, nsrc = [], [], set()
for r in key:
    for kind in ("gross", "placement"):
        ent = blog["images"][f"{kind}/{r['kit_id']}"]
        src = os.path.join(B.MONT, kind, r["qc_id"] + ".png")
        nsrc.add(src)
        if sha256(src) != ent["src_sha256"] or sums.get(f"montages/{kind}/{r['qc_id']}.png") != ent["src_sha256"]:
            bad_src.append(src)
        if sha256(os.path.join(RK, "images", kind, r["kit_id"] + ".png")) != ent["dst_sha256"]:
            bad_dst.append(r["kit_id"])
check(not bad_src, f"all {len(nsrc)} kit source montages in the archive equal build_log.json src_sha256")
check(not bad_dst, f"all {2 * len(key)} kit images equal build_log.json dst_sha256")
sha = sha256(os.path.join(KIT, "_key_DO_NOT_SEND", "kit_key.csv"))
rc = os.path.join(KIT, "README_coordinator.md")
if os.path.exists(rc):
    check(sha in open(rc, encoding="utf-8").read(), f"kit_key.csv SHA-256 {sha} recorded in README_coordinator.md")
else:
    print(f"INFO  README_coordinator.md not written yet; kit_key.csv SHA-256 {sha}")
print("ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED")
sys.exit(0 if ok_all else 1)
