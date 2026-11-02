#!/usr/bin/env python3
"""verify_kit2.py - integrity, sample and blinding checks of human visual-QC kit 2 (read-only; prints a report).

Checks (kit-1 verify_kit.py adapted to kit 2)
  1. rate.html item list == rating_sheet_blank.csv == image folders == kit2_key.csv order; forced-review flag identical
     in all three and equal to mirror/qc/rating_pass1_blank_template.csv; kit-2 IDs unique, 'h' + 6 hex, none used in kit 1.
  2. Sample: 79 new items = every AI grade-0/1 session not among kit 1's primary items; 21 repeats (11 of kit-1 AI-fail
     primaries, 10 of kit-1 AI-pass primaries), each with the same QC ID as its kit-1 item; kit-2 new items and kit-1
     primaries are disjoint and together are exactly the 180 AI-graded sessions; repeats spread (3-6 per block of 20).
  3. Reproducibility: selection, kit IDs and order re-derived from seed 20261028 (build_kit2.select) equal the key.
  4. Blinding: no QC ID, subject code, participant/session label, AI field, item-type or composition text in any
     rater-facing file; PNGs carry no text metadata; every PNG equals its archive montage pixel for pixel outside the
     title band (rows 0-35), and every changed row lies inside the band recorded in build_log_kit2.json.
  5. rate.html equals kit 1's rate.html except for the embedded item list and its storage tag (same interface).
  6. Archive SHA256SUMS verifies (360 PNGs); src/dst SHA-256 equal build_log_kit2.json; key SHA-256 equals the log.
  7. The zip holds exactly rater_kit2/ (byte-identical), nothing else, and make_zip_kit2.py rebuilds it byte-identically.
No bytecode is written (sys.dont_write_bytecode).
"""
import csv
import hashlib
import json
import os
import re
import sys
import zipfile

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

K1 = "/path/to/Penn_LEAD/LEADchi/human_qc_kit"
K2 = os.path.join(K1, "kit2_2026-10-28")
RK = os.path.join(K2, "rater_kit2")
sys.path.insert(0, os.path.join(K2, "analysis"))
import build_kit2 as B2  # noqa: E402
import make_zip_kit2 as MZ  # noqa: E402
B1 = B2.B1

ok_all = True


def check(cond, msg):
    global ok_all
    ok_all &= bool(cond)
    print(("PASS  " if cond else "FAIL  ") + msg)


def sha256(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


rd = lambda p: list(csv.DictReader(open(p, newline="")))
key = sorted(rd(B2.KEY2), key=lambda r: int(r["order"]))
key1 = rd(B2.KEY1)
html = open(os.path.join(RK, "rate.html"), encoding="utf-8").read()
items = json.loads(re.search(r"var ITEMS = (\[.*?\]);", html, re.S).group(1))
sheet = rd(os.path.join(RK, "rating_sheet_blank.csv"))
imgs = {k: sorted(f[:-4] for f in os.listdir(os.path.join(RK, "images", k)) if f.endswith(".png")) for k in ("gross", "placement")}
frev = {r["qc_id"]: r["forced_review"] for r in rd(B1.BLANK1)}

# 1 completeness
ids = [r["kit2_id"] for r in key]
check([i["kit_id"] for i in items] == ids, f"rate.html item list ({len(items)}) == kit2_key order ({len(key)})")
check([r["kit_id"] for r in sheet] == ids, f"rating_sheet_blank.csv ({len(sheet)}) == kit2_key order")
check([i["order"] for i in items] == list(range(1, len(items) + 1)) and [int(r["order"]) for r in key] == list(range(1, 101)),
      "orders are 1..100")
check(imgs["gross"] == sorted(ids) and imgs["placement"] == sorted(ids),
      f"image folders hold exactly the kit-2 IDs (gross {len(imgs['gross'])}, placement {len(imgs['placement'])})")
check(all(str(i["forced_review"]) == r["forced_review"] == s["forced_review"] == frev[r["qc_id"]]
          for i, r, s in zip(items, key, sheet)),
      "forced-review flag identical in rate.html, sheet, key and rating_pass1_blank_template.csv (kit-1 source)")
k1ids = {r["kit_id"] for r in key1}
check(len(set(ids)) == 100 and all(re.fullmatch(r"h[0-9a-f]{6}", k) for k in ids) and not (set(ids) & k1ids),
      "kit-2 IDs unique, format h + 6 hex, none used in kit 1")

# 2 sample
grades = rd(B1.GRADES)
qmap = {(r["code"], r["ses_order"]): r["qc_id"] for r in rd(B1.QCIDS)}
gq = {qmap[(r["code"], r["ses_order"])]: r for r in grades}
prim1 = {r["qc_id"]: r for r in key1 if not r["is_duplicate_of"]}
new = [r for r in key if r["stratum"] == "AI_pass_new"]
rep = [r for r in key if r["is_kit1_repeat"] == "1"]
ai_pass = {q for q, r in gq.items() if r["qc_grade"] in ("0", "1")}
check(len(new) == 79 and {r["qc_id"] for r in new} == ai_pass - set(prim1) and all(r["ai_grade"] in ("0", "1") for r in new),
      f"79 new items = the {len(ai_pass)} AI grade-0/1 sessions minus the {len(ai_pass & set(prim1))} in kit 1")
check(len(rep) == 21 and len({r["qc_id"] for r in rep}) == 21 and
      all(r["kit1_id"] in {p["kit_id"] for p in prim1.values()} and prim1[r["qc_id"]]["kit_id"] == r["kit1_id"] for r in rep),
      "21 repeats, each of a distinct kit-1 primary item with the same QC ID")
nf = sum(prim1[r["qc_id"]]["stratum"] == "AI_fail_all" and r["stratum"] == "repeat_of_kit1_AI_fail" for r in rep)
npas = sum(prim1[r["qc_id"]]["stratum"] == "AI_pass_sample" and r["stratum"] == "repeat_of_kit1_AI_pass" for r in rep)
check(nf == 11 and npas == 10, f"repeat strata: {nf} kit-1 AI-fail + {npas} kit-1 AI-pass")
union = {r["qc_id"] for r in new} | set(prim1)
check(len(gq) == 180 and not ({r["qc_id"] for r in new} & set(prim1)) and union == set(gq),
      f"kit-2 new ({len(new)}) + kit-1 primaries ({len(prim1)}) = all {len(gq)} AI-graded sessions, disjoint")
check(all(r["ai_grade"] == gq[r["qc_id"]]["qc_grade"] and r["ai_placement"] == gq[r["qc_id"]]["placement"] for r in key),
      "key AI grade/placement equal grades_final.csv")
isrep = [r["is_kit1_repeat"] == "1" for r in key]
blocks = [sum(isrep[b:b + 20]) for b in range(0, 100, 20)]
check(all(3 <= c <= 6 for c in blocks), f"repeats spread: per block of 20 = {blocks}")

# 3 reproducibility
ordered, info = B2.select(rd(B1.GRADES), rd(B2.KEY1), {(r["code"], int(r["ses_order"])): r["qc_id"] for r in rd(B1.QCIDS)})
check([o["kit2_id"] for o in ordered] == ids and [o["qc_id"] for o in ordered] == [r["qc_id"] for r in key] and
      [o["kit1_id"] for o in ordered] == [r["kit1_id"] for r in key], f"selection, kit IDs and order reproduce from seed {B2.SEED}")

# 4 blinding (text)
rfiles = [os.path.join(RK, f) for f in os.listdir(RK) if os.path.isfile(os.path.join(RK, f))]
text = "".join(open(f, encoding="utf-8", errors="ignore").read() for f in rfiles)
codes = {c for c, _ in qmap}
check(not [q for q in set(qmap.values()) if q in text], f"no QC ID in rater-facing text files ({len(rfiles)} files)")
lc = [c for c in codes if re.search(r"(?<![0-9a-f])" + re.escape(c) + r"(?![0-9a-f])", text)]
check(not lc, "no subject code in rater-facing text files" + (f" (found {lc[:5]})" if lc else ""))
check(not re.search(r"sub-\d{4,}|ses-\d", text), "no participant/session label in rater-facing text files")
check(not re.search(r"ai_grade|qc_grade|AI_fail|AI_pass|AI-fail|AI-pass|stratum|is_duplicate|is_kit1|kit1_id|repeat_of|"
                    r"kit2_key|kit_key", text, re.I), "no AI-grade, item-type or key field in rater-facing text files")
check(not [k for k in k1ids if k in text], "no kit-1 ID in rater-facing text files")
comp = re.findall(r"[^.\n]*\b(?:79|21|11|89)\b[^.\n]*", open(os.path.join(RK, "README_rater.md"), encoding="utf-8").read())
check(not [c for c in comp if re.search(r"repeat|tekrar|new|yeni|pass|fail|geç|kal", c, re.I)],
      "README_rater.md states no composition count (79 / 21 / 11 / 89)" + (f": {comp}" if comp else ""))

# 4b images
blog = json.load(open(B2.LOG))
meta_bad, pix_bad, band_bad, maxrow = [], [], [], 0
for r in key:
    for kind in ("gross", "placement"):
        dst = os.path.join(RK, "images", kind, r["kit2_id"] + ".png")
        im = Image.open(dst)
        if {k for k in im.info if k not in ("dpi", "gamma", "srgb")} or im.text:
            meta_bad.append(dst)
        a = np.asarray(im.convert("RGB")).astype(int)
        s = np.asarray(Image.open(os.path.join(B1.MONT, kind, r["qc_id"] + ".png")).convert("RGB")).astype(int)
        if a.shape != s.shape:
            pix_bad.append(dst)
            continue
        rows = np.where((a != s).any(axis=(1, 2)))[0]
        if len(rows):
            maxrow = max(maxrow, int(rows.max()))
            if rows.max() > 35:
                pix_bad.append(dst)
            b0, b1 = blog["images"][f"{kind}/{r['kit2_id']}"]["blanked_rows"]
            if rows.min() < b0 or rows.max() > b1:
                band_bad.append(dst)
check(not meta_bad, "kit-2 PNGs carry no text metadata")
check(not pix_bad, f"kit-2 PNGs equal their archive montages outside the title band (rows 0-35); {2 * len(key)} checked; "
                   f"last changed row {maxrow}")
check(not band_bad, "every changed row lies inside the blanked title band recorded in build_log_kit2.json")

# 5 interface identical to kit 1
h1 = open(os.path.join(K1, "rater_kit", "rate.html"), encoding="utf-8").read()
strip = lambda h: re.sub(r'var KIT_TAG = "[0-9a-f]{12}";', "", re.sub(r"var ITEMS = \[.*?\];", "", h, flags=re.S))
check(strip(h1) == strip(html), "rate.html equals kit-1 rate.html apart from the item list and storage tag")
tag2 = re.search(r'var KIT_TAG = "([0-9a-f]{12})"', html).group(1)
tag1 = re.search(r'var KIT_TAG = "([0-9a-f]{12})"', h1).group(1)
check(tag1 != tag2 and not re.search(r"https?://|//cdn|fetch\(|XMLHttpRequest", html),
      f"storage tag differs from kit 1 ({tag2} vs {tag1}); no network reference in rate.html")

# 6 archive and log
arch = os.path.dirname(B1.MONT)
sums = {}
for ln in open(os.path.join(arch, "SHA256SUMS")):
    h, f = ln.rstrip("\n").split(None, 1)
    sums[f.lstrip("*")] = h
bad_arch = [f for f, h in sums.items() if sha256(os.path.join(arch, f)) != h]
check(len(sums) == 360 and not bad_arch, f"archive SHA256SUMS verifies ({len(sums)} PNGs, {len(bad_arch)} mismatches)")
bad_src, bad_dst = [], []
for r in key:
    for kind in ("gross", "placement"):
        ent = blog["images"][f"{kind}/{r['kit2_id']}"]
        src = os.path.join(B1.MONT, kind, r["qc_id"] + ".png")
        if sha256(src) != ent["src_sha256"] or sums.get(f"montages/{kind}/{r['qc_id']}.png") != ent["src_sha256"]:
            bad_src.append(src)
        if sha256(os.path.join(RK, "images", kind, r["kit2_id"] + ".png")) != ent["dst_sha256"]:
            bad_dst.append(r["kit2_id"])
check(not bad_src, "all 200 kit-2 source montages equal build_log src_sha256 and the archive SHA256SUMS")
check(not bad_dst, "all 200 kit-2 images equal build_log dst_sha256")
check(sha256(B2.KEY2) == blog["kit2_key_sha256"], f"kit2_key.csv SHA-256 {blog['kit2_key_sha256']} equals build log")

# 7 zip
zp = MZ.ZIP
z = zipfile.ZipFile(zp)
names = [n for n in z.namelist() if not n.endswith("/")]
local = sorted(os.path.relpath(os.path.join(dp, f), K2) for dp, _, fs in os.walk(RK) for f in fs if not f.startswith("."))
check(sorted(names) == local and all(hashlib.sha256(z.read(n)).hexdigest() == sha256(os.path.join(K2, n)) for n in names),
      f"zip holds exactly rater_kit2/ ({len(names)} files), byte-identical")
check(all(n.startswith("rater_kit2/") for n in z.namelist()) and
      not any(re.search(r"key|analysis|coordinator|PLAN|SHA256|EMAIL|OKU", n, re.I) for n in z.namelist()),
      "zip holds no key, analysis, plan, hash, e-mail or coordinator file")
check(MZ.build_bytes() == open(zp, "rb").read(), "make_zip_kit2.py rebuilds the zip byte-identically")
print(f"INFO  zip SHA-256 {sha256(zp)}")
print(f"INFO  kit2_key.csv SHA-256 {sha256(B2.KEY2)}")
print("ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED")
sys.exit(0 if ok_all else 1)
