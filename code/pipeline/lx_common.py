"""lx_common.py - shared helpers for the LEAD-chi processing pipeline (task P5b, 2026-09-25).

Written, hashed and registered before any QSM is computed. NOT executed on real data before the OSF
registration. Imported by blind_setup.py, qsm_plausibility.py, extract_rois.py, qc_montage.py, build_usable.py,
lock_labels.py and unlock_and_merge.py (the last two use only the standard-library parts).

Contents
  * paths and guards: every write must stay inside $LX; pipeline scripts refuse label files and, after setup,
    raw (uncoded) ds007116 inputs
  * hashing and an append-only, hash-chained run log ($LX/logs/pipeline_runs.jsonl)
  * OpenSSL hybrid encryption (RSA-OAEP key wrap + AES-256-CBC/PBKDF2), plaintext only in memory
  * NIfTI header scrubbing (descrip, aux_file, intent_name, db_name, extensions) without resampling
  * ROI construction in the native MEGRE grid, weighted statistics, Q_phase and head-orientation angles
"""
from __future__ import annotations

import datetime
import getpass
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys

LX_DEFAULT = "/path/to/Penn_LEAD/LEADchi"
PL_DEFAULT = "/path/to/Penn_LEAD"
TE_REGISTERED_S = (0.0075, 0.015, 0.0225, 0.030)   # protocol sec. 12.1 step 0
GAMMA_QSMXT_MHZ_T = 42.576                          # QSMxT internal gamma (E3 finding A1)
SEED = 20260924

# --------------------------------------------------------------------------------------------- paths and guards

def lx_root() -> str:
    return os.path.realpath(os.environ.get("LX", LX_DEFAULT))


def pl_root() -> str:
    return os.path.realpath(os.environ.get("PL", PL_DEFAULT))


def utcnow() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def is_inside(path: str, root: str) -> bool:
    p = os.path.realpath(path)
    r = os.path.realpath(root)
    return p == r or p.startswith(r + os.sep)


def assert_inside(path: str, root: str | None = None) -> str:
    """Abort unless `path` resolves inside `root` (default $LX). Used before every write."""
    root = root or lx_root()
    if not is_inside(path, root):
        sys.exit(f"REFUSED: write target {path} is outside {root}")
    return path


def ensure_parent(path: str, root: str | None = None, dry: bool = False) -> str:
    assert_inside(path, root)
    d = os.path.dirname(os.path.abspath(path))
    if not dry:
        os.makedirs(d, exist_ok=True)
    return path


_LABEL_RE = re.compile(r"(participants\.(tsv|json)$|_sessions\.(tsv|json)$|/phenotype(/|$))")


def forbid_label_inputs(paths) -> None:
    """Pipeline scripts must never read group labels, ages or phenotypes (v1.2 sec. 2.5)."""
    for p in paths:
        if p and _LABEL_RE.search(os.path.abspath(str(p))):
            sys.exit(f"REFUSED: {p} is a label/age/phenotype file; only lock_labels.py and "
                     "unlock_and_merge.py may read it")


def forbid_raw_inputs(paths) -> None:
    """After blind setup, pipeline scripts read only the coded tree (written commitment, E1 sec. 6.4(5))."""
    raw = os.path.join(pl_root(), "ds007116")
    for p in paths:
        if p and is_inside(str(p), raw):
            sys.exit(f"REFUSED: {p} is inside the uncoded dataset {raw}; after setup only $LX/blind is read")


# --------------------------------------------------------------------------------------------- hashing and logs

def sha256_file(path: str, bufsize: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def say(msg: str, dry: bool = False) -> None:
    print(("[DRY-RUN] " if dry else "") + msg, flush=True)


class RunRecord:
    """One JSON line per script run in $LX/logs/pipeline_runs.jsonl, chained by the SHA-256 of the previous line.

    Records host, user, UTC times, the SHA-256 of the script and of lx_common.py, argv, and the SHA-256 of the
    registered input/output tables. Nothing is written in dry-run mode."""

    def __init__(self, script_path: str, dry: bool, log_dir: str | None = None):
        self.dry = dry
        self.log_path = os.path.join(log_dir or os.path.join(lx_root(), "logs"), "pipeline_runs.jsonl")
        self.rec = {
            "utc_start": utcnow(), "host": socket.gethostname(), "user": getpass.getuser(),
            "script": os.path.basename(script_path), "script_sha256": sha256_file(script_path),
            "lx_common_sha256": sha256_file(os.path.abspath(__file__)), "argv": sys.argv[1:],
            "dry_run": dry, "inputs": {}, "outputs": {}, "notes": [],
        }

    def add_input(self, path: str, hash_it: bool = True) -> None:
        self.rec["inputs"][path] = sha256_file(path) if (hash_it and os.path.isfile(path)) else None

    def add_output(self, path: str, hash_it: bool = True) -> None:
        self.rec["outputs"][path] = sha256_file(path) if (hash_it and os.path.isfile(path)) else None

    def note(self, s: str) -> None:
        self.rec["notes"].append(s)

    def finish(self, status: str = "ok") -> None:
        self.rec["utc_end"] = utcnow()
        self.rec["status"] = status
        if self.dry:
            say(f"run record (not written): script={self.rec['script']} status={status}", True)
            return
        ensure_parent(self.log_path)
        prev = None
        if os.path.exists(self.log_path):
            with open(self.log_path, "rb") as f:
                lines = f.read().splitlines()
            if lines:
                prev = sha256_bytes(lines[-1])
        self.rec["prev_line_sha256"] = prev
        with open(self.log_path, "a", encoding="utf-8") as f:
            f.write(json.dumps(self.rec, sort_keys=True) + "\n")


def _json_clean(o):
    """NaN/inf -> null and numpy scalars -> Python, so every JSON file is strict (R jsonlite, P5 gates)."""
    if isinstance(o, dict):
        return {str(k): _json_clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_clean(v) for v in o]
    if hasattr(o, "item") and not isinstance(o, (str, bytes)):
        try:
            o = o.item()
        except (ValueError, AttributeError):
            pass
    if isinstance(o, float) and (o != o or o in (float("inf"), float("-inf"))):
        return None
    return o


def write_json(path: str, obj, dry: bool = False) -> None:
    ensure_parent(path, dry=dry)
    if dry:
        say(f"would write {path}", True)
        return
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_json_clean(obj), f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")


def write_csv(path: str, rows: list[dict], columns: list[str], dry: bool = False) -> None:
    import csv
    ensure_parent(path, dry=dry)
    if dry:
        say(f"would write {path} ({len(rows)} rows x {len(columns)} columns)", True)
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k)) for k in columns})


def read_csv(path: str) -> list[dict]:
    import csv
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t" if path.endswith(".tsv") else ","))


def fmt(v):
    """CSV cell formatting: None/NaN/inf -> NA, booleans -> 0/1, floats rounded to 8 decimals."""
    if v is None:
        return "NA"
    if hasattr(v, "item") and not isinstance(v, (str, bytes)):
        try:
            v = v.item()
        except (ValueError, AttributeError):
            pass
    if isinstance(v, bool):
        return int(v)
    if isinstance(v, float):
        if v != v or v in (float("inf"), float("-inf")):
            return "NA"
        return repr(round(v, 8))
    return v


def load_config(path: str | None = None) -> dict:
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "atlas_config.json")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def lx_path(rel: str) -> str:
    """Resolve a config path: ${VAR} expanded; relative paths are relative to $LX."""
    p = os.path.expandvars(rel)
    return p if os.path.isabs(p) else os.path.join(lx_root(), p)


# --------------------------------------------------------------------------------------------- OpenSSL hybrid

def openssl_bin() -> str:
    b = os.environ.get("OPENSSL") or shutil.which("openssl")
    if not b:
        sys.exit("openssl not found (set $OPENSSL)")
    v = subprocess.run([b, "version"], capture_output=True, text=True).stdout.strip()
    m = re.match(r"OpenSSL (\d+)\.(\d+)\.(\d+)", v)
    if not m or (int(m.group(1)), int(m.group(2)), int(m.group(3))) < (1, 1, 1):
        sys.exit(f"need OpenSSL >= 1.1.1 (for enc -pbkdf2); found: {v!r}")
    return b


def rsa_bits(pub_pem: str) -> int:
    out = subprocess.run([openssl_bin(), "pkey", "-pubin", "-in", pub_pem, "-noout", "-text"],
                         capture_output=True, text=True, check=True).stdout
    m = re.search(r"\((\d+) bit\)", out)
    return int(m.group(1)) if m else 0


def _enc_with_pass(args: list[str], pw: bytes, data: bytes | None) -> bytes:
    """Run `openssl enc` with the pass phrase on a private pipe fd (never on the command line)."""
    r, w = os.pipe()
    try:
        os.write(w, pw + b"\n")
        os.close(w)
        res = subprocess.run(args + ["-pass", f"fd:{r}"], input=data, pass_fds=(r,),
                             stdout=subprocess.PIPE, check=True)
    finally:
        try:
            os.close(r)
        except OSError:
            pass
    return res.stdout


def hybrid_encrypt(plaintext: bytes, pub_pem: str, out_data: str, out_key: str, dry: bool = False) -> dict:
    """Encrypt `plaintext` (held in memory) to `pub_pem`.

    A fresh 256-bit random session secret is wrapped with RSA-OAEP (openssl pkeyutl) and used as the PBKDF2 pass
    phrase for AES-256-CBC (openssl enc -pbkdf2 -iter 200000). Only ciphertext touches the disk."""
    import secrets
    if dry:
        say(f"would encrypt {len(plaintext)} bytes to {out_data} (+ wrapped key {out_key})", True)
        return {}
    ossl = openssl_bin()
    for p in (out_data, out_key):
        ensure_parent(p)
        if os.path.exists(p):
            sys.exit(f"REFUSED: {p} exists (locked files are written once)")
    pw = secrets.token_hex(32).encode()
    subprocess.run([ossl, "pkeyutl", "-encrypt", "-pubin", "-inkey", pub_pem,
                    "-pkeyopt", "rsa_padding_mode:oaep", "-out", out_key], input=pw, check=True)
    ct = _enc_with_pass([ossl, "enc", "-aes-256-cbc", "-salt", "-pbkdf2", "-iter", "200000", "-md", "sha256"],
                        pw, plaintext)
    with open(out_data, "wb") as f:
        f.write(ct)
    del pw
    return {"data": out_data, "data_sha256": sha256_file(out_data),
            "key": out_key, "key_sha256": sha256_file(out_key)}


def hybrid_decrypt(in_data: str, in_key: str, priv_pem: str) -> bytes:
    """Decrypt to memory. openssl itself prompts for the private-key pass phrase on the terminal."""
    ossl = openssl_bin()
    pw = subprocess.run([ossl, "pkeyutl", "-decrypt", "-inkey", priv_pem, "-pkeyopt", "rsa_padding_mode:oaep",
                         "-in", in_key], stdout=subprocess.PIPE, check=True).stdout.strip()
    with open(in_data, "rb") as f:
        ct = f.read()
    pt = _enc_with_pass([ossl, "enc", "-d", "-aes-256-cbc", "-pbkdf2", "-iter", "200000", "-md", "sha256"],
                        pw, ct)
    del pw
    return pt


# --------------------------------------------------------------------------------------------- NIfTI helpers
# numpy / nibabel / scipy are imported lazily so that the stdlib-only scripts can import this module.

SCRUB_TEXT = ("descrip", "aux_file", "intent_name", "db_name", "data_type")


def scrubbed_header(img, dtype):
    import numpy as np
    hdr = img.header.copy()
    for f in SCRUB_TEXT:
        if f in hdr:
            hdr[f] = b""
    for f in ("intent_code", "intent_p1", "intent_p2", "intent_p3", "glmax", "glmin", "toffset",
              "session_error", "extents", "cal_min", "cal_max"):
        if f in hdr:
            hdr[f] = 0
    del hdr.extensions[:]
    hdr.set_data_dtype(np.dtype(dtype))
    hdr.set_slope_inter(1.0, 0.0)
    return hdr


def write_scrubbed(img, data, out_path: str, dry: bool = False, offset_mm=None) -> None:
    """Write `data` with the geometry (qform/sform, zooms, units) of `img` and a scrubbed header. No resampling.

    offset_mm: optional per-session world translation (RAS mm) added to BOTH the qform and the sform (identical for
    every image of a session, so relative geometry, rotations and the B0 direction are unchanged). This defeats a
    header-geometry match between coded and raw files (registration Annex A.3 item 5)."""
    import nibabel as nib
    import numpy as np
    ensure_parent(out_path, dry=dry)
    if dry:
        say(f"would write {out_path} {tuple(data.shape)} {data.dtype}", True)
        return
    hdr = scrubbed_header(img, data.dtype)
    target = img.affine.copy()
    if offset_mm is not None:
        target[:3, 3] += np.asarray(offset_mm, float)
        qc, sc = int(img.header["qform_code"]), int(img.header["sform_code"])
        hdr.set_sform(target, code=sc if sc > 0 else 1)
        hdr.set_qform(target, code=qc if qc > 0 else 1)
    new = nib.Nifti1Image(data, None, header=hdr)
    if not np.allclose(new.header.get_best_affine(), target, atol=1e-4):
        sys.exit(f"geometry changed while writing {out_path}; aborting")
    for k in ("qform_code", "sform_code"):
        if offset_mm is None and int(new.header[k]) != int(img.header[k]):
            sys.exit(f"{k} changed while writing {out_path}; aborting")
    nib.save(new, out_path)


def header_text_blob(path: str) -> bytes:
    """Raw header bytes (up to vox_offset, i.e. including extensions) of a .nii/.nii.gz, for the ID-leak scan."""
    import gzip
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as f:
        head = f.read(352)
        import struct
        vox_offset = int(struct.unpack("<f", head[108:112])[0]) if len(head) >= 112 else 352
        more = f.read(max(0, min(vox_offset, 1 << 22) - 352))
    return head + more


def load_nii(path: str, dtype="float32"):
    import nibabel as nib
    import numpy as np
    img = nib.load(path)
    return img, np.asarray(img.get_fdata(dtype=np.float64), dtype=dtype)


def load_labels(path: str):
    import nibabel as nib
    import numpy as np
    img = nib.load(path)
    a = np.asarray(img.dataobj)
    if not np.issubdtype(a.dtype, np.integer):
        r = np.rint(a)
        if np.nanmax(np.abs(a - r)) > 1e-3:
            sys.exit(f"{path}: label image is not integer-valued")
        a = r.astype(np.int32)
    return img, a.astype(np.int32)


def same_grid(img_a, img_b, tol: float = 1e-4) -> bool:
    import numpy as np
    return img_a.shape[:3] == img_b.shape[:3] and np.allclose(img_a.affine, img_b.affine, atol=tol)


# --------------------------------------------------------------------------------------------- morphology

def _st6():
    from scipy import ndimage
    return ndimage.generate_binary_structure(3, 1)   # 6-connectivity (face neighbours)


def erode(mask, n: int = 1):
    from scipy import ndimage
    if n <= 0:
        return mask.astype(bool)
    return ndimage.binary_erosion(mask, structure=_st6(), iterations=n, border_value=0)


def dilate(mask, n: int = 1):
    from scipy import ndimage
    if n <= 0:
        return mask.astype(bool)
    return ndimage.binary_dilation(mask, structure=_st6(), iterations=n)


def fill_holes(mask):
    from scipy import ndimage
    return ndimage.binary_fill_holes(mask)


def voxel_volume_mm3(img) -> float:
    import numpy as np
    return float(np.prod(img.header.get_zooms()[:3]))


# --------------------------------------------------------------------------------------------- ROI construction
# Registered ROI definitions (protocol v1.1 sec. 6 and 12.2; v1.2 sec. 5), all in the native MEGRE grid:
#   atlas ROIs (CIT168, dentate atlas): w = min(1, sum of component probabilities); support = w >= 0.5 and more
#       than 1 voxel from SynthSeg CSF (CSF plus its 6-neighbourhood excluded); probability-weighted mean.
#       P5b choice: the same 1-voxel CSF buffer is applied to every atlas ROI (dentate, RN, STN), not only SN-VTA.
#   SynthSeg ROIs: each label eroded by 1 voxel (6-connectivity) in the MEGRE grid, NAcc not eroded; equal
#       voxel weights (volume-weighted means of unions).
#   Mask 3 restriction and the implausible-voxel exclusion are applied by roi_stats(), not here.

ATLAS_ROIS = {
    # name: list of (file key, component) summed then clipped at 1
    "SNVTA": [("cit", "SNc"), ("cit", "SNr"), ("cit", "VTA"), ("cit", "PBP")],
    "SNonly": [("cit", "SNc"), ("cit", "SNr")],
    "VTAPBP": [("cit", "VTA"), ("cit", "PBP")],
    "RN": [("cit", "RN")],
    "STN": [("cit", "STH")],
    "DN_R": [("dn", "DN_R")],
    "DN_L": [("dn", "DN_L")],
    "DNIP_R": [("dn", "DN_R"), ("dn", "IP_R")],
}


def seg_rois_def(ss: dict) -> dict:
    """SynthSeg ROIs: name -> (eroded labels, non-eroded labels)."""
    return {
        "GP": ([ss["L_GP"], ss["R_GP"]], []),
        "PUT": ([ss["L_PUT"], ss["R_PUT"]], []),
        "PUT_R": ([ss["R_PUT"]], []),
        "CAU": ([ss["L_CAU"], ss["R_CAU"]], []),
        "THAL": ([ss["L_THAL"], ss["R_THAL"]], []),
        "STR": ([ss["L_CAU"], ss["R_CAU"], ss["L_PUT"], ss["R_PUT"]], [ss["L_ACC"], ss["R_ACC"]]),
        "DGM": ([ss["L_CAU"], ss["R_CAU"], ss["L_GP"], ss["R_GP"]], []),
        "WM": ([ss["L_WM"], ss["R_WM"]], []),
        "CWM": ([ss["L_CWM"], ss["R_CWM"]], []),
        "LV1": ([ss["L_LV"], ss["R_LV"]], []),
    }


def build_rois(reg_dir: str, cfg: dict) -> dict:
    """Return {'rois': {name: {'w': float32, 'sup': bool, 'kind': str}}, 'aseg', 'csf', 'csfbuf', 'img', 'vox'}.

    Reads only registration outputs written by register_rois.sh (atlas warps and SynthSeg in the MEGRE grid)."""
    import numpy as np
    ss = cfg["synthseg"]
    aimg, aseg = load_labels(os.path.join(reg_dir, "aseg_megre.nii.gz"))
    csf = np.isin(aseg, ss["csf_labels"])
    csfbuf = dilate(csf, 1)
    comp = {}
    cimg, cit = load_nii(os.path.join(reg_dir, "cit_megre.nii.gz"))
    for i, nm in enumerate(cfg["cit168"]["subset_order"]):
        comp[("cit", nm)] = cit[..., i]
    dimg, dn = load_nii(os.path.join(reg_dir, "dn_megre.nii.gz"))
    for i, nm in enumerate(cfg["dentate"]["subset_order"]):
        comp[("dn", nm)] = dn[..., i]
    for im in (cimg, dimg):
        if not same_grid(im, aimg):
            sys.exit(f"{reg_dir}: atlas and SynthSeg warps are not on the same MEGRE grid")
    rois = {}
    for name, parts in ATLAS_ROIS.items():
        w = np.clip(sum(comp[p] for p in parts), 0.0, 1.0).astype("float32")
        rois[name] = {"w": w, "sup": (w >= 0.5) & ~csfbuf, "sup_atlas": (w >= 0.5), "kind": "atlas"}
    for name, (er, ne) in seg_rois_def(ss).items():
        m = np.zeros(aseg.shape, bool)
        for lab in er:
            m |= erode(aseg == lab, 1)
        for lab in ne:
            m |= (aseg == lab)
        rois[name] = {"w": m.astype("float32"), "sup": m, "sup_atlas": m, "kind": "seg"}
    pimg, plic = load_labels(os.path.join(reg_dir, "plic_megre.nii.gz"))
    pl = plic > 0          # jhu_plic subset file is binary (1 = PLIC_R or PLIC_L), see register_rois.sh
    rois["PLIC"] = {"w": pl.astype("float32"), "sup": pl & ~csfbuf, "sup_atlas": pl, "kind": "atlas"}
    return {"rois": rois, "aseg": aseg, "csf": csf, "csfbuf": csfbuf, "img": aimg, "vox": voxel_volume_mm3(aimg)}


def support_mask(m3, chi, qsmxt_mask=None):
    """Support S of one reconstruction (registered text sec. 12.1): Mask 3 AND the reconstruction's own
    background-removal/inversion mask (QSMxT `_mask.nii`, when given) AND finite, non-zero chi. The zero rim that
    V-SHARP erosion and MEDI leave inside Mask 3 is therefore never averaged into the reference."""
    import numpy as np
    S = (m3 > 0.5) & np.isfinite(chi) & (chi != 0)
    if qsmxt_mask is not None:
        S &= qsmxt_mask > 0.5
    return S


def two_step_reference(chi_ppm, S, thr_ppm: float):
    """Registered two-step whole-brain reference: ref0 = mean chi over S; implausible voxels = S & |chi - ref0| >
    thr; ref = mean over S minus implausible voxels. Returns (ref_ppm, implausible_mask, ref0_ppm)."""
    import numpy as np
    if not S.any():
        return float("nan"), np.zeros(S.shape, bool), float("nan")
    ref0 = float(chi_ppm[S].mean())
    imp = S & (np.abs(chi_ppm - ref0) > thr_ppm)
    keep = S & ~imp
    ref = float(chi_ppm[keep].mean()) if keep.any() else float("nan")
    return ref, imp, ref0


def roi_stats(x, valid, imp, roi: dict) -> dict:
    """Statistics of map `x` over one ROI.

    valid: the support S (see support_mask; for R2* maps: Mask 3 & finite); imp: implausible voxels (or None).
    Returns the primary probability-weighted mean (implausible excluded), the mean with implausible voxels retained,
    the unweighted median (implausible excluded), the in-support weight share and the implausible weight share
    (schema: inmask_*, fimp_*; "inmask" is measured against S, not against Mask 3)."""
    import numpy as np
    w, sup = roi["w"], roi["sup"]
    wt = float(w[sup].sum())
    inm = sup & valid
    w_in = float(w[inm].sum())
    use = inm & ~imp if imp is not None else inm
    w_use = float(w[use].sum())
    out = {"sumw": wt, "inmask": (w_in / wt) if wt > 0 else float("nan"),
           "fimp": (float(w[inm & imp].sum()) / w_in) if (imp is not None and w_in > 0) else
           (0.0 if w_in > 0 else float("nan")),
           "mean": float((w[use] * x[use]).sum() / w_use) if w_use > 0 else float("nan"),
           "mean_impkeep": float((w[inm] * x[inm]).sum() / w_in) if w_in > 0 else float("nan"),
           "median": float(np.median(x[use])) if use.any() else float("nan"),
           "nvox": int(use.sum())}
    return out


def csf_reference(chi_ref_ppb, valid, imp, aseg, r2s, cfg: dict):
    """Secondary (CSF-only) reference (registered text sec. 12.2). Lateral ventricles (4/43) eroded by 2 voxels
    (fallback 1), restricted to the support S and plausible voxels, choroid plexus removed by the chi-independent rule
    R2* > ventricle median + 2 SD; >= 200 voxels required; valid only if |mean CSF - whole-brain mean| <= 50 ppb.
    Returns (type, mean_ppb_relative_to_whole_brain, n_voxels); type in {'CSF', 'CSF_erode1', 'none'}. There is NO
    PLIC fallback: a session without a valid CSF reference is missing for the CSF-referenced quantities."""
    import numpy as np
    c = cfg["csf_ref"]
    ss = cfg["synthseg"]
    lv = np.isin(aseg, [ss["L_LV"], ss["R_LV"]])
    for n_er, typ in ((c["erode_primary"], "CSF"), (c["erode_fallback"], "CSF_erode1")):
        m = erode(lv, n_er) & valid
        if imp is not None:
            m &= ~imp
        if r2s is not None and m.sum() > 0:
            rv = r2s[m]
            rv = rv[np.isfinite(rv)]
            if rv.size > 1:
                thr = np.median(rv) + c["choroid_r2s_sd"] * np.std(rv, ddof=1)
                m &= ~(np.nan_to_num(r2s, nan=np.inf) > thr)
        n = int(m.sum())
        if n >= c["min_voxels"]:
            mean = float(chi_ref_ppb[m].mean())
            if abs(mean) > c["max_abs_csf_minus_wb_ppb"]:
                return "none", mean, n
            return typ, mean, n
    return "none", float("nan"), 0


def tissue_fractions(aseg, S, imp, cfg: dict) -> dict:
    """SynthSeg WM / GM / CSF fractions of the reference domain S minus implausible voxels (label-free covariates of
    the reference composition; registered text sec. 12.2)."""
    import numpy as np
    ss = cfg["synthseg"]
    dom = S & ~imp if imp is not None else S
    n = int(dom.sum())
    if n == 0:
        return {"ref_frac_wm": float("nan"), "ref_frac_gm": float("nan"), "ref_frac_csf": float("nan")}
    wm = np.isin(aseg, [ss["L_WM"], ss["R_WM"], ss["L_CWM"], ss["R_CWM"], ss["BS"], ss["L_VDC"], ss["R_VDC"]])
    csf = np.isin(aseg, ss["csf_labels"])
    gm = (aseg > 0) & ~wm & ~csf
    return {"ref_frac_wm": float((wm & dom).sum() / n), "ref_frac_gm": float((gm & dom).sum() / n),
            "ref_frac_csf": float((csf & dom).sum() / n)}


def otsu_threshold(values, nbins: int = 256) -> float:
    """Otsu (1979) threshold of a 1-D sample (numpy only)."""
    import numpy as np
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    hist, edges = np.histogram(v, bins=nbins)
    mids = 0.5 * (edges[:-1] + edges[1:])
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids)
    mu0 = m0 / np.maximum(w0, 1)
    mu1 = (m0[-1] - m0) / np.maximum(w1, 1)
    between = w0 * w1 * (mu0 - mu1) ** 2
    return float(mids[int(np.argmax(between[:-1]))])


def trim_partition_edges(mask, n_partitions: int, n_trim: int = 2):
    """Remove the first and last `n_trim` planes along the 3D-slab partition axis (the unique voxel axis whose length
    equals the registered number of partitions, 96) from a mask (registered text Annex B.2)."""
    import numpy as np
    ax = [i for i, s in enumerate(mask.shape[:3]) if s == n_partitions]
    if len(ax) != 1:
        sys.exit(f"cannot identify the partition axis: shape {mask.shape}, partitions {n_partitions}")
    m = mask.copy()
    sl = [slice(None)] * 3
    sl[ax[0]] = slice(0, n_trim)
    m[tuple(sl)] = False
    sl[ax[0]] = slice(mask.shape[ax[0]] - n_trim, None)
    m[tuple(sl)] = False
    return m


def ghost_to_signal(mag, pe_axis: int, dil: int = 3) -> float:
    """Magnitude ghost-to-signal ratio (label-free motion index; registered text sec. 12.2): head = echo-1 magnitude
    above its Otsu threshold, holes filled; ghost band = head shifted by half the FOV along the phase-encode axis,
    minus the head dilated by `dil` voxels; background = the rest outside the dilated head.
    GSR = (mean(ghost) - mean(background)) / mean(head)."""
    import numpy as np
    head = fill_holes(mag > otsu_threshold(mag[np.isfinite(mag)]))
    hd = dilate(head, dil)
    ghost = np.roll(head, mag.shape[pe_axis] // 2, axis=pe_axis) & ~hd
    bg = ~hd & ~ghost
    if not ghost.any() or not bg.any() or not head.any():
        return float("nan")
    return float((mag[ghost].mean() - mag[bg].mean()) / mag[head].mean())


# --------------------------------------------------------------------------------------------- Q_phase

def q_phase(phases, mask) -> float:
    """Registered phase-noise/linearity index (sec. 12.2; NOT a motion index): log of the median absolute multi-echo phase-prediction
    residual in Mask 3; echo 3 predicted as 2*phi2 - phi1 and echo 4 as 3*phi2 - 2*phi1, residuals wrapped to
    (-pi, pi]. Independent of tissue chi because phase is linear in TE for any static field (TE_n = n * dTE)."""
    import numpy as np
    p1, p2, p3, p4 = phases
    r3 = np.angle(np.exp(1j * (p3[mask] - 2.0 * p2[mask] + p1[mask])))
    r4 = np.angle(np.exp(1j * (p4[mask] - 3.0 * p2[mask] + 2.0 * p1[mask])))
    med = float(np.median(np.abs(np.concatenate([r3, r4]))))
    return float(np.log(med)) if med > 0 else float("-inf")


# --------------------------------------------------------------------------------------------- orientation

def read_itk_affine(mat_path: str):
    """ITK/ANTs .mat (MatrixOffsetTransformBase, AffineTransform or Euler3D) -> 4x4 LPS matrix mapping FIXED
    points to MOVING points (ITK convention)."""
    import numpy as np
    import scipy.io as sio
    m = sio.loadmat(mat_path)
    key = [k for k in m if k != "fixed" and not k.startswith("__")][0]
    p = np.asarray(m[key], float).ravel()
    c = np.asarray(m["fixed"], float).ravel()[:3]
    if p.size == 12:
        A = p[:9].reshape(3, 3)
        t = p[9:12]
    elif p.size == 6:   # Euler3D: angles x, y, z then translation; ITK ZXY composition
        ax, ay, az = p[:3]
        t = p[3:6]
        Rx = np.array([[1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])
        Ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
        Rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
        A = Rz @ Rx @ Ry
    else:
        sys.exit(f"{mat_path}: unsupported ITK transform with {p.size} parameters")
    T = np.eye(4)
    T[:3, :3] = A
    T[:3, 3] = t + c - A @ c
    return T


def orientation_angles(rot_chain_lps) -> dict:
    """Head orientation relative to B0 (protocol sec. 12.1 step 13).

    rot_chain_lps: list of 4x4 ITK matrices, each mapping points of the head (reference) frame towards the MEGRE
    (scanner) frame, composed in order [head frame -> ... -> MEGRE]. The registered head frame is the same for every
    session: MNI152NLin2009cAsym through the affine stage of the route-R T1w->template registration. B0 is the
    scanner z axis (the MEGRE world +z; the obliquity of the slab is checked in the header audit). The B0 unit vector is expressed in the head frame (RAS):
      pitch = atan2(b_y, b_z)  rotation about the left-right axis (sign: + = B0 tilted towards anterior)
      roll  = atan2(b_x, b_z)  rotation about the anterior-posterior axis
      tilt  = arccos(b_z)      total angle between B0 and the head superior axis
      yaw   = rotation of the head about B0 (does not change the B0 angle; reported for completeness)
    abs_pitch = |pitch| (S3 covariate). Rotations are taken from the polar decomposition, so small scaling in
    an affine frame does not bias the angles."""
    import numpy as np
    from scipy.linalg import polar
    M = np.eye(3)
    for T in rot_chain_lps:          # x_moving = A x_fixed; compose head -> scanner
        M = np.asarray(T)[:3, :3] @ M
    R, _ = polar(M)
    L = np.diag([-1.0, -1.0, 1.0])
    R_ras = L @ R @ L                # head(RAS) -> scanner(RAS) direction map
    b = np.linalg.inv(R_ras) @ np.array([0.0, 0.0, 1.0])
    b = b / np.linalg.norm(b)
    Rh = np.linalg.inv(R_ras)        # scanner -> head
    return {"pitch_deg": float(np.degrees(np.arctan2(b[1], b[2]))),
            "roll_deg": float(np.degrees(np.arctan2(b[0], b[2]))),
            "b0_tilt_deg": float(np.degrees(np.arccos(np.clip(b[2], -1, 1)))),
            "yaw_deg": float(np.degrees(np.arctan2(Rh[1, 0], Rh[0, 0]))),
            "abs_pitch": float(abs(np.degrees(np.arctan2(b[1], b[2]))))}


def dice(a, b) -> float:
    import numpy as np
    s = int(a.sum()) + int(b.sum())
    return float(2.0 * np.logical_and(a, b).sum() / s) if s > 0 else float("nan")


# --------------------------------------------------------------------------------------------- reconstruction outputs

def variant_outputs(lx: str, variant: str, code: str, ses: int, chi_pattern: str | None = None,
                    mask_pattern: str | None = None) -> dict:
    """Deterministic file resolution for one reconstruction of one coded session.

    QSMxT (default): chi = derivatives/qsmxt/sub-C/ses-K/anat/<lab>*_Chimap.nii* excluding any `desc-` variant;
    mask = the same folder's <lab>*_mask.nii* excluding `desc-` variants (the reconstruction's own final mask);
    localfield = exactly one <lab>*_localfield.nii* under derivatives/qsmxt (ppm; QSMxT runner `local_field_ppm`).
    Non-QSMxT engines pass glob patterns ({lx} {code} {ses} {lab} {variant}). Anything ambiguous returns []."""
    import glob as _g
    lab = sess_label(code, ses)
    base = os.path.join(lx, "qsm", variant, "derivatives", "qsmxt")
    anat = os.path.join(base, f"sub-{code}", f"ses-{ses}", "anat")
    fmt_ = dict(lx=lx, code=code, ses=ses, lab=lab, variant=variant)
    nodesc = lambda xs: [x for x in xs if "desc-" not in os.path.basename(x)]
    chi = sorted(_g.glob(chi_pattern.format(**fmt_))) if chi_pattern else nodesc(sorted(_g.glob(os.path.join(anat, f"{lab}*_Chimap.nii*"))))
    if mask_pattern:
        mask = sorted(_g.glob(mask_pattern.format(**fmt_)))
    elif chi_pattern:
        mask = []
    else:
        mask = nodesc(sorted(_g.glob(os.path.join(anat, f"{lab}*_mask.nii*"))))
    lf = nodesc(sorted(_g.glob(os.path.join(base, "**", f"{lab}*_localfield.nii*"), recursive=True)))
    return {"lab": lab, "chi": chi, "mask": mask, "localfield": lf if len(lf) == 1 else []}


def mask_assertion(qmask, m3) -> dict:
    """Post-run check that the reconstruction used Mask 3 (registered text Annex B.2): >= 99 % of the
    reconstruction's mask lies inside Mask 3, and it covers >= 90 % of Mask 3 eroded by 2 voxels (V-SHARP erodes a
    thin rim). Otherwise the session is invalid for that run (status 'mask_mismatch')."""
    q = qmask > 0.5
    m = m3 > 0.5
    inside = float((q & m).sum() / max(1, q.sum()))
    core = erode(m, 2)
    cover = float((q & core).sum() / max(1, core.sum()))
    return {"mask_inside_m3": inside, "mask_covers_m3core": cover, "mask_ok": inside >= 0.99 and cover >= 0.90}


# --------------------------------------------------------------------------------------------- coded tree

def coded_sessions(lx: str | None = None) -> list[tuple[str, int]]:
    """(code, ses_order) pairs in the randomised processing order written by blind_setup.py."""
    lx = lx or lx_root()
    p = os.path.join(lx, "blind", "processing_order.txt")
    if not os.path.exists(p):
        sys.exit(f"{p} missing: run blind_setup.py first")
    out = []
    with open(p, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                c, s = line.split("\t")
                out.append((c, int(s)))
    return out


def sess_label(code: str, ses: int) -> str:
    return f"sub-{code}_ses-{ses}"


def pipeline_choices(lx: str | None = None) -> dict:
    p = os.path.join(lx or lx_root(), "prereg", "pipeline_choices.json")
    if not os.path.exists(p):
        sys.exit(f"{p} missing: run `run_qsmxt.sh --stage config` first")
    with open(p, encoding="utf-8") as f:
        return json.load(f)
