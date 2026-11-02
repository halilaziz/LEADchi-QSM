"""LEAD-chi POST HOC recovery test (2026-10-07): inject a known susceptibility change into the SN-VTA ROI.

For each coded session, a uniform susceptibility change `level` (ppb) is placed in the SN-VTA ROI voxel set actually
averaged by extract_rois.py (CIT168 SNc+SNr+VTA+PBP, summed probability >= 0.5, outside the 1-voxel CSF buffer;
lx_common.build_rois()["rois"]["SNVTA"]["sup"]). Its field is computed with the k-space dipole kernel
D = 1/3 - (k.b)^2/|k|^2 (B0 direction from the affine; zero-padded 2x; D(0) = 0), converted to phase per echo
(dphi = sign * 2*pi * f0 * dB[ppm] * 1e-6 * TE, f0 = MagneticFieldStrength * 42.576e6 Hz as in the coded sidecar) and
added to the measured phase (re-wrapped). Magnitudes, sidecars and Mask 3 are linked unchanged, so the re-run differs
from the primary run only by the injected field. Nothing outside $LX/posthoc and $LX/qsm/<new variant> is written.
"""
import argparse, json, os, sys
from concurrent.futures import ProcessPoolExecutor

LX = os.environ["LX"]
sys.path.insert(0, os.environ.get("LX_SCRIPTS", os.path.join(LX, "scripts", "pipeline")))
import lx_common as L  # noqa: E402


def dipole_field(dchi, vox, bdir):
    import numpy as np
    sh = np.array(dchi.shape)
    pad = [int(s) for s in sh]                      # pad to 2x in every dimension
    big = np.zeros(tuple(sh + pad), dtype=np.float64)
    big[: sh[0], : sh[1], : sh[2]] = dchi
    ks = [np.fft.fftfreq(n, d=v) for n, v in zip(big.shape, vox)]
    kx, ky, kz = np.meshgrid(*ks, indexing="ij")
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    kb = kx * bdir[0] + ky * bdir[1] + kz * bdir[2]
    with np.errstate(invalid="ignore", divide="ignore"):
        D = 1.0 / 3.0 - kb ** 2 / k2
    D[k2 == 0] = 0.0
    f = np.real(np.fft.ifftn(D * np.fft.fftn(big)))
    return f[: sh[0], : sh[1], : sh[2]]


def one(job):
    import numpy as np, nibabel as nib
    code, ses, level, sign, out_bids, src_bids, cfg = job
    lab = L.sess_label(code, ses)
    rd = os.path.join(LX, "reg", f"sub-{code}", f"ses-{ses}")
    src = os.path.join(src_bids, f"sub-{code}", f"ses-{ses}", "anat")
    dst = os.path.join(out_bids, f"sub-{code}", f"ses-{ses}", "anat")
    if not os.path.exists(os.path.join(rd, "cit_megre.nii.gz")):
        return {"code": code, "ses": ses, "status": "no_registration"}
    R = L.build_rois(rd, cfg)
    sup = R["rois"]["SNVTA"]["sup"]
    p1 = os.path.join(src, f"{lab}_echo-1_part-phase_MEGRE.nii.gz")
    img = nib.load(p1)
    if img.shape[:3] != sup.shape:
        return {"code": code, "ses": ses, "status": "grid_mismatch"}
    aff = img.affine
    vox = np.sqrt((aff[:3, :3] ** 2).sum(axis=0))
    rot = aff[:3, :3] / vox
    bdir = np.linalg.solve(rot, np.array([0.0, 0.0, 1.0]))   # scanner z expressed in voxel axes
    bdir = bdir / np.linalg.norm(bdir)
    dchi = np.zeros(sup.shape, dtype=np.float64)
    dchi[sup] = level * 1e-3                                  # ppb -> ppm
    dB = dipole_field(dchi, vox, bdir)                        # ppm
    os.makedirs(dst, exist_ok=True)
    rec = {"code": code, "ses": ses, "status": "ok", "n_sup": int(sup.sum()), "bdir": [round(float(b), 4) for b in bdir],
           "dB_max_ppb": float(np.abs(dB).max() * 1e3), "dB_mean_in_roi_ppb": float(dB[sup].mean() * 1e3)}
    for e in (1, 2, 3, 4):
        for part in ("mag",):
            for ext in ("nii.gz", "json"):
                f = f"{lab}_echo-{e}_part-{part}_MEGRE.{ext}"
                os.symlink(os.path.relpath(os.path.join(src, f), dst), os.path.join(dst, f))
        fj = f"{lab}_echo-{e}_part-phase_MEGRE.json"
        os.symlink(os.path.relpath(os.path.join(src, fj), dst), os.path.join(dst, fj))
        js = json.load(open(os.path.join(src, fj)))
        te, b0 = float(js["EchoTime"]), float(js["MagneticFieldStrength"])
        f0 = b0 * 42.576e6
        pimg = nib.load(os.path.join(src, f"{lab}_echo-{e}_part-phase_MEGRE.nii.gz"))
        ph = np.asarray(pimg.dataobj, dtype=np.float64)
        dphi = sign * 2.0 * np.pi * f0 * dB * 1e-6 * te
        new = np.angle(np.exp(1j * (ph + dphi))).astype(np.float32)
        hdr = pimg.header.copy()
        hdr.set_data_dtype(np.float32)
        hdr.set_slope_inter(1.0, 0.0)
        nib.save(nib.Nifti1Image(new, pimg.affine, hdr), os.path.join(dst, f"{lab}_echo-{e}_part-phase_MEGRE.nii.gz"))
        rec[f"dphi_max_e{e}"] = float(np.abs(dphi).max())
    return rec


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--level-ppb", type=float, required=True)
    a.add_argument("--sign", type=int, default=1, choices=[1, -1])
    a.add_argument("--out-bids", required=True)
    a.add_argument("--sessions", default=os.path.join(LX, "blind", "processing_order.txt"))
    a.add_argument("--limit", type=int, default=0)
    a.add_argument("--jobs", type=int, default=8)
    args = a.parse_args()
    out = os.path.abspath(args.out_bids)
    if not out.startswith(os.path.join(LX, "posthoc") + os.sep):
        sys.exit("refusing to write outside $LX/posthoc")
    if os.path.exists(out):
        sys.exit(f"{out} exists")
    src_bids = os.path.join(LX, "blind", "bids_all")
    os.makedirs(out)
    import shutil
    shutil.copy(os.path.join(src_bids, "dataset_description.json"), out)
    os.symlink(os.path.relpath(os.path.join(src_bids, "derivatives"), out), os.path.join(out, "derivatives"))
    cfg = L.load_config()
    sess = [(l.split()[0], int(l.split()[1])) for l in open(args.sessions) if l.strip()]
    if args.limit:
        sess = sess[: args.limit]
    jobs = [(c, s, args.level_ppb, args.sign, out, src_bids, cfg) for c, s in sess]
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        recs = list(ex.map(one, jobs))
    json.dump({"level_ppb": args.level_ppb, "sign": args.sign, "n": len(recs), "sessions": recs},
              open(out + "_injection_log.json", "w"), indent=1)
    print(f"{out}: {sum(r['status'] == 'ok' for r in recs)}/{len(recs)} sessions injected")


if __name__ == "__main__":
    main()
