"""LEAD-chi POST HOC (2026-10-07): I/O around the QSM.jl second engine.

  prepare  : sessions.tsv (code, ses, lab, B0, TE1-4, B0 direction in voxel axes, source dir, Mask 3 path) for every
             coded session in blind/processing_order.txt that has a Mask 3
  finalize : raw Julia outputs -> NIfTI in the QSMxT-like layout read by the registered scripts:
             $LX/qsm/<variant>/derivatives/qsmxt/sub-C/ses-K/anat/<lab>_{Chimap,mask,localfield}.nii
             (geometry copied from the echo-1 phase image; chi and local field in ppm)
"""
import argparse, json, os, sys

LX = os.environ["LX"]
sys.path.insert(0, os.environ.get("LX_SCRIPTS", os.path.join(LX, "scripts", "pipeline")))
import lx_common as L  # noqa: E402


def prepare(args):
    import numpy as np, nibabel as nib
    bids = args.bids
    mask_bids = os.path.join(LX, "blind", "bids_all")
    rows = []
    for line in open(os.path.join(LX, "blind", "processing_order.txt")):
        if not line.strip():
            continue
        c, s = line.split()[0], int(line.split()[1])
        lab = L.sess_label(c, s)
        src = os.path.join(bids, f"sub-{c}", f"ses-{s}", "anat")
        mp = os.path.join(mask_bids, "derivatives", "mask3", f"sub-{c}", f"ses-{s}", "anat", f"{lab}_mask.nii.gz")
        if not os.path.exists(mp) or not os.path.exists(src):
            continue
        img = nib.load(os.path.join(src, f"{lab}_echo-1_part-phase_MEGRE.nii.gz"))
        aff = img.affine
        vox = np.sqrt((aff[:3, :3] ** 2).sum(axis=0))
        b = np.linalg.solve(aff[:3, :3] / vox, np.array([0.0, 0.0, 1.0]))
        b = b / np.linalg.norm(b)
        tes, b0 = [], None
        for e in (1, 2, 3, 4):
            js = json.load(open(os.path.join(src, f"{lab}_echo-{e}_part-phase_MEGRE.json")))
            tes.append(float(js["EchoTime"])); b0 = float(js["MagneticFieldStrength"])
        rows.append([c, str(s), lab, f"{b0:.4f}"] + [f"{t:.6f}" for t in tes] + [f"{x:.6f}" for x in b] + [src, mp])
    if args.limit:
        rows = rows[: args.limit]
    with open(args.tsv, "w") as f:
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"{args.tsv}: {len(rows)} sessions")


def finalize(args):
    import numpy as np, nibabel as nib
    bids = os.path.join(LX, "blind", "bids_all")
    n = 0
    for v, chi_file in ((f"{args.prefix}_rts", "chi_rts.f32"), (f"{args.prefix}_tv", "chi_tv.f32")):
        base = os.path.join(LX, "qsm", v, "derivatives", "qsmxt")
        if os.path.exists(base) and not args.resume:
            sys.exit(f"{base} exists")
        os.makedirs(base, exist_ok=True)
        json.dump({"Name": f"LEAD-chi post hoc engine 2 ({v}; QSM.jl 0.5.4)", "BIDSVersion": "1.8.0",
                   "DatasetType": "derivative"}, open(os.path.join(LX, "qsm", v, "derivatives", "qsmxt",
                                                                    "dataset_description.json"), "w"))
        for line in open(args.tsv):
            f = line.rstrip("\n").split("\t")
            c, s, lab = f[0], int(f[1]), f[2]
            w = os.path.join(args.work, lab)
            if not os.path.exists(os.path.join(w, "done")):
                continue
            nx, ny, nz = [int(x) for x in open(os.path.join(w, "meta.txt")).read().split()[:3]]
            ref = nib.load(os.path.join(bids, f"sub-{c}", f"ses-{s}", "anat", f"{lab}_echo-1_part-phase_MEGRE.nii.gz"))
            rd = lambda fn, dt: np.fromfile(os.path.join(w, fn), dtype=dt).reshape((nx, ny, nz), order="F")
            anat = os.path.join(base, f"sub-{c}", f"ses-{s}", "anat")
            os.makedirs(anat, exist_ok=True)
            hdr = ref.header.copy(); hdr.set_data_dtype(np.float32); hdr.set_slope_inter(1.0, 0.0)
            nib.save(nib.Nifti1Image(rd(chi_file, np.float32), ref.affine, hdr), os.path.join(anat, f"{lab}_Chimap.nii"))
            nib.save(nib.Nifti1Image(rd("lf.f32", np.float32), ref.affine, hdr), os.path.join(anat, f"{lab}_localfield.nii"))
            mh = ref.header.copy(); mh.set_data_dtype(np.uint8); mh.set_slope_inter(1.0, 0.0)
            nib.save(nib.Nifti1Image(rd("mask.u8", np.uint8), ref.affine, mh), os.path.join(anat, f"{lab}_mask.nii"))
            n += 1
    print(f"finalized {n} session-variants")


def main():
    a = argparse.ArgumentParser()
    a.add_argument("stage", choices=["prepare", "finalize"])
    a.add_argument("--tsv", required=True)
    a.add_argument("--work", required=True)
    a.add_argument("--limit", type=int, default=0)
    a.add_argument("--resume", action="store_true")
    a.add_argument("--bids", default=os.path.join(LX, "blind", "bids_all"), help="input BIDS tree (prepare)")
    a.add_argument("--prefix", default="e2_lap_vsharp", help="output variant prefix (finalize): <prefix>_rts/_tv")
    args = a.parse_args()
    posthoc = "/path/to/Penn_LEAD/LEADchi/posthoc"     # real post hoc root (LX may be the shadow tree)
    for p in (args.tsv, args.work):
        if not os.path.realpath(p).startswith(posthoc + os.sep):
            sys.exit("refusing to read/write outside LEADchi/posthoc")
    if args.stage == "finalize" and not os.path.realpath(os.path.join(LX, "qsm")).startswith(posthoc + os.sep):
        sys.exit("finalize must run with LX pointing to the post hoc shadow tree (outputs stay inside posthoc/)")
    {"prepare": prepare, "finalize": finalize}[args.stage](args)


if __name__ == "__main__":
    main()
