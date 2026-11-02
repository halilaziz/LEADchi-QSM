#!/usr/bin/env python3
"""Unit tests (synthetic arrays only) for the registered reference and support rules (registration sec. 12.1, 12.2).

  python3 tests/test_reference.py
Checks: (1) the support S excludes the zero rim left by background removal inside Mask 3; (2) the two-step reference
ignores planted +/-2 ppm outliers and equals the tissue mean; the naive Mask-3 mean is biased; (3) the ROI in-support
share is measured against S; (4) the CSF-only reference is 'none' when |CSF - WB| > 50 ppb or voxels are too few, with
no PLIC fallback; (5) Otsu, slab-partition trimming, the ghost-to-signal index and the post-run mask assertion.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lx_common as L  # noqa: E402

ok = []
def check(name, cond):
    ok.append(bool(cond)); print(f"{name:66s} {'PASS' if cond else 'FAIL'}")

rng = np.random.default_rng(1)
shape = (40, 44, 36)
X, Y, Z = np.indices(shape)
m3 = (((X - 20) / 17) ** 2 + ((Y - 22) / 19) ** 2 + ((Z - 18) / 15) ** 2) <= 1
qmask = L.erode(m3, 1)                                   # reconstruction mask: one-voxel V-SHARP erosion
chi = np.where(qmask, 0.020 + rng.normal(0, 0.01, shape), 0.0)   # tissue mean 20 ppb; zero rim inside Mask 3
out = qmask & (rng.random(shape) < 0.01)
chi[out] = rng.choice([-2.0, 2.0], out.sum())          # planted +/-2 ppm failed-fit voxels
S = L.support_mask(m3.astype(float), chi, qmask.astype(float))
check("support S excludes the zero rim", not (S & ~qmask).any() and S.sum() > 0.9 * qmask.sum())
ref, imp, ref0 = L.two_step_reference(chi, S, 1.0)
tissue = chi[qmask & ~out].mean()
check("two-step reference equals the tissue mean (|diff| < 0.5 ppb)", abs(ref - tissue) * 1000 < 0.5)
check("every planted outlier is flagged implausible", imp[out].all() and not imp[qmask & ~out].any())
naive = chi[m3].mean()
check("the naive Mask-3 mean is biased by the rim and outliers (> 1 ppb)", abs(naive - tissue) * 1000 > 1.0)
roi = {"w": np.zeros(shape, np.float32), "sup": np.zeros(shape, bool)}
roi["sup"][18:23, 20:25, 16:21] = True; roi["w"][roi["sup"]] = 0.8
roi["sup"][2:4, 20:25, 16:21] = True; roi["w"][2:4, 20:25, 16:21] = 0.8     # part of the ROI in the rim / outside
st = L.roi_stats((chi - ref) * 1000, S, imp, roi)
exp_in = roi["w"][roi["sup"] & S].sum() / roi["w"][roi["sup"]].sum()
check("ROI in-support share is measured against S", abs(st["inmask"] - exp_in) < 1e-9 and st["inmask"] < 1)
cfg = L.load_config()
aseg = np.zeros(shape, np.int16); aseg[m3] = 2; aseg[14:26, 16:28, 12:24] = 4
chi_ppb = (chi - ref) * 1000
typ, mean, n = L.csf_reference(chi_ppb, S, imp, aseg, None, cfg)
check("CSF reference valid for a plausible ventricle", typ == "CSF" and n >= 200 and abs(mean) <= 50)
chi_bad = chi_ppb.copy(); chi_bad[aseg == 4] += 80.0
typ2, _, _ = L.csf_reference(chi_bad, S, imp, aseg, None, cfg)
check("CSF reference 'none' when |CSF - WB| > 50 ppb (no PLIC fallback)", typ2 == "none")
aseg_small = np.zeros(shape, np.int16); aseg_small[m3] = 2; aseg_small[19:22, 21:24, 17:20] = 4
typ3, _, n3 = L.csf_reference(chi_ppb, S, imp, aseg_small, None, cfg)
check("CSF reference 'none' when fewer than 200 voxels", typ3 == "none")
v = np.r_[rng.normal(0, 1, 5000), rng.normal(10, 1, 5000)]
check("Otsu threshold separates two modes", 3 < L.otsu_threshold(v) < 7)
mk = np.ones((20, 24, 96), bool)
t = L.trim_partition_edges(mk, 96, 2)
check("slab trimming removes 2 partitions at each end", t.sum() == 20 * 24 * 92 and not t[:, :, :2].any() and not t[:, :, -2:].any())
mag = np.zeros((40, 48, 30)); mag[10:30, 14:34, 5:25] = 100.0
g0 = L.ghost_to_signal(mag + rng.normal(0, 1, mag.shape), 0)
ghost = mag.copy(); ghost[np.roll(mag > 0, 20, axis=0) & ~(mag > 0)] = 10.0
g1 = L.ghost_to_signal(ghost + rng.normal(0, 1, mag.shape), 0)
check("ghost-to-signal index rises with a half-FOV ghost", g1 > g0 + 0.05)
ma = L.mask_assertion(qmask.astype(float), m3.astype(float))
mb = L.mask_assertion(np.ones(shape), m3.astype(float))
check("mask assertion passes for Mask 3 minus a rim and fails for a foreign mask", ma["mask_ok"] and not mb["mask_ok"])
print(f"\n{sum(ok)} of {len(ok)} reference/support tests passed")
sys.exit(0 if all(ok) else 1)
