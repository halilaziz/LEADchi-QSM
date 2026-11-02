# LEAD-χ Annex C.4: R2* ROI voxel set, hashed code vs registered text (2026-10-07)

- **Registered text** (prereg §12.2, line 253): R2\*_ROI is the mean over the *same voxel sets as the matching χ ROI*: the support S, implausible voxels excluded and the same weights, restricted to Mask 3.
- **Hashed code** (`extract_rois.py`): it averaged R2\* over Mask 3 ∩ finite R2\* inside the ROI, without S and without the implausible-voxel exclusion. This was found by the post hoc local analyses (`../validation_local/`, declaration 7).
- **Both versions** were recomputed from the primary maps by `posthoc/scripts/r2s_textconform.py` (sha256 `7c03eb03…`; read-only, using the frozen `lx_common` functions). The output is `r2s_roi_code_vs_text.csv`. The code version reproduces the frozen table to 5e-9.

| Quantity (first-usable n = 80) | Hashed code | Registered text |
|---|---|---|
| Mean difference text − code, SN–VTA R2\* | — | +0.0004 s⁻¹ |
| K2: β_PS(R2\*) [95% CI] | +1.098 [0.080, 2.117] | +1.098 [0.078, 2.118] |
| K2: one-sided p (lower) | .9825 | .9824 |
| K2 met? | no | no |
| PC3: r(χ_GP, R2\*_GP) (criterion ≥ 0.30) | 0.679 | 0.679 |

No decision changes; both results are shown, as Annex C.4 requires. The difference is numerically negligible.

Code: `code/k2_pc3_textconform.R`. Numbers: `results.json`.
