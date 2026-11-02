# LEAD-χ POST HOC: synthetic SN–VTA recovery test (2026-10-07)

**Status: post hoc.** The registration did not include this analysis. It answers a methods objection raised after unlock: "Can the pipeline see an SN–VTA deficit at all?" It changes no registered result (outcome B stands; amendment1_log #57).

## What was done
1. **Injection.**
   - Region: the SN–VTA ROI voxel set that `extract_rois.py` actually averages (CIT168 SNc+SNr+VTA+PBP, summed probability ≥ 0.5, outside the 1-voxel CSF buffer).
   - Change: a uniform susceptibility change L, at L = −5, −10, −14.24, −20 and +10 ppb, in all 180 registered coded sessions.
   - Field: k-space dipole convolution, zero-padded 2×, B0 direction from the affine.
   - Phase: the field was converted to phase per echo (Δφ = 2π·f0·ΔB·TE) and added to the measured phase. Magnitudes, sidecars and Mask 3 were linked unchanged.
   - Script: `posthoc/recovery/scripts/inject_sn.py`, sha256 `c6c389fd…`.
2. **Re-run.**
   - QSMxT v9.22.0 with the registered primary settings: ROMEO, V-SHARP 12 → 1.5 mm, MEDI λ = 0.002, Mask 3, the same TOML (`run_rec.sh`, `5216e7a3…`).
   - Then the unchanged registered `qsm_plausibility.py` and `extract_rois.py`.
3. **Analysis** (`code/recovery_analysis.R`, using the frozen helper functions unchanged).
   - Recovery factor per session: k = (χ_SNVTA,rec − χ_SNVTA,primary)/L, computed on the raw ROI means.
   - Linearity across L.
   - Whether k depends on group or quality metrics.
   - A counterfactual: the recovered change was added to the PS youths only, and M1 and the decision rule were refitted.

## Results (first-usable n = 80 unless stated)

**Recovery is about half, linear on average, with no detectable group difference.**

| L (ppb) | k mean | SD | median | min to max |
|---|---|---|---|---|
| −5 | 0.470 | 0.088 | 0.474 | 0.16 to 0.96 |
| −10 | 0.467 | 0.059 | 0.472 | 0.30 to 0.74 |
| −14.24 | 0.473 | 0.052 | 0.473 | 0.37 to 0.68 |
| −20 | 0.469 | 0.047 | 0.470 | 0.37 to 0.63 |
| +10 | 0.477 | 0.058 | 0.480 | 0.33 to 0.70 |

- **All sessions.** Across all 180 sessions and five levels, the change regresses on L with slope 0.474 and intercept −0.05 ppb (R² 0.96, 900 session-levels). On average the response is linear and sign-symmetric. Per session it is noisier: k ranges 0.07–1.72 at −5 ppb, and one non-analysed session (`sub-20494_ses-2`) has k = −0.19 at +10, i.e. a wrong-sign change. An earlier test at +200 ppb in 3 sessions gave 0.50–0.57.
- **Group.** There is no detectable group difference in k at L = −14.24. The model is the M1 form without the ROI-volume and Q_phase terms (df 70): PS +0.0006 (95% CI −0.021 to 0.022; p .95), ADHD +0.003 (p .79). Cell means are 0.46 to 0.49. The CI limits a PS difference to about −4% to +5% of k, so a material differential sensitivity is unlikely for this synthetic, phase-only change.
- **QC metrics.** k shows no clear association with registration Dice (ρ .02), motion (−.09), ROI volume (−.01) or Q_phase (.04); with n = 80, correlations up to about ±0.22 cannot be ruled out. With age, ρ = .26: p .02 uncorrected, about .10 after Holm over the 5 correlations (exploratory).
- **Off-target effects are small on average.**
  - Whole-brain reference: mean shift 0.001–0.007 ppb across levels.
  - Thalamus: mean shift −0.17 ppb at L = −14.24, about −1.2% of L at every level.
  - GP: mean shift −0.27 to +0.57 ppb across levels, but the per-session changes are not negligible (SD 1.6–1.9 ppb, range about ±6 ppb at every level, including +10). Because this spread does not scale with L, it points to run-to-run variability of the re-reconstruction.

**Post hoc descriptive interpretation for outcome B.** If a real deficit were recovered like the synthetic one, a true change of D ppb in PS youths would show up as about 0.47·D. Treating k as known, the observed β_PS = +4.96 ppb (SE 3.125, df 68) is incompatible with such changes at these one-sided p values:

| True deficit D (ppb) | Expected observed shift | One-sided p for H0: β_PS = 0.473·D against β_PS > 0.473·D |
|---|---|---|
| −5 | −2.36 | .011 |
| −10 | −4.73 | .0014 |
| −14.24 (numerically equal to the adult raw difference) | −6.73 | .0002 |
| −20 | −9.46 | 9e-6 |

- **True-scale translation** (k = 0.473, treated as exact; over the 95% CI of mean k, 0.461–0.484, the bound becomes −12.7 to −13.3 ppb). This descriptive post hoc translation divides the registered quantities by k; it is not a registered result:
  - the registered B bound (−6.13 ppb, observed scale) corresponds to a true deficit of −13.0 ppb;
  - the observed one-sided 95% lower limit (−0.25 ppb) corresponds to −0.53 ppb true;
  - the two-sided 95% lower limit (−1.27 ppb) corresponds to −2.7 ppb true.
- **Counterfactual** (recovered change added to the PS youths of this sample, frozen M1 and rule refitted):

  | True deficit | β_PS | Outcome |
  |---|---|---|
  | −5 | +2.61 | B |
  | −10 | +0.22 | B |
  | −14.24 | −1.80 | C |
  | −20 | −4.46 | C |
  | +10 | +9.78 | B, and "disconfirming" by the registered rule (p .003) |

  At the observed SE and σ_B, B holds for any β_PS above about −0.92 ppb (−6.132 + 1.668 × 3.125). Added to what was observed, an extra true deficit of 10 ppb would still have given B, whereas 14.24 ppb or more would have given C. The pipeline therefore recovers about half of a uniform SN–VTA change, and that is the sensitivity against which B should be read.

## Caveats
- **The attenuation is not specific to MEDI** (corrected 2026-10-07 after the engine-2 run; an earlier draft of this README guessed it was). Engine 2 (QSM.jl V-SHARP + RTS or TV, which use no magnitude prior) recovers even less: k = 0.36 (RTS) and 0.25 (TV) against 0.47 for MEDI (`../engine2/README.md`). Possible shared causes (not tested) are SMV/V-SHARP high-pass filtering and regularised inversion acting on a small structure, 3–4 voxels thick, at 1.5 mm. Whether a real iron change, which also alters magnitude/R2\*, would be recovered better is untested.
- **Partial volume is not modelled.** The change is uniform inside the ROI voxel set, and no partial-volume effect enters.
- **Other pipelines.** The adult reference (−14.24 ppb) was measured with another pipeline that has its own, unknown recovery factor. The true-scale numbers above apply to this pipeline only.
- **Not a calibration.** k is not used to change any registered result; the true-scale translation above is a descriptive post hoc aid for translating ppb across pipelines.
- **Injection geometry.** The B0 direction was computed from each affine; it is (0, 0, 1) in every session, consistent with the non-oblique acquisitions (header audit: obliquity 0°).

## Files
- `raw/roi_qsmxt_{primary,rec_m05,rec_m10,rec_m14,rec_m20,rec_p10}.csv`: coded ROI tables copied from server, with `SHA256SUMS`.
- `raw/injection_log_rec_m14.json`: the per-session injection log.
- `recovery_long.csv`: one row per session and level.
- `results.json`: all the numbers above.
- `code/recovery_analysis.R`: the analysis script.
- Server side: `posthoc/recovery/` (scripts, injected BIDS trees, logs); `qsm/rec_*` (reconstructions).
