# LEAD-χ amendment 1a: G5 registration gate (decided by the user on 2026-10-02)

**Frozen before** the full QSM run, before any group information, and before the revised registration has been run. The only χ that exists is the 10-session pilot, which is the registered exception. The decision rests only on label-free, pre-named G5 metrics, shown as aggregates (`gates/g5_diag_v1.json`, `g5_diag2_v1.json`, `g5_diag3_v1.json`), and on coded, group-blind anatomical overlays of 3 failing and 1 passing session (echo-1 magnitude, warped T1w; no χ).

## Observation (registered chain v1, 180 registered sessions)

- **G5 fail share: 46%** (83/180). The MEGRE↔T1w side accounts for nearly all of it: 78 sessions fail on that side only, 4 on both sides, and 1 on the T1w→MNI side only.
- **Gross failures: 20 sessions (11%)**, with brainstem Dice < 0.5. The overlays show acceptable alignment but **motion-corrupted echo-1 magnitude** (ringing and ghosting). On these images SynthSeg finds ≈ 0 brainstem voxels, so they are image-quality failures, not registration failures. FOV coverage is ≈ 1.0 in all of them.
- **Cerebellar-WM Dice:** median 0.806, with 45% below 0.80. Among the non-gross sessions, 41 lie at 0.75–0.80. This is consistent with lower cerebellar WM/GM contrast in 1.5 mm GRE magnitude for SynthSeg.
- **Without the cerebellar-WM criterion,** 15% of the non-gross sessions still fail (marginal brainstem or VDC).
- **The registered trigger cannot be met.** Under the registered "> 10% → revise once → if > 10% still fail → methods report" rule, the image-quality failures alone (11%) would force the methods-report path whatever the registration does.

## Changes (Stage II technical change, Annex C.3: a pre-named label-free metric revealed the problem; applied to all sessions)

1. **Registered menu revision**, applied once, blind, to all sessions:
   - **N4 bias correction of the T1w** before brain masking and registration;
   - **brain-mask dilation 2 → 4 voxels**.
   
   The third menu item, the rigid metric MI → CC, is **not** chosen. MEGRE echo-1 magnitude and T1w differ in contrast, and MI is the appropriate cross-contrast metric. This choice was made before any revised result existed.
2. **Image-quality failure (new, label-free, registration-independent).** IQ fail = MEGRE-SynthSeg brainstem volume < 50% of the T1w-aseg brainstem volume, both on the MEGRE-segmentation grid and inside the MEGRE FOV.
   - IQ-fail sessions fail first-usable criterion (e) (excluded per session, as registered).
   - They are **not** counted toward the chain trigger.
3. **G5 split by analysis:**
   - **G5-H1** (used for first-usable criterion (e) and therefore H1, K1–K4, S1, S2 and T1):
     - T1w→MNI brainstem + VDC Dice ≥ 0.80;
     - MEGRE↔T1w brainstem Dice ≥ 0.80;
     - MEGRE↔T1w VDC Dice ≥ the VDC threshold (pilot-calibrated rule unchanged, re-applied to the revised chain).
   - **G5-DN** = G5-H1 **and** MEGRE↔T1w cerebellar-WM Dice ≥ 0.80. It is required for every dentate-based quantity (S3, X1, PC4-DN, G2-DN): dentate χ / R2\* columns are set to NA where G5-DN fails.
4. **Chain trigger after the revision.** If > 10% of the **non-IQ-fail** registered sessions fail G5-H1, G5 is recorded as failed and the registered methods-report path applies.
   - Pre-revision value for comparison: 15% (24/160) under the G5-H1 definition.
   - Both before and after aggregates are reported.
5. **Downstream of the registration.**
   - The registered steps are recomputed on the revised chain: `register_rois.sh --stage all`, `--stage harmref`, `extract_rois.py --stage geom`, `pilot_checks.py --stage g5cal` and `extract_rois.py --stage g5`.
   - Version-1 outputs are kept as `*_v1` for the before/after report.
   - The pilot sign check (pass; all four pooled medians > 0) and the λ rule (λ = 0.002) are not repeated: neither depends on the T1w chain beyond ROI placement, and the sign result has a wide margin.

## Unchanged

H1, K1–K4, S1, S2, S3 (exploratory per amendment 1), T1, all thresholds other than the G5 split above, multiplicity, power statements, and the decision rules.
