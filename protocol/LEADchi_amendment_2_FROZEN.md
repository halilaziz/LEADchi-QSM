# LEAD-χ amendment 2 (FROZEN): pre-unlock freeze

Frozen 2026-10-05, **before the unlock**. No group label, phenotype or code↔ID key has been opened. All numbers below are label-free aggregates produced by the hashed scripts through `audit_run.py`.

**Freeze manifest:** `freeze/freeze_manifest.json`, 158 files, SHA-256 `efc4b909a394e94a6aed5d4c23fa18506a91ddb560bab31e08b8dc6a7614445f` (2026-10-05T07:01:52Z). It covers the coded tables, flags, QC sheets and κ, usable lists, setup/lock manifests and ciphertexts, frozen configs, the gate JSON, the analysis code and the registration documents.

## 1. Gate outcomes and branches (Q24 order)

| Gate | Outcome | Detail |
|---|---|---|
| G2 / G2-DN | **pass / pass** | within-QSMxT branch (A2 declined); SN–VTA ICC(3,1) 0.889 [0.78, 0.95]; dentate 0.771 [0.50, 0.93] |
| Engine / pipeline | QSMxT / primary (V-SHARP + MEDI, λ 0.002) | no fallback needed |
| PC1 | **pass** | GP median highest; GP > thalamus in 100% (n 107) |
| PC2w | **INDETERMINATE** | returner GP/putamen change − thalamus change: mean +2.6 ppb, one-sided p = .16, n = 26 pairs. Not a fail → no pipeline switch; PC2x (age → putamen/GP) is evaluated after unlock as registered |
| PC3 | **pass** | r(χ, R2\*) = 0.68 (n 80) |
| G3 | **pass → ROI = SN–VTA** | ICC(3,1) 0.82 [0.77, 0.87] over the 9 multiverse nodes (n 80) |
| PC4 | **pass** | split-echo r 0.924, R_SB 0.961, **r̂_L = 0.938** → outcome-B bound δ_B = −0.6·√0.938 = **−0.581** |
| PC4-DN | pass | r̂_L 0.705 |
| G4 | **pass** | linear κ 0.664 (36 re-rated; amendment 1c) |
| G5 | **per-session only** (amendment 1b) | fail share all 0.244, non-IQ 0.117; IQ fail 26; the registered trigger would have read "fail"; the chain trigger is removed (1b) |
| PC5 | **fail** | S3 exploratory (with n_min = none, #10d) |
| n_min | none | S3 exploratory |

## 2. Samples (first-usable rule; primary branch, SN–VTA)

- Usable sessions 107/181; **H1 first-usable subjects n = 80**; returner pairs with two usable sessions 26.
- Fail reasons (session level, overlapping): visual grade 2 → 62; G5 → 45; no/failed T1w → 6; auto grade 2 → 7; SN–VTA plausibility → 8; placement → 1; not registered → 1.
- The registration expected ≈ 100–106. The shortfall comes from the visual grade (33% grade 2, AI rater) and the image-quality G5 failures. The registered **automated-exclusions-only** sensitivity (no visual grade/placement) is co-reported.
- Realised precision: with n = 80 the H1 SE is ≈ 0.24σ (registered 0.206σ at n 112); power ≈ 0.70 at δ = 0.6. Stated in the manuscript.

## 3. Pooled label-free constants (first-usable, SN–VTA)

mean χ_SNVTA 53.0 ppb, SD 11.88 (= sd_total); mean χ_STR 16.6 ppb, SD 5.92; ρ(SN–VTA, STR) 0.19; SD(Δref) 6.1 ppb (0.51 of SD) → the K3 reference-drift rule is live; cross-session r 0.66 / ICC 0.66 (n 26; T1 informative); boundary distance median 0.81 mm (p95 0.97).

## 4. QC (amendment 1c)

Calibration κ 0.83 (3 instances × 15). Pass 1: 180 sessions, grades 0/1/2 = 7/114/59; placement 180 pass. Pass 2 (36, fresh instances): κ 0.664. Final (conservative rule): 7/112/61; dental D flagged 13. Both readings of "at the midbrain or dentate level" (structure vs slice) occurred across instances and are reported.

## 5. Deviations since amendment 1 (all logged in `amendment1_log.md` #28–46)

1a (G5 menu revision), 1b (revert to v1; G5 per session; trigger removed), 1c (AI rater; 14-day rule waived), sheet-validation and G5-gate implementation patches, DWI bundle QC deferred (no pyAFQ PNGs; S3 exploratory).

## 6. Amendment-2 modules

Tier 3–4 modules (E-dev, E-z, X1–X7, node-level max-T, DWI multiverse) are **not** part of this freeze; they are run after unlock as exploratory, labelled as such. D10 (ΔR2\*) not run.

## 7. What happens at unlock

`unlock_and_merge.py` (Mac, user, private key) joins the frozen coded tables to the locked labels once; the merged tables are then analysed by the frozen `run_all.R`. Any change after this point is a post-unlock deviation (Annex C.4).
