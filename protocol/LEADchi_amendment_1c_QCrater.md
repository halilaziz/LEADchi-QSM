# LEAD-χ amendment 1c: visual QC rater (decided by the user on 2026-10-04)

**Frozen before any QC grade exists** and before any group information. χ exists for the full sample (primary, split-echo and multiverse runs), but no ROI value, no group label and no code key has been opened by anyone.

## Registered design (qc_rubric.md §1, registration Q13)

Single rater: the first author grades every session from coded montages; a random 20% is re-rated ≥ 14 days later; G4 = linear-weighted κ ≥ 0.60.

## Change: the rater is an AI assistant (Claude, model Fable 5.1), not the first author

**Reason (user's decision).** The first author's prior viewing of raw MEGRE images is declared "unsure" (Q18). An AI rater instantiated in this session has seen no raw ds007116 MEGRE image, no χ value, no ROI table, no group label and no phenotype; it sees only the coded montages (random QC IDs), the rubric and the forced-review flag. The user also judged the AI better placed to apply the written rubric consistently across 180 sessions.

**Procedure.**
1. **Material and blinding** are unchanged: `qc/montages/{gross,placement}/<qc_id>.png`, rating sheets by QC ID, no χ numbers, no colour-bar values, random order (seed 20260924 + pass).
2. **Rater instances.** Each montage is graded by an independent AI agent instance with a fresh context. An instance receives: the rubric text verbatim, the two montages of its sessions, the forced-review flag, and nothing else. It is forbidden to read any other LEADchi file. It outputs, per session, `grade` (0/1/2), `failure_codes`, `placement` (pass/fail) and a one-line `notes`.
3. **Calibration (added).** Before pass 1, 15 randomly chosen sessions are graded by 3 independent instances. Linear-weighted κ across instances is reported. If κ < 0.60, the rubric prompt is clarified once (wording only; the rubric's criteria are unchanged), logged with its hash, and calibration is repeated once.
4. **Pass 1.** All 180 registered sessions, split across instances in the sheet's random order; each instance grades ≤ 15 sessions.
5. **Pass 2 (20% re-rating).** The registered 14-day interval exists to break a human rater's memory of pass 1. For AI instances, memory independence is obtained by construction: pass 2 is graded by **fresh instances in a separate run that has no access to the pass-1 sheet**. The 14-day wait is therefore waived; the script's 14-day check is bypassed with an explicit `--independent-rater` flag that writes this fact into the sheet header. This is a deviation from Annex A.5 and is declared as such.
6. **κ and G4** are computed by the registered `qc_montage.py kappa` (linear-weighted κ of the 0/1/2 grade; G4 ≥ 0.60), unchanged. The stop rule (at most one rubric refinement) is unchanged.
7. **Finalization** by `qc_montage.py finalize`, unchanged (conservative rule for re-rated items).
8. **Declared limitation.** An AI's visual grading has not been validated against expert radiological QC on this sequence. Mitigations: (a) the calibration κ (item 3); (b) the registered **automated-exclusions-only sensitivity analysis** (plausibility, G5, T1w QC; no visual grade), which is co-reported with H1 and S1–S2 in all cases; (c) the registered QC-by-cell test at unlock (Fisher); (d) the user may audit any montage after unlock without changing frozen grades; any disagreement is reported, not applied.
9. **Automatic grade 2** (f_imp > 1% or ≥ 2 iron-rich ROI failures, from `qsm_plausibility.py`) is applied as registered and cannot be overridden by the rater.

## Unchanged

Everything else in the registration and in amendments 0, 1, 1a and 1b.
