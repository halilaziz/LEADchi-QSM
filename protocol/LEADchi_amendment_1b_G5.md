# LEAD-χ amendment 1b: G5 after the menu revision (decided by the user on 2026-10-03)

**Frozen before** the full QSM run and before any group information. The only χ that exists is still the 10-session pilot. The evidence is label-free and aggregate-only: `gates/g5_cmp_v1v2.json`, `gates/g5_diag3_v2.json`, `flags_coded/g5_summary.json`.

## Outcome of the amendment-1a menu revision (v2a = N4 on T1w + mask dilation 4; MI kept)

- **G5-H1 fail share among non-IQ-fail sessions: v1 15% (24/160) → v2a 27% (41/153).** IQ fails: 20 → 27.
- **MEGRE↔T1w Dice is unchanged by the revision.** Median change is 0 for brainstem, VDC and cerebellar WM; at most 1 session changes by more than 0.01. The MEGRE-side failures are therefore not registration failures. They track echo-1 magnitude image quality and SynthSeg agreement on it, which no registration change can fix.
- **The worsening comes entirely from the T1w→MNI step.** Brainstem + VDC Dice median 0.898 → 0.861; share < 0.80 rose from 2.8% to 21.1%; 105 sessions worse and 21 better. The dilated 4-voxel mask also feeds the T1w→MNI SyN, so the revision degraded a step it was not aimed at.
- Sessions passing G5-H1: v1 136/180 versus v2a 112/180. 24 pass in v1 and fail in v2a; none the other way.

## Changes

1. **Revert to the registered v1 chain.**
   - The active `register_rois.sh` is the registered version (sha256 `ebe8b6d1…`, kept in `code/pipeline_deviations/register_rois_v1.0_registered.sh`).
   - The v1 outputs (`reg_v1`) become `reg` again; the v2a outputs are kept as `*_v2a` for the report.
2. **The chain trigger is removed.**
   - G5 operates as a **per-session exclusion only**, using the 1a definitions: G5-H1 for first-usable criterion (e); G5-DN for dentate quantities; IQ fail excluded.
   - The "> 10% → methods report" rule is removed. Reason: the remaining failures are shown to be image-quality-driven (MEGRE Dice invariant to the registration revision, IQ index), not a pipeline failure, so the rule's premise does not hold.
   - This is declared as a deviation from the registered path, and the manuscript reports the full sequence: v1 → 1a revision → 1b revert.
3. **Bias safeguards** (the registered Dice-by-group comparison is kept and extended):
   - At unlock, G5-H1 pass/fail and IQ fail are compared by PS × ADHD cell (Fisher), and Dice by group (Wilcoxon).
   - If any comparison gives p < .05, H1 is co-reported with inverse-probability-of-inclusion weights from a label-free model of G5 failure (age, sex, mean FD of the same-session DWI, motion index), as a sensitivity analysis.
   - The registered "automated exclusions only" sensitivity analysis is kept.
4. **The G5 VDC threshold** is recalibrated on the v1 geometry under the unchanged pilot rule (expected 0.80, as in v1).

## Unchanged

Everything else in amendment 1a (IQ rule, G5-H1 / G5-DN split, dentate NA rule) and in the registration.
