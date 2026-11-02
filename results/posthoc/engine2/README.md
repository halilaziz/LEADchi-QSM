# LEAD-χ POST HOC: engine-2 check with an independently coded QSM chain (2026-10-07)

**Status: post hoc.** The registration declined a second engine (A2), so this analysis is new (amendment1_log #57, #59). Outcome B is not re-decided. The engine-2 "B" entries below are a **B pattern, not a registered outcome**: the bound borrows δ_B from the primary engine's split-echo reliability, and engine 2's own precision gates were not assessed.

## What was done
- **Engine 2** is QSM.jl 0.5.4 (Kames), run under Julia 1.10.9 in `LEADchi/tools/julia*`. FFTW.jl is pinned to 1.8.0 because QSM.jl 0.5.4 fails to load with FFTW.jl 1.10. Its chain is QSM.jl's documented default:
  1. Laplacian phase unwrapping (multigrid-preconditioned conjugate gradient (MGPCG), Dirichlet on Mask 3);
  2. magnitude-weighted multi-echo linear fit with a phase offset, giving the field in ppm;
  3. V-SHARP (QSM.jl defaults: radii 27 → 3 mm, thr 0.05);
  4. dipole inversion: **RTS** (main; defaults) and **TV** (sensitivity; defaults).
- NLTV diverged to NaN on these data and was replaced by TV. RTS and TV share steps 1–3.
- **What is independent:** the code. No QSMxT/QSM.rs component is used. The method is only partly independent: both chains use V-SHARP-type background removal and share the coded inputs, Mask 3, the ROI warps, the plausibility rules and the two-step whole-brain reference.
- **How it was run:** the registered `qsm_plausibility.py` and `extract_rois.py` were run unchanged, with `LX` pointed to the shadow tree `posthoc/shadow/`, so nothing was written into the frozen folders. Engine 2 was also run on the −14.24 ppb injected tree to obtain its own recovery factor k.
- **Scripts** (on server): `posthoc/engine2/scripts/e2_qsmjl.jl` `4a013bf7…`, `e2_io.py` `30bcc527…`, `e2_chain.sh` `fb952f1f…`. Local analysis: `code/engine2_analysis.R`, which uses the frozen helper functions unchanged.
- **Sample:** the 80 primary first-usable sessions.
  - One session (`sub-20866_ses-1`) has no engine-2 value (run status `mask_mismatch`).
  - In 16 more sessions, engine 2's SN–VTA ROI fails the registered plausibility rule (in-support share 0.82–0.95, below 0.95). All 80 primary values pass (minimum 0.969).
  - **Main analysis (n = 63):** registered rule §12.4 applied, so a failing engine reconstruction is missing for that engine only.
  - **Sensitivity (n = 79):** the failures are kept.

## Results

**Agreement with the primary (SN–VTA, n = 79 with an engine-2 value)**

| | RTS | TV |
|---|---|---|
| ICC(3,1) | 0.68 | 0.69 |
| Pearson / Spearman | 0.73 / 0.71 | 0.73 / 0.71 |
| Mean χ, ppb (primary 53.0) | 27.1 | 24.0 |
| SD, ppb (primary 12.0) | 8.1 | 8.7 |
| GP Pearson with primary | 0.89 | 0.90 |

Engine 2 recovers about half of the primary's SN–VTA contrast (51% RTS, 45% TV). Across other deep-grey-matter ROIs the ratio is 44–97%. Youths are ranked moderately consistently (Spearman 0.71).

**H1: frozen M1 and decision rule** (engine-2 rows: main analysis, n = 63)

| | Primary QSMxT (n = 80) | Engine 2 RTS | Engine 2 TV |
|---|---|---|---|
| β_PS, ppb [95% CI] | +4.96 [−1.27, 11.20] | +2.54 [−2.94, 8.01] | +3.34 [−2.37, 9.05] |
| p (two-sided) | .117 | .36 | .25 |
| δ = β/σ_resid | +0.44 | +0.31 | +0.38 |
| Bound −δ_B·σ_B, ppb; p_B | −6.13; .0004 | −4.29; .008 | −4.71; .003 |
| Outcome | B (registered) | B pattern | B pattern |
| Huber β | +3.94 | +1.63 | +2.07 |
| PS × ADHD interaction, ppb [95% CI] | −6.39 [−18.8, 6.0], p .31 | +1.37 [−8.3, 11.0] | +0.67 [−9.6, 11.0] |
| PS-only − TD | +8.16 (p .044) | +1.85 (p .59) | +3.01 (p .40) |
| PS+ADHD − ADHD-only | +1.77 (p .71) | +3.22 (p .41) | +3.67 (p .38) |
| Per-stratum p_B (PS-only; comorbid), descriptive | .0003; .052 | .038; .029 | .017; .024 |

The sensitivity analysis with engine-2 plausibility failures kept (n = 79) gives:

| | β_PS, ppb [95% CI] | p | δ | p_B |
|---|---|---|---|---|
| RTS | +3.81 [−0.51, 8.14] | .083 | +0.50 | .0003 |
| TV | +4.37 [−0.25, 8.99] | .063 | +0.53 | .0002 |

Strata in the sensitivity analysis:
- RTS: +3.56 and +4.07;
- TV: +4.52 and +4.22;
- per-stratum p_B ≤ .007.

**Recovery factor k at −14.24 ppb (first-usable)**

| Chain | k mean (SD) |
|---|---|
| Primary (MEDI) | 0.47 (0.05) |
| Engine 2 RTS | 0.36 (0.04) |
| Engine 2 TV | 0.25 (0.06) |

Treating k and −14.24 ppb as fixed, a true deficit numerically equal to the adult raw difference would be expected to appear as −6.7 (primary), −5.2 (RTS) or −3.5 ppb (TV). Each engine's observed estimate is incompatible with it (post hoc, one-sided):

| Analysis | Primary | RTS | TV |
|---|---|---|---|
| Main (n = 63 for engine 2) | p = .0002 | p = .003 | p = .010 |
| Engine-2 sensitivity (n = 79) | — | p = 5e-5 | p = 5e-4 |

Under TV, the registered-style bound (−4.71 ppb) is less stringent than the adult deficit as TV would record it (−3.52 ppb). TV's B pattern therefore does not by itself exclude an adult-magnitude deficit; only the k-scaled test does. For the primary (−6.13 vs −6.73) and RTS (−4.29 vs −5.16), the bound is the stricter of the two.

## Reading (post hoc, descriptive)
- **Engine-2 check.** An independently coded chain gives a point estimate on the same (positive) side and the same B pattern. This is consistent with, but does not demonstrate, insensitivity of the registered result to the QSM engine:
  - only one alternative engine was tested;
  - its two inversions share most of the chain;
  - all CIs include 0, so no direction is established.
- **Size of the engine-2 estimates.** They are smaller once the registered plausibility rule is applied (δ +0.31/+0.38 vs +0.44).
- **ADHD strata.** Engine 2 does not reproduce the primary's contrast between the two ADHD strata: its stratum estimates are similar and its interaction is near zero. The primary's PS-only − TD contrast should therefore not be interpreted.
- **Comparing engines.** δ is more comparable across engines than ppb, but it still depends on each engine's recovery k and its noise. Relative to the primary, k is ×0.77 (RTS) and ×0.52 (TV), and σ_resid is ×0.72 and ×0.78 (main analysis). The k-scaled comparison above addresses this.
- **Agreement.** Between-engine agreement is moderate (ICC ≈ 0.68). This is a measurement limitation and is reported as such.

## Files
- `raw/`: coded ROI tables from the shadow tree and the plausibility summaries, with SHA256SUMS.
- `agreement_first_usable_{rts,tv}.csv`: per-session primary vs engine-2 values, run status and plausibility flag.
- `results.json`: every number above, including the primary's Huber estimate, bound, per-stratum p_B and k-scaled adult test.
- `run_log.txt`: the analysis run log.
- `code/engine2_analysis.R`: the analysis script.
- `verification/VERIFY.md`: independent Python re-derivation of the first version. All numbers matched; its description and wording corrections are applied here.
- Server side: `posthoc/engine2/` (scripts, logs, raw outputs) and `posthoc/shadow/` (variant NIfTI and tables).
