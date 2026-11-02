# LEAD-χ group-blind visual QC rubric and κ procedure

Task P5b, 2026-09-25; revised 2026-09-26 for the final registration (code D, the single-rater design, the stop rule,
the DWI bundle-QC rubric). Registered with the pipeline code (SHA-256 in `SHA256SUMS`). It implements:
- protocol v1.1 §11.4 and §6 criteria (c) and (d);
- v1.2 delta §2.6 and §5.4;
- gate G4.

The material is produced by `qc_montage.py`, and the frozen grades feed `build_usable.py`.

## 1. What the rater sees and does not see

| Sees | Never sees before unlock |
|---|---|
| Montages named by a random QC ID (`q` + 8 hex), in random order | Group, cell, diagnosis, age, sex, any phenotype |
| The automated **forced-review** flag (plausibility or slab coverage) | The code↔QC-ID map (`qc/qc_ids.csv`), subject codes and within-subject session order |
| Magnitude, field and χ images with fixed display windows | Any χ, R2*, ROI, Dice, orientation or coverage **number**; colour bars carry no values |
| This rubric | Raw (uncoded) ds007116 images; round-1 figures or tables |

Rules:
- **Single-rater design (declared).** One rater, the first author, grades every session. He is also the person whose
  prior viewing of raw MEGRE images is "unsure" (registration Q18). G4 therefore measures **intra-rater** consistency
  only; no inter-rater agreement is claimed.
- An AI assistant may prepare the sheets but does not grade.
- Grades are entered only in `qc/rating_pass<k>.csv`.
- A sheet is never edited after it is closed.
- An automated-exclusions-only sensitivity (plausibility, G5 and T1w QC; no visual grade or placement rating) is
  registered, so that the visual grades are never the only route to a primary-sample decision.

## 2. Gross-failure grade (montage `gross`; criterion (c))

Each gross montage shows three axial levels: dentate, midbrain (SN–VTA) and basal ganglia (GP). The panels are:
- echo-1 and echo-4 magnitude;
- total field and local field (when QSMxT writes them);
- χ, shown in greyscale over −0.2 to +0.2 ppm, with the Mask-3 outline in yellow.

χ is used **only** to find gross failures, never to judge anatomy or ROI placement.

| Grade | Meaning | Examples |
|---|---|---|
| **0 pass** | No visible artefact at the three levels | — |
| **1 minor** | An artefact is visible but spares the midbrain and dentate levels, or is faint | Mild ghosting; a small residual wrap in cortex; faint streaks outside the posterior fossa; minor air artefact at the frontal sinus |
| **2 fail** | Any failure code below, at the midbrain or dentate level or across the whole slice | See the failure codes |

**Failure codes.** Record them in `failure_codes`, separated by `;`.

| Code | Failure (grade 2 when present at the midbrain or dentate level) |
|---|---|
| M | Motion ghosting or blurring across the brainstem or cerebellum on magnitude |
| W | Residual phase wraps or unwrapping discontinuities in the total or local field inside Mask 3 |
| S | Streaking on χ that crosses the whole midbrain or dentate level |
| B | Background-field residue: a smooth large-scale gradient in the local field or χ, especially at the posterior fossa or skull base |
| K | Mask 3 excludes the midbrain or a dentate, or cuts through the brainstem |
| F | R–L phase wrap-around (75% phase FOV) that overlaps the brain |
| A | Air or susceptibility artefact inside the SN–VTA or dentate |
| C | Slab coverage: the midbrain or a dentate is not fully in the field of view |
| D | Dental or metal artefact (braces, retainers, implants): a signal void or distortion anywhere in the slab. **Record D for every session in which it is visible, whatever the grade.** It forces grade 2 when the void or distortion reaches the midbrain or posterior-fossa level, or when the residual-background index flags the session |

Two further grading rules:
- **Forced review.** When the flag is 1, the rater must inspect every panel at full size before grading. The flag itself does not change the grade.
- **Automatic grade 2.** `qsm_plausibility.py` assigns grade 2 when f_imp > 1% or when ≥ 2 iron-rich ROIs fail. The rater cannot override it.

## 3. ROI placement (montage `placement`; criterion (d))

Placement is judged on **T1w (resampled to the MEGRE grid) and echo-1 magnitude only**, never on χ.

The reason: lower SN χ, the hypothesised effect, makes the SN less conspicuous on χ, so χ-based judgements could differ by group.

The montage shows the SN–VTA outline (CIT168, w ≥ 0.5) in axial, coronal and sagittal planes, plus the dentate outlines.

**Pass** requires all of the following:
1. The SN–VTA outline lies in the ventral midbrain, between the cerebral peduncle and the red-nucleus/tegmentum region.
2. It does not extend into the interpeduncular cistern or the third ventricle.
3. It is left–right symmetric within about one voxel.
4. The dentate outlines lie in the deep cerebellar white matter lateral to the fourth ventricle.

**Fail** if any of the following holds:
- the outline is displaced by more than about 2 voxels (3 mm) from these landmarks;
- it is visibly misregistered, i.e. the brainstem outline on T1w and magnitude disagree;
- it overlaps CSF over more than about a third of its extent in any plane.

Write `pass` or `fail`, and give the reason in `notes`.

## 4. Procedure and κ (gate G4)

1. `qc_montage.py montage --kind gross` and `--kind placement`, then `sheet --pass 1`. Every coded session is graded. The order is random with seed 20260924 + pass.
2. When pass 1 is complete, run `finalize --close-pass1`, which time-stamps closure.
3. **At least 14 days later**, run `sheet --pass 2`. This produces a random 20% re-rating sheet. The script refuses to run earlier. The rater does not look at pass 1 while rating.
4. `kappa` computes:
   - the linear-weighted Cohen κ of the 0/1/2 grade, which is the decision statistic;
   - the quadratic-weighted κ and the unweighted placement κ, which are co-reported.
   It writes aggregate JSON only.
5. **G4 passes if linear-weighted κ ≥ 0.60.** If it fails (stop rule: at most ONE refinement):
   - the rubric wording is refined blind, once, and the change is logged with a new rubric hash;
   - every session is re-rated in pass 3 (`sheet --pass 3`, then `finalize --close-pass3`);
   - at least 14 days later a new random 20% re-rating of pass 3 is drawn (`sheet --pass 4`) and κ is recomputed
     (pass 3 vs pass 4);
   - if G4 fails again, the pass-3 grades are final, G4 is recorded as failed and every H1 statement carries the
     caveat "QC reliability not established". No further refinement is made.
6. `finalize` writes `qc/grades_final.csv` (code, ses_order, qc_grade, placement, re_rated, failure_codes, dental_D)
   and `qc/g4_final.json`.
   - The final grade is the reference-pass grade (pass 1, or pass 3 after the refinement).
   - For re-rated items, the worse of the two grades is used, placement fails if either rating failed, and the
     failure codes are the union of both ratings.
7. The QC tables and κ are hashed into the freeze manifest (`build_usable.py --stage manifest`) and posted in amendment 2 before unlock. After unlock, QC failure is tested by cell (Fisher's exact test; P5 `qc_by_group`).

## 5. Registered implementation rules

- Linear κ weights are the G4 decision statistic.
- The conservative rule is used for re-rated items.
- Raters work from random QC IDs rather than subject codes.
- Visual grades are shared across the fallback and SEPIA branches, which are graded on the primary montages. Their
  automated flags are branch-specific.
- χ display window: greyscale −0.2 to +0.2 ppm, no colour-bar values; magnitude panels windowed at their 1st–99th
  percentile, field panels symmetric at the 99th percentile of |field| inside Mask 3.

## 6. DWI bundle QC (S3 inclusion; registration Q13)

- **Material.** The pyAFQ per-session bundle figures (HTML/PNG) of SCP_R and MCP_R, copied under a separate random
  DWI code set; the DWI code↔ID file is plaintext in `LEADchi/dwiqc/`, hashed, and not opened until the DWI-QC freeze.
- **Rater.** The first author. He is **not blind in practice**: he has seen the Penn LEAD DWI group results and the
  round-1 20-subject AFQ-QC list (declared).
- **Pass rule (per bundle).** Pass unless any of: the bundle is absent or has fewer than 20 streamlines; it does not
  connect the dentate-hilum region to the midbrain (SCP_R) or the pons to the cerebellar hemisphere (MCP_R); more than
  about a third of its streamlines are visibly spurious (leaving the peduncle); the tract is displaced into CSF over
  more than a third of its length. S3 requires an SCP_R pass (and an MCP_R pass for NC2).
- **Freeze.** Grades are entered once, mapped back to real IDs at the DWI-QC freeze (`bundle_qc.csv`: participant_id,
  session_id, scp_r_pass, mcp_r_pass), hashed, and read by `dwi_geometry.py --stage table`.
