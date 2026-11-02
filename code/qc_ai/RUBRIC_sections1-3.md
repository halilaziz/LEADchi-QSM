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

