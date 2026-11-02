# LEAD-χ visual QC re-rating: instructions for the rater

Thank you for rating these images. Please read this page once before you start. It takes about 10 minutes.

## 1. Purpose

The quality of every brain scan in this study was graded by an **automated rater**, using a written rubric and two image montages per scan. We need an experienced human rater to grade a set of the same montages **independently**, with the same rubric, so that the automated grades can be checked.

You are blind to everything except the images:
- you do not know any participant's group, diagnosis, age or sex;
- you do not know the automated grade of any item.

Please keep it that way until you have returned your ratings (see section 6).

## 2. What you received

| File | What it is |
|---|---|
| `rate.html` | The rating tool. It is one offline file: it makes no network connection and sends nothing anywhere. |
| `images/gross/<kit_id>.png` | The gross-failure montage of each item |
| `images/placement/<kit_id>.png` | The ROI-placement montage of each item |
| `rating_sheet_blank.csv` | A paper-style fallback sheet, in the same order as the tool |
| `README_rater.md` | This page |

There are **111 items**. Each item is named by a random kit ID (`h` followed by 6 characters). The kit ID tells you nothing about the participant, the session or the group. The montage titles show only the kit ID.

## 3. What the montages show

**Gross montage** (`images/gross/`): a grid of 3 rows × 5 columns.
- **Rows** are three axial levels, from top to bottom: **dentate**, **midbrain** (SN–VTA level) and **basal ganglia** (GP level). The level name is written at the left of each row.
- **Columns**, from left to right:
  1. echo-1 magnitude ("mag e1");
  2. echo-4 magnitude ("mag e4");
  3. total field;
  4. local field;
  5. χ, in greyscale over −0.2 to +0.2 ppm ("chi (+/-0.2 ppm)"), with the Mask-3 brain-mask outline in **yellow**.
- The field panels are shown in grey with a symmetric window. No panel carries a colour bar or a number.

**Placement montage** (`images/placement/`): a grid of 2 rows × 4 columns.
- **Rows**: **T1w** (resampled to the MEGRE grid) on top, **echo-1 magnitude** below.
- **Columns**: axial, coronal and sagittal sections through the SN–VTA region, then an axial section at the dentate level.
- **Red outlines** mark the SN–VTA region of interest (CIT168 atlas, weight ≥ 0.5) in the first three columns, and the left and right dentate nuclei in the fourth column.

**Forced-review flag.** Some items carry a forced-review flag, set by an automated check. The tool shows it as an orange box; the sheet shows it in the `forced_review` column. When the flag is shown, inspect every panel at full size before you grade. The flag itself does not change the grade.

## 4. Rating rules

For each item, record four things.

1. **Gross-failure grade**: 0 (pass), 1 (minor) or 2 (fail), following rubric section 2 below.
2. **Failure codes**: every code that applies, from **M, W, S, B, K, F, A, C, D**.
   - A grade 2 needs at least one code.
   - Record **D** (dental or metal artefact) whenever you see it, whatever the grade.
3. **Placement**: `pass` or `fail`, following rubric section 3 below.
   - Judge placement **only on the T1w and echo-1 magnitude rows of the placement montage**.
   - **Never** use χ, or the gross montage, to judge placement.
   - If placement fails, give the reason in the note.
4. **Note** (optional): a short free-text remark.

Further rules:
- **Rate each item on its own**, against the rubric. You may go back to an item, but do not re-calibrate by comparing items side by side.
- **Repeats.** Some items may appear more than once, under different kit IDs. Rate every item as if you were seeing it for the first time. Do not search for repeats.
- **No expected proportions.** The items are a selection of the study's scans in random order. Do not assume any particular share of passes or fails.
- **Grade only what you see.** Some automated exclusion rules (for example, the automatic grade 2 in the rubric) are applied separately by software. You do not need to apply them, and you cannot see them.
- **Do not try to identify participants, sessions or groups.** Do not open any other file of this project (QC sheets, grade tables, data tables, manuscript drafts) before your ratings are returned. Do not discuss the items with the first author, or with anyone who has seen the automated grades, until then.

## 5. How to use `rate.html`

1. **Unzip the whole kit** and keep its folder structure: `rate.html` must stay next to the `images/` folder.
2. Open `rate.html` in a desktop browser (Chrome, Edge, Firefox or Safari), on a normal-sized screen.
   - The header should read **"Kit check: all 222 images found"**.
   - If images are missing, the kit was not fully unzipped.
3. Type your initials in the header.
4. For each item, the two montages appear side by side, with the rating panel on the right.
   - **Zoom**: use the buttons above each montage (Fit, 100 %, 150 %, 200 %, 300 %). Drag the image to pan.
   - **Full screen**: use the "Full screen" button, or press `G` (gross) or `L` (placement). Press `+` and `−` to zoom, and `Esc` to close.
5. **Keyboard shortcuts** (inactive while you type in the note box):

   | Key | Action |
   |---|---|
   | `0` `1` `2` | grade |
   | `M` `W` `S` `B` `K` `F` `A` `C` `D` | toggle a failure code |
   | `P` | placement: switch between pass and fail |
   | `←` `→` | previous / next item |
   | `N` | next unrated item |
   | `R` | show the rubric |
   | `?` | list the shortcuts |

6. The strip of numbered buttons at the bottom shows your progress:
   - green: grade and placement both entered;
   - amber: started but incomplete.
7. **Saving.** Progress is saved automatically in this browser, for this file location. You can stop and continue later in the same browser.
   - As a backup, press **Export CSV** at each break.
   - If you change browser or computer, or move the folder, use **Import CSV** to restore your last export.
   - If the header says that browser storage is unavailable (for example, in a private window), export the CSV before you close the page.

**Expected time:** about **2 hours** in total, possibly more if you open many items at full size. Breaks are fine, and there is no time limit: please do not rush.

## 6. How to return your ratings

1. When every item is green, press **Export CSV**.
   - A file named `LEADchi_humanQC_<initials>_<date-time>.csv` is saved, usually in your Downloads folder.
   - Its columns are `kit_id, grade, failure_codes, placement, note, timestamp, rater` (`rater` repeats the initials you typed in the header).
   - If automatic saving fails, the CSV text is shown in a box. Copy it into a plain-text file with the same name.
2. Send **only this CSV file** to the study coordinator who sent you the kit. Use the channel the coordinator gave you.
3. **Fallback without the tool.** Fill in `rating_sheet_blank.csv` in a spreadsheet and return it as CSV, keeping the `order` and `kit_id` columns unchanged. Save it as **CSV UTF-8 (comma-separated)**, not as an Excel workbook:
   - `grade`: 0, 1 or 2;
   - `failure_codes`: letters separated by `;`, for example `M;W`;
   - `placement`: `pass` or `fail`;
   - `note`: free text.
4. After you have sent the file, please do not change your ratings. If you notice an error, tell the coordinator, who will record it as a correction.

Questions about the rubric or the tool, but not about specific items, can go to the coordinator at any time.

---

## 7. The rubric (registered text, verbatim)

The two sections below are the operative rubric for this re-rating. They are copied unchanged from the study's registered QC rubric, and the tool's "Rubric" button shows the same text. Section 1 of the registered rubric describes the blinding of the original grading round; it is reproduced in the Appendix for context only. For this re-rating, sections 1 to 6 of this page replace its procedural rules (for example, where grades are entered and who grades).

<!-- BEGIN VERBATIM: qc_rubric.md sections 2-3 -->

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

<!-- END VERBATIM -->

---

## Appendix: registered rubric section 1 (context only; verbatim)

This section describes the blinding of the original grading round. A later amendment replaced its single human rater by the automated rater that this re-rating evaluates. Its procedural rules (who grades, where grades are entered) do not apply to you; sections 1 to 6 of this page do.

<!-- BEGIN VERBATIM: qc_rubric.md section 1 -->

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

<!-- END VERBATIM -->
