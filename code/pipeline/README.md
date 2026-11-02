# LEAD-χ processing pipeline (registered attachment)

*Status on 2026-09-26: final version for the registration (task F). Written, tested on synthetic data only, and hashed.
Nothing has been run on the server, nothing was written there and nothing was downloaded. No MEGRE image was read.*

These scripts produce the coded MEGRE copy, the QSM, R2* and T2* maps, the T1w-only registrations, the coded ROI
tables, every pre-unlock gate and decision constant, the real-ID DWI tables for S3, the group-blind QC, the frozen
first-usable lists and the single unlock. They implement the registered text `LEADchi_prereg_v1.0_FINAL.md`
(Annexes A and B; Annex H maps each registered rule to its function). **The registered text governs; where a script
deviates from it, that is a coding error handled under Annex C.4.**

They will run on server inside `/path/to/Penn_LEAD/LEADchi` only after the approvals in section 4.
The confirmatory R code (`../analysis/`) reads the tables written by `unlock_and_merge.py`
(schemas `../analysis/schema/roi_table_schema.csv` and `multiverse_schema.csv`).

## 1. Files

| File | Runs on | Role |
|---|---|---|
| `env.sh` | server (sourced) | Redirects HOME, XDG, MATLAB, TemplateFlow, Apptainer, TMPDIR, matplotlib and the QSMxT caches into `$LX/tools/`; `nice -n 10`, ≤ 16 jobs × 2 threads; `ANTS_RANDOM_SEED=20260924`; helpers `lx_init_tree`, `lx_run` (dry-run), `lx_apx`, `lx_home_check`, `lx_assert_inside`. Creates nothing when sourced |
| `lx_common.py` | server + Mac | Write guard, label-file and raw-dataset guards, SHA-256, hash-chained run log, OpenSSL hybrid encryption, header scrubbing with the per-session translation offset, ROI construction, the support S, the two-step reference, the CSF-only reference, reference-domain tissue fractions, Otsu, slab-edge trimming, the ghost-to-signal index, the post-run mask assertion, Q_phase, head orientation |
| `atlas_config.json` | — | Frozen atlas, label and plausibility configuration (CIT168 SHA-256, SUIT dentate, JHU PLIC for description only, SynthSeg labels, envelopes, thresholds, phase-encode axis, slab partitions) |
| `lock_labels.py` | server or Mac (user, once) | Locks `participants.tsv/.json`, every `sub-*_sessions.tsv` and the phenotype TSVs used (PRIME, BDI, ARI, Tanner, substance, CPT, RSAS) into one encrypted file |
| `blind_setup.py` | server (once) | Coded copies of MEGRE, T1w/T2w and the QSIPrep anatomy images with scrubbed headers and a per-session random world translation; phase as float32 radians; JSON reduced to TE/TR/FA/B0 (4 dp); encrypted key (code, IDs, epoch, f0 cluster, offsets); visible flags; draws. `--header-audit` (headers and JSON only): scaling fields, obliquity, sequence parameters, uniqueness of the visible channels |
| `run_qsmxt.sh` | server | QSMxT v9.22.0: install (A1), config (MEDI SMV off; PDF/LBV tolerances from `init`), subsets, Mask 3 (Otsu Mask 2, slab-edge trim), plan, runs (pilot, L-curve, G2, primary, split-echo with λ×2 for {2,4}, multiverse, fallback), sign negation branch. Single pass, ROMEO template mode, obliquity resampling disabled |
| `register_rois.sh` | server | T1w-only registration, route R for every session; atlas warps by one interpolation; G5 inputs; FOV on the T1w grid and in 2009cAsym; `harmref` builds the coverage-harmonised reference region |
| `qsm_plausibility.py` | server | Registered plausibility rules on the support S, the two-step reference, the post-run mask assertion, session indices, CSF-only reference, tissue fractions; flags and aggregate summaries |
| `extract_rois.py` | server | `geom` (volumes, orientation in one frame, Dice, motion index, brain-in-FOV); `g5` (visible G5 flags); `roi` (per-variant coded ROI tables); `multiverse` (9 nodes on their common support) |
| `qc_montage.py` | server | Group-blind montages under random QC IDs, rating sheets (passes 1–4), weighted κ, final grades with failure codes (code D) |
| `qc_rubric.md` | — | Visual QC rubric, single-rater declaration, stop rule, DWI bundle-QC rubric |
| `build_usable.py` | server | Frozen code-keyed first-usable lists per branch; freeze manifest for amendment 2 |
| `pilot_checks.py` | server | Sign check, the L-curve λ rule and the G5 ventral-DC threshold calibration (amendment 1) |
| `gates_preunlock.py` | server | `g2` (engine gate, amendment 1) and `final` (PC1, PC2w, PC3, fallback, G3, PC4, r_L, PC4-DN, pooled constants, cross-session reliability, SD(dref), G4, G5, PC5; amendment 2) → aggregate gate JSON |
| `audit_run.py` | server | Audited execution; rejects any aggregate output holding a coded ID or a per-session numeric vector |
| `dwi_geometry.py` | server | Real-ID DWI node geometry and pSCP_R segment values; `geometry` fixes n_min before posting; `t1t2`; `table` (S3 table after the DWI-QC freeze) |
| `pc5_dwi.py` | server | Gate PC5 on the DWI tables |
| `dwiqc.py` | server | Coded DWI bundle QC: renamed figures, key file, sheet, freeze |
| `s3_candidates_preqc.csv` | — | The 108 pre-QC S3 candidate sessions from public metadata (first T1-eligible MEGRE session with a same-session DWI scalar pass and an SCP_R profile); used only to fix n_min |
| `unlock_and_merge.py` | Mac (user, once) | The single unlock; writes `unlocked/merged_<engine>_<pipeline>_<roi_def>.csv` and `unlocked/multiverse_<engine>_<roi_def>.csv` |
| `tests/test_reference.py` | Mac | Unit tests of the support, the two-step reference, the CSF rule and the helpers (synthetic arrays) |
| `SHA256SUMS` | — | SHA-256 of every file above |

**Written and hashed later (amendment modules, registration A.6 step 6).** `sepia_g2.m` + `sepia_config.m` (only if
A2 is approved; implements Annex B.5 settings verbatim; amendment 1, before any SEPIA χ exists); the D10 ΔR2* job
(only if A3-5 is approved; amendment 1).

## 2. How blinding is implemented

- **Labels.** Group labels, ages and phenotypes exist in the analysis environment only as
  `locked/labels_locked.tar.enc` (`lock_labels.py`). The code↔ID key (with epoch, f0 cluster and the translation
  offsets) exists only as `blind/key.csv.enc` (`blind_setup.py`). Both are encrypted to the user's `lock_pub.pem`.
  **The labels remain public on OpenNeuro and in the first author's own earlier files; the lock is procedural.**
- **Header geometry.** Every coded image of a session carries the same random world translation, so a coded header
  cannot be matched to its raw header; voxel data are unchanged (no resampling).
- **Guards in the code.** No pipeline script opens `participants.tsv`, `*_sessions.tsv` or `phenotype/`
  (`forbid_label_inputs`); after setup no script reads the uncoded MEGRE (`forbid_raw_inputs`); every write is
  refused outside `$LX`.
- **Hidden tables.** Code-keyed continuous values sit in `tables_coded/`; only hashed scripts read them, through
  `audit_run.py`, which rejects any aggregate output that holds a coded ID or a per-session vector.
- **Visible outputs.** Pass/fail flags and aggregates: `flags_coded/`, `gates/`, `freeze/usable_summary.json`, `qc/kappa.json`.
- **Unlock order.** Possible only after amendment 2 (freeze manifest), logged with UTC time, host, user and every hash.

## 3. Order of execution and approvals

`$LX = /path/to/Penn_LEAD/LEADchi`. Every server step begins with `source $LX/scripts/pipeline/env.sh`
and can be previewed with `--dry-run`. Scripts that read code-keyed tables are run as
`$PY audit_run.py --aggregate <json> -- $PY <script> ...`.

**Before posting (no χ, R2* or local field is computed in this phase).**

| # | Where / who | Command | Approval |
|---|---|---|---|
| 0 | Mac, user | E1 §6.1–6.3 actions: delete the uncoded local copy of sub-20139 ses-1, remove or seal the two PNGs, seal `megre_coverage.tsv`; record dates | user action |
| 1 | Mac, user | `openssl genpkey -algorithm RSA -pkeyopt rsa_keygen_bits:3072 -aes-256-cbc -out lock_priv.pem`; `openssl pkey -in lock_priv.pem -pubout -out lock_pub.pem` | user action |
| 2 | Mac → server | Copy this folder to `$LX/scripts/pipeline/`; `lock_pub.pem`, PennLINC `EF_desc-T1_qc.csv` and `EF_desc-diffusion_qc.csv` to `$LX/prereg/inputs/`; check `SHA256SUMS` | A4 |
| 3 | server | `lx_init_tree` | A4 |
| 4 | server | Tool checks (openssl, nibabel/scipy/dipy/matplotlib, ANTs 2.5.3, container inspect) | A6 |
| 5 | server | `$PY blind_setup.py --header-audit` → scaling fields, obliquity, sequence parameters, visible-channel uniqueness, PE axis | read-only |
| 6 | server, user | `$SYSPY lock_labels.py` (once) | user action |
| 7 | server | `$PY blind_setup.py --dry-run`, then `$PY blind_setup.py` (once) | A4, A6 |
| 8 | server | `run_qsmxt.sh --stage install` | **A1** |
| 9 | server, user | Read the captured help/init; write `qsmxt_quality_map_cmd.txt`; sign `qsmxt_flags_verified.txt`; `run_qsmxt.sh --stage config` | A1 |
| 10 | server | Atlas downloads (CIT168, SUIT probseg) | **A3** |
| 11 | server | Copy the two TemplateFlow templates and the pyAFQ SCP_R intermediate ROI into `$LX/rois/` | **A5** |
| 12 | server | `register_rois.sh --stage templates` | A4–A6 |
| 13 | server | `$PY dwi_geometry.py --stage geometry` (fixes n_min; real-ID DWI only) | A4 |
| 14 | server | `run_qsmxt.sh --stage subsets`, `--stage plan` | A1, A4 |
| 15 | Mac | Hash everything → OSF registration | A8 |

**After posting.** register (`--stage all`, `--stage harmref`) → `extract_rois.py --stage geom` → Mask 3 → pilot and
L-curve runs → `pilot_checks.py` (sign, lcurve, g5cal) → `extract_rois.py --stage g5` → G2 runs (`--variant g2`, and `g2alt` if A2 was declined; SEPIA only
if A2) → `gates_preunlock.py --stage g2` → **amendment 1** → primary, split, multiverse → plausibility (indices,
flags) and ROI tables per variant → `extract_rois.py --stage multiverse` → QC passes and κ → `dwiqc.py` setup/freeze
and `dwi_geometry.py --stage table` (and `t1t2`) → `pc5_dwi.py --covariates <age/sex table>` → `build_usable.py --stage usable` → `gates_preunlock.py --stage
final` → `build_usable.py --stage manifest --gates … --include <analysis code> <DWI S3 table>` → **amendment 2** →
mirror to the Mac → `unlock_and_merge.py` (user) → `Rscript ../analysis/run_all.R` and `verify_statsmodels.py`.

## 4. Approvals

| Approval | Covers |
|---|---|
| **A1** | Download of the pinned QSMxT musl tarball (sha256 `6593438c…ca49`) into `$LX/tools/qsmxt/` |
| **A2** (optional) | SEPIA v1.3.0.1 + Cornell MEDI toolbox for G2 (the user completes the Cornell form) |
| **A3** | Atlases: CIT168 v1.1 prob and the SUIT probseg |
| **A4** | Every write and compute step inside `$LX` (coded copy ≈ 15 GB; registrations ≈ 10 GB; QSM 30–110 GB; ≤ 16 nice jobs) |
| **A5** | Reading `~/.cache/templateflow` (two 1 mm templates) and the pyAFQ template ROI, copying them into `$LX/rois/` |
| **A6** | Execute-only use of the qsiprep-26.0.0 container and the AFQ_env/qsirecon python |

## 5. Items to verify at install (header-only or help-text only; logged in amendment 1)

1. Every QSMxT flag used by `build_args` appears in `qsmxt run --help` (the script aborts otherwise), including
   `--no-romeo-individual --romeo-template 2 --obliquity-threshold -1 --medi-data-weighting --medi-max-iter
   --medi-cg-max-iter --medi-tol --lbv-tol`.
2. The `qsmxt quality-map` command line for Mask 2.
3. The TOML key that switches MEDI SMV off, and the PDF/LBV tolerance keys in the `init` file.
4. That `--use-custom-masks mask3` finds the session masks (`--dry`), and the post-run mask assertion on the pilot.
5. The output names (`*_Chimap.nii`, `*_mask.nii`, `*_R2starmap.nii`, `*_localfield.nii`).
6. The ANTs 2.5.3 layout; the SUIT probseg layout (volume = label − 1); the pyAFQ bundle and figure names.
7. The phase value range (`--phase-mode signed12`), the PE axis, the slab partition axis (length 96).
8. `dwi_geometry.py` transform check (template brainstem points inside the session's brainstem) on one session.

## 6. Registered implementation rules (formerly "P5b choices")

B0 4 dp; phase float32 radians; route R for every session; one orientation frame (MNI152NLin2009cAsym affine); the
1-voxel SynthSeg-CSF buffer for every atlas ROI; SynthSeg ROIs eroded per label (6-connectivity), NAcc not eroded;
`roi_vol` = Σw over atlas voxels with w ≥ 0.5 (no mask, no buffer); implausible voxels excluded from the CSF means;
failed ROIs other than the primary ROI written as NA; R2* of every variant from the primary run; slab-coverage flag
forces a review only; G5 VDC = 28 ∪ 60, MEGRE↔T1w Dice inside the MEGRE FOV, SynthSeg `--robust` on echo-1
magnitude; pilot and G2 draws among sessions with a usable T1w, the G2 set disjoint from the pilot with distinct
subjects; the L-curve uses the first 2 pilot sessions; linear κ; worse grade for re-rated items; random QC IDs;
visual grades shared by the fallback and SEPIA branches.

## 7. Tests performed (synthetic data only, session scratchpad)

- **End-to-end on a fake tree** (16 subjects, 26 MEGRE sessions, fake IDs sub-999xx, random arrays):
  header audit → blind setup (offsets identical within a session; geometry preserved; B0 4 dp) → lock → fabricated
  registration and QSM outputs for all 15 variants → geom → plausibility and ROI tables (planted streak flagged) →
  multiverse (234 rows) → montages, sheets, 14-day refusal, κ, finalize with code D → usable lists → pilot sign and
  L-curve (λ file written) → G2 (within-QSMxT branch) → final gates → freeze manifest → unlock (dry-run and real,
  with a throwaway key) → the merged and multiverse tables pass the analysis schema validation in R (95 of 95 schema
  columns present). `audit_run.py` rejected a planted coded ID and a planted 30-element vector.
- **DWI side** on fake bundles/maps (dipy): geometry, t1t2 and table stages; n_min selection; PC5; DWI-QC setup/freeze.
- `tests/test_reference.py`: 12/12 (zero rim, planted ±2 ppm outliers, CSF rule, helpers).
- `run_qsmxt.sh` and `register_rois.sh`: every stage dry-run via `LX_TEST_PL` (refused on server).
- **Not tested:** anything that needs QSMxT, ANTs, the container, real atlases or server.
