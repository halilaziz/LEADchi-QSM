# LEAD-χ: midbrain susceptibility in psychosis-spectrum youth (QSM, OpenNeuro ds007116)

Analysis code, the internally time-stamped protocol, compact results and the visual-QC grades of the study

> Velioglu HA, Altunc AT, Celik MC, Aslan D, Cicek Y, Sakiroglu MA. *Midbrain susceptibility in psychosis-spectrum youth: a procedurally blinded quantitative susceptibility mapping study with ADHD modeled factorially.* (submitted)

## The study in brief

- **Design.** Secondary analysis of the Penn LEAD dataset (OpenNeuro **ds007116**). Quantitative susceptibility mapping (QSM) from multi-echo GRE (MEGRE) with **QSMxT v9.22.0**. Primary outcome: substantia nigra–ventral tegmental area (SN–VTA) magnetic susceptibility χ (ppb). Groups form a **2 × 2 factorial design**: psychosis spectrum (PS) × ADHD, with cells TD, ADHD-only, PS-only and PS+ADHD.
- **Blinding.** All imaging steps ran on coded copies. Labels and the code key were locked with hybrid OpenSSL encryption until every gate, QC grade and freeze had been fixed. The **unlock was on 2026-10-05** (`results/confirmatory/UNLOCK_RECORD.json`).
- **Protocol.** The protocol v1.0 FINAL, amendments 0, 1, 1a, 1b, 1c and 2, the running deviation log and the freeze manifests were frozen internally with SHA-256 and UTC time stamps. The study was not posted on OSF; see "Verifying the protocol" below.
- **Registered analysis** (`results/confirmatory/`): frozen `run_all.R` with the AI visual-QC rater; n = 80; H1 final outcome **B**.
- **Exploratory analyses** (`results/exploratory/`): tier 3 (estimation) and tier 4 (exploratory), plus an adversarial statistical review.
- **Post hoc validation** (`results/posthoc/`):
  - `recovery/`: synthetic-deficit injection and recovery;
  - `engine2/`: a second QSM engine, QSM.jl 0.5.4;
  - `validation_local/`: local validation analyses;
  - `r2s_textconform/`: an R2\* text-conformance check;
  - `humanqc_primary/`: **the analysis reported in the manuscript**. The frozen `run_all.R` was re-run unchanged on the complete human visual-QC grades of the first author (180 sessions graded, **n = 85** analysed).
- **Human visual QC.** Kit 1 (101 sessions plus 10 hidden repeats) and kit 2 (79 sessions plus 21 hidden repeats) were graded blind to group and to item-level AI grades, using random kit IDs. The rater was the first author, so the grading is not independent of the study. See `docs/human_qc/RATER_CHANGE_NOTE.md` and `docs/human_qc/PLAN_kit2.md`.

## Folder map

| Folder | Content |
|---|---|
| `protocol/` | Byte-identical frozen protocol files: `LEADchi_prereg_v1.0_FINAL.md`, `LEADchi_amendment_{0_route,1_FROZEN,1a_G5,1b_G5,1c_QCrater,2_FROZEN}.md` with their original `.sha256` freeze records, `LEADchi_prechi_freeze.md` (pre-χ freeze), `freeze_manifest.json` (data freeze before unlock, amendment 2), `MANIFEST.sha256` (freeze manifest of protocol v1.0), `amendment1_log.md` (running deviation log), `SHA256SUMS`, `verify_protocol.sh` |
| `code/pipeline/` | Processing on the compute server: environment, coding and blinding (`blind_setup.py`, `lock_labels.py`, `unlock_and_merge.py`), QSMxT runs (`run_qsmxt.sh`), registration (`register_rois.sh`), plausibility, ROI extraction, QC montages and κ (`qc_montage.py`, `qc_rubric.md`), gates (`gates_preunlock.py`), usable lists and freeze (`build_usable.py`), DWI steps, audit. `deviations_registered_versions/` holds the registered versions of scripts patched later; each patch is in `protocol/amendment1_log.md` |
| `code/analysis/` | Frozen registered analysis: `run_all.R`, `R/` helpers, `verify_statsmodels.py` (Python cross-check), schema, synthetic demo inputs, Type-I/power simulation. `deviations/` holds the registered versions of `R/secondary.R` and `verify_statsmodels.py` (deviation C4) |
| `code/config/` | Frozen QSMxT configuration (`qsm_primary.toml`, defaults), pipeline choices, software versions, λ and G5 thresholds, ROI hashes, header audit |
| `code/protocol_design/` | Design and power scripts used when writing the protocol (protocol v1.2 delta, E2) |
| `code/qc_ai/` | AI visual-QC rater workflows (passes 1–2), sheet assembly, rubric |
| `code/exploratory/` | Tier 3 / tier 4 scripts and review re-derivations |
| `code/posthoc/` | Post hoc code: `recovery/`, `engine2/` (R analysis plus Python verification), `validation_local/`, `r2s_textconform/`, `humanqc_primary/` (`build_input.py`, `supplement_humanqc.R`, `make_results.py`, `outlier_check_verify_v8.R`) |
| `code/human_qc_kit/` | Human visual-QC kits 1 and 2: build, HTML rating page, verification, simulated-rating test harness (code only), analysis (`analyze_human_qc.py`, `analyze_kit2.py`, `refit_m1_human_qc.R`), rater instructions |
| `code/figures/` | Figure scripts: `figs_humanqc.py` (manuscript Figures 1–3, human-QC sample), `figs_pcn.py` and `figs_pcn_data.R` (same figures, AI-QC sample), `fig2_forest.py` |
| `results/confirmatory/` | Registered run: `results.json`, `report.md`, `manifest.json`, logs, `verify_statsmodels.json`, `UNLOCK_RECORD.json`, `gates/`, `usable_summary.json` |
| `results/exploratory/` | Tier 3, tier 4 and review tables (CSV / JSON / Markdown) |
| `results/posthoc/` | Recovery, engine 2, validation_local, r2s_textconform, and humanqc_primary (`results.json`, `README.md`, frozen re-run, supplement, sanity re-runs) |
| `results/human_qc/` | Kit-1 and kit-2 analysis outputs (agreement, κ, scenario tables, M1 refits) and the kit-1 pre-rating sweep summary |
| `results/qc_grades/` | QC grades (next table) |
| `docs/human_qc/` | Rater-change note and kit-2 analysis plan, both written before rating (English, with a Turkish summary) |

### QC grades (`results/qc_grades/`)

| File | Content |
|---|---|
| `ai/rating_pass1.csv`, `ai/rating_pass2.csv` | AI rater grades per QC montage ID (`qc_id`), with failure codes and free-text notes |
| `ai/grades_final.csv` | Final AI grades per coded session (`code`, `ses_order`) |
| `ai/qc_ids.csv` | `code`, `ses_order` → `qc_id` |
| `ai/kappa.json`, `ai/calibration_kappa.json`, `ai/g4_final.json` | Agreement and the G4 gate |
| `human/kit1_returned_*.csv`, `human/kit2_returned_*.csv` | Human grades as returned by the rater, per kit ID. Kit 2 has two exports; the analysed one is `...T2010Z.csv` (sha256 `fd709028…b767`). It differs from `...T2005Z` in one completed item |
| `human/kit1_key.csv`, `human/kit2_key.csv` | Kit keys: kit ID → `code`, `ses_order`, `qc_id`, AI grade, stratum, sampling weight, repeat links |
| `session_code_map.csv` | `code`, `ses_order` → public OpenNeuro `participant_id`, `session_id` (derived from the unlocked blinding key; see `SANITIZE_LOG.md`) |

To link a human grade to an OpenNeuro session, go kit ID → kit key → (`code`, `ses_order`) → `session_code_map.csv`.

### Path map (original working tree → this repository)

The scripts keep their original relative paths. Some also contain an absolute root, replaced by a placeholder such as `/path/to/Penn_LEAD/LEADchi`. To re-run them, either rebuild the original layout or edit the root constants.

| Original | Here |
|---|---|
| `prereg/*.md`, `prereg/*.sha256`, `mirror/freeze/freeze_manifest.json` | `protocol/` |
| `mirror/scripts/pipeline/` (= `prereg/code/pipeline/`) | `code/pipeline/` |
| `prereg/code/pipeline_deviations/` | `code/pipeline/deviations_registered_versions/` |
| `mirror/scripts/analysis/` | `code/analysis/` |
| `mirror/prereg/` (configs) | `code/config/` |
| `prereg/work/e2_scripts/` | `code/protocol_design/` |
| `prereg/qc_ai/` | `code/qc_ai/` (code), `results/qc_grades/ai/` (grades) |
| `mirror/qc/` | `results/qc_grades/ai/` |
| `results_confirmatory/`, `mirror/gates/` | `results/confirmatory/` |
| `results_exploratory/{code,tier3_estimation,tier4_exploratory,review}` | `code/exploratory/`, `results/exploratory/` |
| `results_posthoc/<x>/code` and `/verification/*.py,*.R` | `code/posthoc/<x>/` |
| `results_posthoc/<x>/` (outputs) | `results/posthoc/<x>/` |
| `human_qc_kit/analysis/`, `human_qc_kit/kit2_2026-10-28/analysis/` | `code/human_qc_kit/kit1/`, `code/human_qc_kit/kit2/` (code); `results/human_qc/` (outputs) |
| `human_qc_kit/**/returned/`, `human_qc_kit/**/_key_DO_NOT_SEND/` | `results/qc_grades/human/` |
| `manuscript/code/{figs_pcn.py,figs_pcn_data.R,fig2_forest.py}`, `results_posthoc/humanqc_primary/code/figs_humanqc.py` | `code/figures/` |

## Software

| Step | Software |
|---|---|
| QSM (primary) | QSMxT 9.22.0 (qsm-core v0.38.0; musl tarball sha256 `6593438c…ca49`): ROMEO, V-SHARP 12 → 1.5 mm, MEDI λ = 0.002 with SMV off, Mask 3. Settings are in `code/config/qsm_primary.toml` and `pipeline_choices.json` |
| QSM (engine 2, post hoc) | QSM.jl 0.5.4 under Julia 1.10.9 (FFTW.jl pinned to 1.8.0) |
| Registration and segmentation | ANTs 2.5.3; SynthSeg and SynthStrip (in the qsiprep 26.0.0 container); FSL 6.0.7.11 |
| Atlases and templates | CIT168 v1.1.0 (probabilistic), SUIT probseg, JHU labels, TemplateFlow MNI152NLin2009cAsym / MNI152NLin6Asym, pyAFQ SCP ROI |
| Python | Python 3 with numpy 2.4.6, nibabel 5.3.3, scipy, DIPY 1.11.0, matplotlib, statsmodels 0.14.5; OpenSSL 3.5.3 for the label lock |
| R | R 4.4.1 with nlme 3.1-164, lme4 1.1-38, sandwich 3.1-1, lmtest 0.9-40, boot 1.3-30, MASS 7.3-61, jsonlite 1.8.9, digest 0.6.37 (see `results/confirmatory/sessionInfo.txt`) |

## How to reproduce

### Data
Get **OpenNeuro ds007116** (CC0), for example `datalad install https://github.com/OpenNeuroDatasets/ds007116.git` followed by `datalad get` on the MEGRE, T1w, DWI and phenotype files. The QC inputs of the PennLINC `EF_desc-T1_qc` / `EF_desc-diffusion_qc` tables are listed in `code/pipeline/README.md`.

### Registered analysis (n = 80, AI visual QC)
1. **Processing.** Follow `code/pipeline/README.md`, section 3, on a Linux server:
   1. `source code/pipeline/env.sh`, after editing `/path/to/...`;
   2. generate your own RSA lock key pair and run `lock_labels.py`;
   3. `blind_setup.py` (coding and header scrubbing);
   4. `run_qsmxt.sh` (install, config, runs);
   5. `register_rois.sh`, `qsm_plausibility.py`, `extract_rois.py`;
   6. `qc_montage.py`, then visual grading;
   7. `build_usable.py`, `gates_preunlock.py --stage final`;
   8. freeze, then `unlock_and_merge.py`.

   Because subject codes are drawn at random, a re-run produces different coded IDs. Our grades can still be applied by session, through `results/qc_grades/session_code_map.csv`.
2. **Analysis** (in `code/analysis/`):
   ```
   Rscript run_all.R --input unlocked/merged_qsmxt_primary_SNVTA.csv --gates gates_preunlock.json \
                     --multiverse unlocked/multiverse_qsmxt_SNVTA.csv --out results/confirmatory
   python3 verify_statsmodels.py --input unlocked/merged_qsmxt_primary_SNVTA.csv --gates gates_preunlock.json \
                     --results results/confirmatory/results.json --out results/confirmatory/verify_statsmodels.json
   ```
   Compare the output with `results/confirmatory/results.json` and `report.md`. The registered draw counts are the defaults; the `--n-*` options are for testing only. `bash run_demo.sh` runs the whole analysis on the synthetic tables in `code/analysis/synthetic/`.

### Human-QC analysis (manuscript, n = 85)
1. Kit analyses:
   - Kit 1: `code/human_qc_kit/kit1/analyze_human_qc.py` with `results/qc_grades/human/kit1_returned_*.csv` and `kit1_key.csv`; outputs in `results/human_qc/kit1_real_run_20261028_*`.
   - Kit 2: `code/human_qc_kit/kit2/analyze_kit2.py --kit2-ratings results/qc_grades/human/kit2_returned_LEADchi_humanQC_rater_20261028T2010Z.csv --kit1-ratings … --key1 … --key2 … --out …`. Its output includes `kit2_scenario_usable.csv` (scenario `H_all180`).
2. `code/posthoc/humanqc_primary/build_input.py --scenario H_all180 --out merged_H_all180.csv` replaces only the columns `usable`, `first_usable`, `fail_reason` and `qc_grade` of the frozen merged table. Its check: rebuilding with scenario `AI_frozen` reproduces the frozen table byte for byte.
3. Run the frozen `code/analysis/run_all.R`, unchanged, on `merged_H_all180.csv` (same gates and multiverse files) → `results/posthoc/humanqc_primary/frozen_run_H_all180/`.
4. `code/posthoc/humanqc_primary/supplement_humanqc.R` → `supplement_H_all180.json`. Then `code/figures/figs_humanqc.py` → Figures 1–3, and `make_results.py` → `results.json` and `README.md`. That README lists every number, with its code path and sanity checks.

### Verifying the protocol
The protocol was frozen and time-stamped internally, not on a public registry. Each file's SHA-256 and UTC freeze time are in its `.sha256` record and in `LEADchi_prechi_freeze.md` / `MANIFEST.sha256`. To check them:
```
bash protocol/verify_protocol.sh     # or: cd protocol && shasum -a 256 -c SHA256SUMS
```
The script checks:
1. every original `.sha256` record;
2. protocol v1.0 FINAL against `MANIFEST.sha256` and against the pre-χ freeze;
3. the repository checksum list;
4. the frozen registered analysis code in `code/analysis/` against `MANIFEST.sha256`;
5. `freeze_manifest.json` against the hash recorded in `UNLOCK_RECORD.json`.

`amendment1_log.md` is a running log. Its post-unlock and post hoc entries were added after the freezes and are labelled as such. The pipeline scripts were sanitized for publication (paths and hostnames only), so their hashes no longer match the freeze manifests; see `SANITIZE_LOG.md`.

## What is not included
- Imaging data and derived maps (NIfTI), QC montages and rater-kit images, figure image files.
- The lock private key and the encrypted label and key files.
- Manuscript text and its build scripts, cover letters, and submission files.
- E-mail drafts, internal working notes, coordinator, statistician and time-evidence packages, AI session transcripts.
- Simulated test-run outputs and superseded or checker folders.
- Per-session merged ROI tables and files larger than 5 MB (for example the human-QC extrapolation draws).
- The server-side scripts of the post hoc engine-2 run (`e2_qsmjl.jl`, `e2_io.py`, `e2_chain.sh`) and of the recovery injection (`inject_sn.py`, `run_rec.sh`). They live only on the compute server; their hashes are given in `results/posthoc/engine2/README.md` and `results/posthoc/recovery/README.md`.

## Licences
- **Code:** MIT (`LICENSE`).
- **Protocol, documentation and result files:** CC BY 4.0 (`LICENSE-docs`).
- **Data:** the source dataset, OpenNeuro ds007116, is released under **CC0**. The QC grades and result tables derived from it are shared here under CC BY 4.0. Please also cite the ds007116 dataset and its original publications.

## Citation
If you use this code, cite the article (`CITATION.cff`):

Velioglu HA, Altunc AT, Celik MC, Aslan D, Cicek Y, Sakiroglu MA. Midbrain susceptibility in psychosis-spectrum youth: a procedurally blinded quantitative susceptibility mapping study with ADHD modeled factorially. (submitted) *[journal, year, DOI to be added]*

Code archive DOI: *[Zenodo DOI to be added]*

## AI assistance
The code in this repository was written with Claude Code (Anthropic) under the first author's direction. The first author specified the analyses, reviewed the code and outputs, and takes responsibility for them. The registered AI visual-QC rater (`code/qc_ai/`) is part of the registered protocol (amendment 1c) and is described there.
