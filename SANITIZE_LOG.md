# Sanitize log

Copies in this repository were sanitized on 2026-11-02 so that no server paths, hostnames, usernames or e-mail addresses are published. Only literal strings were replaced; no program logic was changed. Paths become neutral placeholders under `/path/to/...`; set them to your own locations before running anything.

**Not sanitized:** the files in `protocol/` are byte-identical copies of the internally time-stamped originals, so their SHA-256 values still verify (`protocol/verify_protocol.sh`). They still contain the original server paths and hostnames.

**Consequence for frozen code hashes.** The frozen pipeline scripts in `code/pipeline/` and `code/pipeline/deviations_registered_versions/`, and `code/analysis/R/synthetic.R`, differ from their frozen originals only by the substitutions listed below, so their SHA-256 no longer matches `protocol/MANIFEST.sha256` or `protocol/LEADchi_prechi_freeze.md`. The registered analysis code (`code/analysis/R/*.R` except `synthetic.R`, `run_all.R`, the schema, simulation and synthetic files, and the registered versions of `secondary.R` and `verify_statsmodels.py` in `code/analysis/deviations/`) contained no such strings and still matches the frozen manifest byte for byte; `protocol/verify_protocol.sh` checks it. The per-folder `SHA256SUMS` files of the original working tree were not copied, because they list files that are not shared here (images, data) or files changed here.

## Rules (applied in this order)

The original strings are shown in redacted form so that this log does not itself republish them; the unredacted originals remain only in the byte-identical `protocol/` files.

| rule | original string (shown redacted) | replacement | what it was | replacements |
|---|---|---|---|---|
| R1 | `<tmp>/claude-<uid>/<project>/<session-uuid>/scratchpad` | `/path/to/scratch` | AI-session scratch directory | 7 |
| R2 | `<server-volume>/<user>/<dir>/Penn_LEAD` | `/path/to/Penn_LEAD` | server project root | 9 |
| R3 | `<server-volume>/<user>/apptainers` | `/path/to/apptainers` | server container/env directory | 3 |
| R4 | `<server-volume>/<user>` (any other) | `/path/to/server_home` | server user directory | 0 |
| R5 | `<mac-home-root>/<user>/Downloads/Penn_LEAD` | `/path/to/Penn_LEAD` | local (Mac) project root | 143 |
| R6 | `<mac-home-root>/<user>` (any other) | `/path/to/home` | local home directory | 0 |
| R7 | `/misc/opt` | `/path/to/opt` | server software prefix | 2 |
| R8 | the account username (whole word) | `<user>` | username | 1 |
| R9 | `<name>.local` machine hostname | `<local-host>` | local machine hostname | 1 |
| R10 | compute-server hostname 1 (lower case) | `server` | server hostname | 61 |
| R11 | compute-server hostname 1 (capitalised) | `Server` | server hostname | 0 |
| R12 | compute-server hostname 2 | `server2` | server hostname | 0 |
| R13 | gateway hostname | `gateway` | gateway hostname | 0 |
| R14 | VPN product name (any case) | `VPN` | VPN product name | 0 |
| R15 | personal e-mail addresses | `<email>` | e-mail address | 0 |

Total: 227 replacements in 107 files.

## Every substitution, per file

| file | rule → replacement | count |
|---|---|---|
| `code/analysis/R/synthetic.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/exploratory/review_check_statsmodels.py` | R5 → `/path/to/Penn_LEAD` | 4 |
| `code/exploratory/review_rederive.R` | R5 → `/path/to/Penn_LEAD` | 4 |
| `code/exploratory/review_split_echo.R` | R5 → `/path/to/Penn_LEAD` | 4 |
| `code/exploratory/t3a_tier3_estimation.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/exploratory/t3a_verify_statsmodels.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/exploratory/tier4_exploratory.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/figures/figs_pcn_data.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit1/analyze_human_qc.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit1/build_kit.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit1/build_log.json` | R5 → `/path/to/Penn_LEAD` | 3 |
| `code/human_qc_kit/kit1/make_html.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/human_qc_kit/kit1/make_readme_rater.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/human_qc_kit/kit1/make_zip.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit1/simulate_ratings.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit1/verify_kit.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/human_qc_kit/kit2/analyze_kit2.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit2/build_kit2.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit2/build_log_kit2.json` | R5 → `/path/to/Penn_LEAD` | 6 |
| `code/human_qc_kit/kit2/make_html_kit2.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/human_qc_kit/kit2/make_readme_rater_kit2.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/human_qc_kit/kit2/make_zip_kit2.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit2/run_tests_SIMULATED.sh` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit2/simulate_kit2_ratings.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/human_qc_kit/kit2/verify_kit2.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/pipeline/README.md` | R2 → `/path/to/Penn_LEAD` | 2 |
| `code/pipeline/README.md` | R10 → `server` | 33 |
| `code/pipeline/atlas_config.json` | R3 → `/path/to/apptainers` | 1 |
| `code/pipeline/atlas_config.json` | R2 → `/path/to/Penn_LEAD` | 1 |
| `code/pipeline/atlas_config.json` | R10 → `server` | 1 |
| `code/pipeline/blind_setup.py` | R10 → `server` | 1 |
| `code/pipeline/build_usable.py` | R10 → `server` | 1 |
| `code/pipeline/deviations_registered_versions/blind_setup_v1.0_registered.py` | R10 → `server` | 1 |
| `code/pipeline/deviations_registered_versions/extract_rois_v1.0_registered.py` | R10 → `server` | 1 |
| `code/pipeline/deviations_registered_versions/extract_rois_v1.0b_before_1a.py` | R10 → `server` | 1 |
| `code/pipeline/deviations_registered_versions/qc_montage_v1.0_registered.py` | R10 → `server` | 1 |
| `code/pipeline/deviations_registered_versions/qsm_plausibility_v1.0_registered.py` | R10 → `server` | 1 |
| `code/pipeline/deviations_registered_versions/register_rois_v1.0_registered.sh` | R10 → `server` | 1 |
| `code/pipeline/deviations_registered_versions/register_rois_v2a_amendment1a.sh` | R10 → `server` | 1 |
| `code/pipeline/deviations_registered_versions/run_qsmxt_v1.0_registered.sh` | R10 → `server` | 1 |
| `code/pipeline/dwi_geometry.py` | R10 → `server` | 1 |
| `code/pipeline/env.sh` | R3 → `/path/to/apptainers` | 2 |
| `code/pipeline/env.sh` | R2 → `/path/to/Penn_LEAD` | 3 |
| `code/pipeline/env.sh` | R7 → `/path/to/opt` | 2 |
| `code/pipeline/env.sh` | R10 → `server` | 4 |
| `code/pipeline/extract_rois.py` | R10 → `server` | 1 |
| `code/pipeline/lock_labels.py` | R10 → `server` | 1 |
| `code/pipeline/lx_common.py` | R2 → `/path/to/Penn_LEAD` | 2 |
| `code/pipeline/qc_montage.py` | R10 → `server` | 1 |
| `code/pipeline/qsm_plausibility.py` | R10 → `server` | 1 |
| `code/pipeline/register_rois.sh` | R10 → `server` | 1 |
| `code/pipeline/run_qsmxt.sh` | R10 → `server` | 1 |
| `code/pipeline/unlock_and_merge.py` | R10 → `server` | 1 |
| `code/posthoc/engine2/engine2_analysis.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/humanqc_primary/build_input.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/humanqc_primary/outlier_check_verify_v8.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/humanqc_primary/supplement_humanqc.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/r2s_textconform/k2_pc3_textconform.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/recovery/recovery_analysis.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/recovery/recovery_analysis.R` | R10 → `server` | 1 |
| `code/posthoc/recovery/verification/verify_recovery.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/validation_local/validation_local.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/validation_local/verification/check_readme.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/posthoc/validation_local/verification/diag_nlme_optimizer.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/validation_local/verification/rerun_2026-10-07/check_new_outputs.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/posthoc/validation_local/verification/rerun_2026-10-07/check_readme_v2.py` | R5 → `/path/to/Penn_LEAD` | 4 |
| `code/posthoc/validation_local/verification/rerun_2026-10-07/diag_nlminb_tolerance.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/validation_local/verification/rerun_2026-10-07/diag_reml_criterion.R` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/validation_local/verification/verify_independent.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/posthoc/validation_local/verify_statsmodels.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/protocol_design/e2_design_v12.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/protocol_design/e2_design_v12.py` | R1 → `/path/to/scratch` | 1 |
| `code/protocol_design/e2_hc2_vs_hc3.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/protocol_design/e2_hc3_afterqc.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/protocol_design/e2_hc3_power_check.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/protocol_design/e2_hc3_power_check.py` | R1 → `/path/to/scratch` | 1 |
| `code/protocol_design/e2_labels_hc3.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/protocol_design/e2_labels_hc3.py` | R1 → `/path/to/scratch` | 1 |
| `code/protocol_design/e2_lit_ranges.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/protocol_design/e2_pc2w_states.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/protocol_design/e2_power_v2_rerun.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/protocol_design/e2_power_v2_rerun.py` | R1 → `/path/to/scratch` | 1 |
| `code/protocol_design/e2_robust_bound_sim.py` | R5 → `/path/to/Penn_LEAD` | 1 |
| `code/protocol_design/e2_robust_bound_sim.py` | R1 → `/path/to/scratch` | 1 |
| `code/protocol_design/e2_s1_hc3.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/protocol_design/e2_scpr_reliability.py` | R5 → `/path/to/Penn_LEAD` | 2 |
| `code/qc_ai/pass1_workflow.js` | R1 → `/path/to/scratch` | 1 |
| `code/qc_ai/pass2_workflow.js` | R1 → `/path/to/scratch` | 1 |
| `results/confirmatory/UNLOCK_RECORD.json` | R5 → `/path/to/Penn_LEAD` | 7 |
| `results/confirmatory/UNLOCK_RECORD.json` | R9 → `<local-host>` | 1 |
| `results/confirmatory/UNLOCK_RECORD.json` | R8 → `<user>` | 1 |
| `results/confirmatory/manifest.json` | R5 → `/path/to/Penn_LEAD` | 3 |
| `results/confirmatory/report.md` | R5 → `/path/to/Penn_LEAD` | 3 |
| `results/confirmatory/results.json` | R5 → `/path/to/Penn_LEAD` | 3 |
| `results/confirmatory/run_all.log` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/confirmatory/run_all_attempt1_halted.log` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/exploratory/review/split24_server_check.txt` | R2 → `/path/to/Penn_LEAD` | 1 |
| `results/exploratory/review/split24_server_check.txt` | R10 → `server` | 1 |
| `results/exploratory/tier3_estimation/t3a_console.txt` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/exploratory/tier3_estimation/t3a_log.txt` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/exploratory/tier3_estimation/tier3_estimation.json` | R5 → `/path/to/Penn_LEAD` | 9 |
| `results/exploratory/tier4_exploratory/results_tier4.json` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/exploratory/tier4_exploratory/run.log` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/human_qc/kit1_real_run_20261028_noexcl/human_qc_report.md` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/human_qc/kit1_real_run_20261028_noexcl/human_qc_results.json` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/human_qc/kit1_real_run_20261028_primary_exclA/human_qc_report.md` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/human_qc/kit1_real_run_20261028_primary_exclA/human_qc_results.json` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/human_qc/kit2_real_run_20261028/kit2_results.json` | R5 → `/path/to/Penn_LEAD` | 2 |
| `results/posthoc/engine2/README.md` | R10 → `server` | 1 |
| `results/posthoc/engine2/verification/VERIFY.md` | R10 → `server` | 1 |
| `results/posthoc/humanqc_primary/frozen_run_H_all180/manifest.json` | R5 → `/path/to/Penn_LEAD` | 3 |
| `results/posthoc/humanqc_primary/frozen_run_H_all180/report.md` | R5 → `/path/to/Penn_LEAD` | 2 |
| `results/posthoc/humanqc_primary/frozen_run_H_all180/results.json` | R5 → `/path/to/Penn_LEAD` | 3 |
| `results/posthoc/humanqc_primary/sanity_frozen_run_AI/manifest.json` | R5 → `/path/to/Penn_LEAD` | 3 |
| `results/posthoc/humanqc_primary/sanity_frozen_run_AI/report.md` | R5 → `/path/to/Penn_LEAD` | 2 |
| `results/posthoc/humanqc_primary/sanity_frozen_run_AI/results.json` | R5 → `/path/to/Penn_LEAD` | 3 |
| `results/posthoc/humanqc_primary/sanity_supplement_AI/supplement_AI_frozen.json` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/posthoc/humanqc_primary/supplement_H_all180.json` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/posthoc/recovery/README.md` | R10 → `server` | 1 |
| `results/posthoc/validation_local/results.json` | R5 → `/path/to/Penn_LEAD` | 3 |
| `results/posthoc/validation_local/run.log` | R5 → `/path/to/Penn_LEAD` | 1 |
| `results/posthoc/validation_local/verification/VERIFY.md` | R5 → `/path/to/Penn_LEAD` | 1 |

## Other changes to copies

- `results/qc_grades/session_code_map.csv` is not a copy. It was derived from the unlocked blinding key (`mirror/unlocked/key_plain.csv` in the working tree) by keeping only the columns `code`, `ses_order`, `participant_id` and `session_id`, so that the coded QC grades can be linked to public OpenNeuro ds007116 subject and session IDs. The per-session blinding offsets and B0-cluster columns were dropped.
- Files were moved into the repository layout (see the path map in `README.md`). Apart from the substitutions above, their contents are unchanged.

## Added 2026-11-02: server-side post hoc scripts
Copied read-only from the compute server (original SHA-256 matched the recorded hashes). In `e2_chain.sh`, `e2_io.py` and `run_rec.sh`, the absolute server project root was replaced by `/path/to/Penn_LEAD` (3 replacements). No logic changed.
