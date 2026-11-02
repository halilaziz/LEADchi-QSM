# LEAD-χ amendment 1 (FROZEN): pilot, λ, sign, G5, G2, S3 status, Stage I/II technical log

Frozen 2026-10-04T16:40:59Z, **before the full QSM run**. The only χ in existence: pilot (10 sessions), L-curve (2 sessions × 4 λ), G2 set (20 sessions × 3 pipelines). No group label or code key has been opened.

## Registered decisions fixed by this amendment
- Sign check: **pass** (no phase negation).
- MEDI λ (QSMxT scale): **0.002** (L-curve corner rule); split {2,4}: 2λ = 0.004.
- G5: v1 chain; per-session G5-H1 / G5-DN / IQ-fail (amendments 1a + 1b); VDC threshold 0.80.
- G2: **pass** (engine QSMxT; primary pipeline V-SHARP + MEDI); G2-DN: **pass**.
- S3: **exploratory** (n_min has no solution; Holm = {S1, S2}).

## Full technical log

| # | Step | Result | Deviation? |
|---|---|---|---|
| 1 | Pipeline and analysis code copied to hebb; `sha256sum -c SHA256SUMS` | all OK; prereg FINAL `139a99fb…`, amendment 0 `01206e9e…` match | no |
| 2 | `lx_init_tree`; tool checks; home check | tree created; DIPY 1.11.0, nibabel 5.3.3, OpenSSL 3.5.3, ANTs 2.5.3, qsiprep 26.0.0 (SynthSeg, SynthStrip present); home clean | no |
| 3 | `blind_setup.py --header-audit` (headers and JSON only) | 181 complete sessions, 0 partial. int16 mag/phase; `scl_slope`/`scl_inter` None; 120×160×96 at 1.5 mm; qform = sform = 1; 0 NIfTI extensions; `descrip` contains Time in all (it is scrubbed by blind setup). Obliquity 0.0° in all; PE direction i (ROW); Adaptive Combine; PF 1; GRAPPA 2; `*fl3d4r`; syngo MR E11 in all. B0 clusters 2.8929 T (22) / 2.8935 T (44) / 2.8951 T (115). Visible-channel uniqueness: 0 unique sessions, 12 combinations. Protocol names: `QSM_SWI_1.5mm` 16, `qsm_acq-1.5mm_GRE` 165. `header_audit.json` sha256 `093ab496…6d5c` | no |
| 4 | QSMxT install (A1) | v9.22.0, qsm-core v0.38.0, tarball SHA-256 matches the pin `6593438c…ca49`; help, init and quality-map help captured and hashed (`frozen.sha256`) | no |
| 5 | Flag verification | every QSMxT flag used by `build_args` is present in `qsmxt_run_help.txt`. **Sign-off by Claude on the user's behalf** under the standing approval (amendment 0 item 3), instead of by the user in person | **yes, procedural (C.2)** |
| 6 | Quality-map command | `qsmxt quality-map --magnitude {MAG1} --phase2 {PHASE2} --te1 0.0075 --te2 0.015 --output {OUT} {PHASE1}` (from the captured help; TE from the header audit) | no (fills a VERIFY) |
| 7 | `--stage config` | `qsm_primary.toml` = defaults with `[inversion.medi] smv = false` (single key). The wrapper's flat-key regex found 0 PDF/LBV tolerance keys and aborted. `pipeline_choices.json` was then written with a section-aware parser: `pdf_tol = 1e-5` (`[bg_removal.pdf] tol`), `lbv_tol = 1e-6` (`[bg_removal.lbv] tol`), the values in the captured init file as registered. The rest of the content is identical to the wrapper's; hashed in `frozen.sha256` | **yes, technical (C.2)**; the values are as registered |
| 8 | Atlases (A3) | CIT168 v1.1.0 prob: SHA-256 matches the pin `7a91c367…0f66`; 193×229×193×16, 1 mm. SUIT `atl-Anatom_space-MNI_probseg.nii`: sha256 `bbc3b018…d259`, 153×103×84×**34** uint8, which **confirms volume = label − 1** (VERIFY resolved); `atl-Anatom.tsv` labels 29 L_Dentate, 30 R_Dentate, 31 L_Interposed, 32 R_Interposed | no |
| 9 | Template and ROI copies (A5) | TemplateFlow 1 mm T1w and brain masks for MNI152NLin2009cAsym and MNI152NLin6Asym (non-empty), plus pyAFQ `SCP_R_inter_prob.nii.gz`, copied read-only from outside into `rois/`; hashes in `prereg/rois_sha256.txt` | no |
| 10a | `register_rois.sh --stage templates` | rc 0, 10:02–10:51Z. SynthStrip/SynthSeg on the 2009c template; NLin6→2009c SyN (ANTs seed 20260924); CIT168 subset (SHA check passed); SUIT probseg subset (sha256 `bbc3b018…d259`); JHU labels sha256 `3c3f5d2f…0ec16`, PLIC voxels 7,506 (labels 19/20 present, VERIFY resolved); home clean | no |
| 10b | `dwi_geometry.py --stage geometry` (real-ID DWI only; no MEGRE, χ or labels) | rc 0, 11:16Z. 217 sessions, all "ok" (every transform-direction check passed). 108/108 pre-QC S3 candidates processed. Split-half Spearman–Brown for SCP_R tw-ICVF 0.983. **n_min has no solution:** sessions keeping ≥ n retained nodes: n = 10: 6/108, n = 8: 6/108, n = 6: 9/108, against the ≥ 90% rule → `fallback_required: true` | **yes, registered-rule outcome** |
| 10c | Why (label-free node table, SCP_R nodes 0–30, 217 scans) | Rule (ii), mapped-dentate distance < 4 mm, removes 99–87% of nodes 0–16 (median distance 0.5–2.1 mm) and ≥ 69% up to node 19. Rule (iii), ISOVF > 0.30, removes ≥ 51% from node 19 and 82–92% from node 24 onward (median ISOVF 0.44–0.48, next to the fourth ventricle). 4% of node observations pass both; median retained nodes per scan 0 (mean 1.2). Rule (i) was not limiting (n_int ≈ 36–40, so the limit stays 30) | — |
| 10d | **Decision (pre-χ, label-free): S3 is exploratory.** Holm family = {S1, S2} at .025 / .05 | The registered text foresees S3 as exploratory when the geometry step fails (§12.3 "Fallback") and when PC4-DN, G2-DN or PC5 fail. The present case (the job ran but no n_min satisfies the rule) is not named explicitly; the closest registered consequence is applied. S3 is estimated as registered, with n_min = 6 as the DWI-multiverse level (§24.5) and the "fixed 0–30 with ISOVF only" node rule as a further level. It is reported as exploratory, with the reason: the proximal right SCP cannot be sampled free of dentate proximity and fourth-ventricle CSF at 1.7 mm DWI. H1, K1–K4, S1, S2 and T1 are unaffected | **yes, decision recorded before any χ** |
| 11 | Lock key pair | Generated by the user on the Mac (2026-09-29 13:05 local). `lock_pub.pem`: RSA 3072, sha256 `09b15301…8065`. The private key was never seen by the AI | no |
| 12 | `lock_labels.py` | Run by Claude under the standing approval (registered: run by the user). 282 files locked: participants 1 file, sessions.tsv 132, phenotype 16. Hashes: `labels_locked.tar.enc` `6f0af846…ce74`; `labels_locked.sessionkey.enc` `ca9bf0e0…7fe2`; manifest `fbca1864…c889c` | **yes, procedural (C.2)** (who ran it) |
| 13 | `blind_setup.py`, first run | Stopped after 60/181 sessions. For sub-21038 ses-1 the echo-1 phase range did not match `--phase-mode signed12`. The partial output (`blind/anat_all`, `blind/bids_all`, 3.6 GB, no key written) was removed, as the script instructs | — |
| 14 | Phase-encoding check (min/max only, the same technical QC as round-1 M6; no group information) | 700 files signed 12-bit (−4096…4094); **24 files (6 sessions) unsigned 12-bit (0…4095)**: sub-20303 ses-3, sub-20712 ses-1, sub-20818 ses-1, sub-20823 ses-1, sub-20866 ses-1, sub-21038 ses-1. Every echo of a session uses the same encoding | — |
| 15 | `blind_setup.py` patched with `--phase-mode auto12` | Per file: signed12 if min < 0, else unsigned12. The existing registered conversions are used: signed v·π/4096; unsigned (v − 2048)·π/2048. The two encodings are mutually exclusive, all echoes of a session must share one, and no per-session mode is written to any visible table. Registered file sha256 `76923c53…cd6a` (kept as `code/pipeline_deviations/blind_setup_v1.0_registered.py`); patched `6e292e91…820a`. Re-run started 10:14:45Z | **yes, technical (C.2)**; the phase conversion rules themselves are as registered |

## 2026-09-29 → 10-02: post-freeze, pre-full-run steps

The pre-χ freeze is `LEADchi_prechi_freeze.md`, sha256 `e909deec…a319`, frozen 2026-09-29.

| # | Step | Result | Deviation? |
|---|---|---|---|
| 16 | `blind_setup.py --phase-mode auto12` | rc 0 (10:14–10:24Z); 181 sessions coded (9.7 GB). `key.csv.enc` `e33f641b…2f51`; `key.sessionkey.enc` `6825cef1…7a0f`; pilot 10 / L-curve 2 / G2 20 draws | as in #15 |
| 17 | `run_qsmxt.sh --stage subsets`, `--stage plan` | subsets built; plan (QSMxT `--dry`) OK; B0 read from the coded sidecar (e.g. 2.8935 T; the plan display rounds to 2.9) | no |
| 18 | Registration chain v1: `register_rois.sh --stage all`, `harmref`, `geom`, `mask3` | First launch failed because the `.sh` exec bit was lost in transfer (fixed with chmod inside LEADchi; nothing written). Re-run 09-29 13:44Z → 09-30 00:33Z: 180 registered, 0 failed; harmref 1,881,049 template voxels; geom 181; Mask 3 for 181 | technical, trivial |
| 19 | `run_qsmxt.sh` obliquity flag | QSMxT (clap) rejected `--obliquity-threshold -1`; changed to `--obliquity-threshold=-1` (same registered value), and `check_flags` now strips `=value` | **yes, technical (C.2)**; value unchanged |
| 20 | Q_phase ordering | `extract_rois --stage roi` read Q_phase from the primary-run indices, which do not exist before amendment 1. Q_phase depends only on the phase and Mask 3, so it is now also computed for the pilot/g2 variants and read from the variant's own indices when the primary file is absent | **yes, technical (C.2)** |
| 21 | Pilot (10 sessions), L-curve (2 × 4 λ) | Sign check **pass**: pooled medians GP 92.2, SN–VTA 54.0, RN 34.9, DN_R 32.9 ppb. L-curve: both sessions chose **λ = 0.002** (default 0.001); ρ is non-monotonic in λ (possible incomplete MEDI convergence at max_iter 30), and the registered corner rule was applied as written. G5 VDC calibration: 20% of the pilot below 0.80 (not > 20%) → **threshold 0.80** | branch taken |
| 22 | `extract_rois --stage g5` | Summary bug: NA in `g5_pass` for the unregistered session; fixed to count NA as 0. The partial flag file from the crashed run was replaced | technical (C.2) |
| 23 | G5 outcome (v1) | **46% fail** (83/180). Diagnosis (aggregate-only, via `audit_run`; plus coded overlays of 3 gross and 1 good session, with no χ): 20 sessions (11%) are motion-corrupted echo-1 magnitude where SynthSeg fails (image quality, not registration); cerebellar-WM Dice median 0.81 | — |
| 24 | **Amendment 1a (G5), user decision 2026-10-02**, `LEADchi_amendment_1a_G5.md` sha256 `cb24cad7…9faf` (09:53:47Z) | Menu revision: N4 on T1w + mask dilation 4 (MI kept); IQ-fail rule; G5-H1 / G5-DN split; chain trigger on non-IQ-fail sessions. The revised chain was launched 10-02 09:55Z; v1 outputs kept as `*_v1` | **yes, Stage II (C.3)**; registered |

| 25 | **Amendment 1b (G5), user decision 2026-10-03**, `LEADchi_amendment_1b_G5.md` sha256 `164e9092…c45c` (09:05:23Z) | The 1a revision (v2a) worsened G5: T1w→MNI Dice fell because the dilated mask also feeds the SyN, while the MEGRE↔T1w Dice were unchanged. Reverted to the registered v1 chain (`register_rois.sh` `ebe8b6d1…`; v2a kept as `reg_v2a` and `register_rois_v2a_amendment1a.sh` `2dc98410…`). The 1a G5 rules apply per session; the chain trigger is removed; bias safeguards added. Applied 10-04 16:29–16:30Z (connection outage 10-03, during which nothing ran) | **yes, deviation (C.3)** |
| 26 | v1 G5 under the 1a rules (`flags_coded/g5_summary.json`) | Registered 180; **G5-H1 pass 136**; IQ fail 26; non-IQ fail share 11.7% (informative only; no trigger after 1b); **G5-DN pass 97**; VDC threshold 0.80 (pilot rule re-applied: 20% below 0.80, not > 20%) | branch taken |
| 27 | **G2 engine gate** (A2 declined → within-QSMxT branch; `gates/g2.json`) | SN–VTA ICC(3,1) across V-SHARP+MEDI / V-SHARP+RTS / PDF+iLSQR = **0.889** [0.783, 0.950], n = 20 → **G2 pass**. Dentate ICC(3,1) = **0.771** [0.501, 0.925], n = 11 complete (dentate NA where G5-DN fails) → **G2-DN pass**. Engine: QSMxT, λ = 0.002 | branch taken |


- n_min and the transform check (step 10).
- JHU PLIC label check; SynthSeg version.
- Lock and blind setup: needs `lock_pub.pem` from the user.
- Pilot λ (L-curve), sign check, G5 calibration and the G2 result: these come after blind setup and the pilot, as registered.

## Hashes of the gate outputs
```
c01fcb1415f0b6acf4e64372313979c08c066b743d80a3bce94c090b4598709c  gates/g2.json
37d633020da688349308a419e00ccf53eae49577e48d78d656f4faf9d2fe7ea3  gates/g5_cmp_v1v2.json
9aff6fa514c726e59f75a0d51255487ef5e19313f42eb50a8be306cd128639b5  gates/g5_diag2_v1.json
c7cd9134584ad828f2cc796d9520c4fab2eec936e0f4523d910dcbca49fabff4  gates/g5_diag3_v1.json
97a67f96b933fb1c0daec7a5c170bcb30fde4b6da805eddd31daf26aaa409017  gates/g5_diag3_v2.json
abc161e4e7768563dc1a1cec990595a7ecb577c2c6e99ed179bdfee6aacd8c31  gates/g5_diag_v1.json
dd6c0b8dd2e8fa917eded4eddce1dc2e5e674781f86628dc13e6a27f01eb2342  gates/pilot_g5cal_qsmxt.json
83c384eb0995c5bcd4bf6f4166b1ae991f7f1fe92fb11815d15720eb0229bbea  gates/pilot_g5cal_qsmxt_v1.json
5f7f28309ad2c4dac20d045f4bf8b5ca216cc126f053dc32bedfbfbd87554a1c  gates/pilot_g5cal_qsmxt_v2a.json
dec34f384c03396f04a092e3239919792ca19de921d0a813c25ea451fc787bfe  gates/pilot_lcurve_qsmxt.json
ded7156c491a31a18a0c0478fdc5a796af74ece96f72b52f0a9686501018b500  gates/pilot_sign_qsmxt.json
16528de27a77ece8e56974a0fc49daf07a1905fff99cedc76480f9608f5532af  flags_coded/g5_summary.json
ee8eef6273e51ddf6aacc8fb3f6199c865f5bfef5936afbbcc4b737a2ac28fb6  prereg/amendment1_lambda_qsmxt.txt
e38430d35e1bee73ba0393cb4a499583fe38646227eae1ff8728bdbcda9cf9d1  prereg/g5_vdc_threshold.txt
9cfb8c7b438f8e9d42667ca1b0541e5022f9505f23db0836e48fb8c576e4cfe7  prereg/frozen.sha256
dd080b0d4a73f1441b454b62243edf65f86c1a18f942360860f7718692f458a4  tables_coded/session_geom.csv
528a4b28959cc3ec43a5211ab7668a69f145773c5c6ee08621e9294dd8e625fe  scripts/pipeline/audit_run.py
6e292e915413964e3a49a32c2fc668fdc4a35add3345d2e5a2f14f8c4913820a  scripts/pipeline/blind_setup.py
47a17cd6bf9c5a1c9359141fb4c09c6c7d53eb1819c1213d96fff07f13ded5b1  scripts/pipeline/build_usable.py
d5ab6b43afc43a8f070fedc771f8ceeab5b4d37fa01d410eb4036a8d2a0eabc8  scripts/pipeline/dwi_geometry.py
1e15f70cd7738ad492ddc22fedb5cee2480a8bb6b9eb6546f1e41a3262057f61  scripts/pipeline/dwiqc.py
23592ebb19c69ee3425851d29c4847379952d520c1054ac112d72a15bf9ca693  scripts/pipeline/extract_rois.py
4fbf220cbdcce0009a4c9d37e9b956e37da74cb7a1ba7e5471700297a785beb3  scripts/pipeline/gates_preunlock.py
46fcd22f8f2ed5d527f1f7ffefb1be8a0b5992ef53e248b8f47f550ab9357de5  scripts/pipeline/lock_labels.py
ebad4364ee9ae071b2d929b12aa40eb4177c16d9bf16e6bd25fdbef326147689  scripts/pipeline/lx_common.py
b795b10f4d3403928da07b1c6c1691f6e3287c5c8814e86f4860e5a694df82d4  scripts/pipeline/pc5_dwi.py
bd11f9112094829e156270cbed5cdf5a3d1d77863d925bd0eaf2a6c9b6db2b84  scripts/pipeline/pilot_checks.py
bf6fb6b3c57e4dba11b47b77d4b420b0a44504a2c8d039df9adccee164e27450  scripts/pipeline/qc_montage.py
9426b310c4697ed549c82abe448d1666dc69d9ba1331b0b2e5ad2e915c19550d  scripts/pipeline/qsm_plausibility.py
7453913955e5474c8529d8e8f372b8bcf81cb5cd8a3446a801b0d05d0d667d76  scripts/pipeline/unlock_and_merge.py
ff8482b2d14b6b98201b6892287369515b7f8db1288e2f7bd0bc5d0dcf2f7809  scripts/pipeline/env.sh
ebe8b6d128fbcb907d13c685f5e56930a52cc38656402caefa300fa252240be2  scripts/pipeline/register_rois.sh
7a2e642b108e30e85c2daeac5fdda7b2ad0a670e7018145d516539ee18138478  scripts/pipeline/run_qsmxt.sh
```
