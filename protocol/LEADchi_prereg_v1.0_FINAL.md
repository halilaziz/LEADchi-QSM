# LEAD-χ: Preregistration of a secondary data analysis (v1.0 FINAL)

> **Registration metadata**
> - **OSF form:** "Secondary Data Preregistration" (OSF registration schema v3, id `64775783798e08000a70407e`), implementing the template of Van den Akker et al. (2021): 25 questions in five parts plus a statement of integrity. Numbering Q1–Q25 follows the template. Item headings are paraphrased; the finaliser checks each heading against the live form when pasting [UNVERIFIED wording].
> - **Type:** embargoed registration [FILL: embargo end date, at most 4 years; recommended: lifted at the preprint]. COS category: **"Registration prior to analysis of the data"** (Q18).
> - **Version 1.0 FINAL**, prepared 2026-09-26 to 2026-09-29. It consolidates protocol v1.1, the v1.2 change set (E2), the prior-exposure audit (E1/E4), the technical pinning file (E3), and the fixes required by three independent reviews of the v1.0 draft (editor, QSM physics, biostatistics; attachments R1–R3). `CHANGELOG.md` lists every change from the v1.0 draft. Annex E gives the hashes of every attachment.
> - **Precedence (registered).** *The registered text governs. Hashed code that deviates from it is a coding error, handled under Annex C.4: both results are shown and the text-conforming result is decisive.* Annex H maps every registered rule to the file and function that implements it.
> - **State of the data at posting.** No susceptibility (χ) map, local-field map or ROI value (χ or R2\*) has been computed for any Penn LEAD session. Whole-head voxelwise R2\* was fitted in memory for feasibility in 4 sessions and summarised as whole-head percentiles only (Q18, rows M8 and M10). The first χ map is computed only after the OSF timestamp.
> - **Placeholders** marked `[FILL: …]` are completed by the finaliser before posting (Annex G). No other text changes between this version and the posted one.
> - **Amendments under embargo.** Amendments 1 and 2 are posted as updates of this OSF registration (or, if the form does not allow updates, as time-stamped files in the same embargoed OSF project, linked from the registration). The OSF timestamps are the evidential record; commit dates of the private GitHub repository are author-controlled and are not relied on.
> - **Cross-references.** "Q5"–"Q25" are template items. The four H1 claim qualifiers are called **K1–K4** (renamed from "Q1–Q4" of the draft to avoid a clash with the template numbering). Sections inside an item are cited as §4.1, §12.4 and so on.
> - **Structure.** Parts 1–6 answer the template items. Annexes A–D (blinding, pinned processing, deviations policy, timeline), H (rule → code crosswalk), I (registered implementation rules) and J (module status) are part of the registered plan. Annexes E–G hold attachments, references and the finaliser checklist.

---

## Part 1. Study information

### Q1. Title
**LEAD-χ: Is the midbrain susceptibility (iron) deficit of established psychosis present in help-seeking psychosis-spectrum youth? A procedurally blinded quantitative susceptibility mapping study of Penn LEAD (OpenNeuro ds007116 v1.0.6) with ADHD modelled factorially.**

### Q2. Authors
- Halil Aziz Velioglu [FILL: ORCID; affiliation: Feinstein Institutes for Medical Research, Northwell Health] et al. [FILL: co-authors, affiliations, ORCIDs].
- The first author is the sole analyst, the lock custodian and the single visual-QC rater (Annex A).
- **AI disclosure.** The protocol, the design simulations, this document and all hashed scripts were prepared with an AI system (Claude, Anthropic) under the first author's direction. The AI system will also run the hashed scripts under the rules in Annex A. It is not an author. Its session transcripts are kept as order-of-work logs.

### Q3. Research questions
- **RQ1 (primary, confirmatory).** The design is an effect-coded 2 × 2 of psychosis spectrum (PS) × ADHD. Averaged over ADHD status (the unweighted mean of the two ADHD-stratum simple effects), is SN–VTA χ lower in PS youth than in non-PS youth at each subject's first usable MEGRE session? The answer is one of three registered outcomes:
  - **A:** a deficit is present;
  - **B:** a group-level deficit of adult magnitude is excluded;
  - **C:** inconclusive.
- **RQ2 (key secondary, confirmatory; Holm).**
  - **S1:** Is the PS effect relatively nigral, and the ADHD effect relatively striatal?
  - **S2:** Does SN–VTA χ covary with PRIME psychosis-spectrum symptoms, adjusted for depression, irritability and ADHD?
  - **S3:** Is right dentate χ coupled to tissue-weighted neurite density of the right proximal superior cerebellar peduncle (pSCP_R)?
- **RQ3 (temporal stability; conditional on the final H1 outcome A).** In the follow-ups held out from group comparison until unlock (about 2 years later), did the PS deficit decline by at least 0.4 total SD?
- **RQ4 (estimation and exploratory; never claimed).** Developmental rates, cell means and the interaction, the ADHD striatal effect, TD-referenced deviations, a threshold-psychosis anchor, and seven named exploratory analyses (Q25).

### Q4. Hypotheses

**Rationale.**
- **Adult evidence.** In psychosis, midbrain susceptibility χ, a tissue-iron marker, is lower:
  - substantia nigra QSM meta-analytic g = −0.61 across mixed illness stages (Vano et al., 2026, *AJP*);
  - SN–VTA d = −0.66 in early-course schizophrenia, an adjusted difference of −14.24 ppb, with lower χ going with higher striatal [18F]-DOPA uptake (Vano et al., 2025).
- **The gap.** No iron or susceptibility study exists below first-episode psychosis (Pacchiano et al., 2026). The nearest youth study imaged neuromelanin, not iron (Tubiolo et al., 2026).
- **Why ADHD is modelled.** In Penn LEAD, 41 of 69 PS youths also have ADHD, and ADHD iron findings are mixed (Tang et al., 2022; Schulze et al., 2025; Morandini et al., 2026). The design is therefore factorial. The PS × ADHD interaction is **not testable** here (minimum detectable interaction ≈ 1.1 SD) and separability from ADHD is addressed only by the underpowered K1b label.
- **S3.** Human cerebello-VTA fibres arise mainly from the interposed nucleus and paravermal cortex (Hoffman et al., 2024), and the proximal SCP carries the dentate and interposed outflow. Oligodendrocytes need iron to myelinate (Vano et al., 2026, *Mol Psychiatry*). To our knowledge, dentate iron has never been related to outflow microstructure.
- **Directions.** All directions were fixed before any MEGRE was processed.

**Statistical framework.**
- Frequentist decision tests on OLS or mixed models; HC3 heteroscedasticity-consistent standard errors for every OLS decision test (t reference, df = n − p).
- Equivalence-type bound tests (Lakens, 2017); Holm family-wise control for the key secondary tests; intersection–union logic for claim labels (no extra α).
- Bayes factors, permutation p-values, bootstrap intervals and robust fits are reported but never decisive, except where a label is explicitly defined on them (K4 "robust"; S3b).
- Notation: δ = β/σ (residual-SD units), σ the residual SD of the fitted model. Negative means lower in PS. β_PS is the **unweighted** mean of the two ADHD-stratum simple effects (effect coding ±0.5), not a sample-weighted average.

**4.1 Primary hypothesis H1 and the decision rule.** Averaged over ADHD strata, PS youth have lower SN–VTA χ than non-PS youth: **β_PS < 0**. Benchmark: the adult magnitude, δ ≈ −0.6; null of interest: δ = 0.

*Step 1 — M1 outcome* (model M1, Q19; whole-brain reference):
1. **A, deficit present:** β̂_PS < 0 and two-sided HC3 p < .05.
2. **B, adult-magnitude group-level deficit ruled out:** the one-sided test of H0: δ ≤ −δ_B is rejected at α = .05, i.e. t_B = (β̂_PS + δ_B·σ_B)/SE_HC3 > t₀.₉₅,df, where
   - **δ_B = 0.6·√min(1, r̂_L)**, and r̂_L is the Spearman–Brown-transformed lower limit of the 95% Fisher-z CI of the split-echo correlation of χ_primary (PC4; Q24). This is a **heuristic, conservative shrinkage of an observed-scale adult bound**: δ_B ≤ 0.6 always, so false B against an observed −0.6 stays ≤ .05 (0.038 at the bound; Q21). Split-echo reliability is an **upper bound** on measurement reliability (the halves share acquisition, motion, background field, mask, registration, ROI and reference), and the adult benchmark is itself an observed, attenuated effect; δ_B is therefore not an attenuation correction. The B wording refers to the ppb bound δ_B·σ_B, which is always reported.
   - **σ_B = min(σ̂_resid, σ̂_rob)**, with σ̂_rob the Huber proposal-2 residual scale of M1 (`MASS::rlm`, psi.huber, k = 1.345, maxit = 200). If rlm fails or does not converge, σ̂_rob = 1.4826 × MAD of the OLS residuals, flagged in the report. This stops a few dominant scans from inflating the SD and making B easier (Paper 0 lesson iv).
   - **B is evaluable only if PC4 passes (split-echo R_SB ≥ 0.70) and r̂_L ≥ 0.25.** Otherwise B is not declared, the outcome is A or C, H1 is labelled **"low-precision measure"**, and the bound test is reported descriptively.
3. **C:** neither A nor B.

*Step 2 — final H1 outcome* = the M1 outcome after the K3 reference-agreement rule (§4.2): if K3 triggers, A requires A under both the whole-brain and the CSF reference, B requires B under both, and anything else is C. Every downstream rule uses the **final** outcome (T1 "conditional on A", the K4 "robust" label, the registered sentences), except the K4 fragility and influence refits, the ses-1-only and the exposure-sensitivity refits, which compare **M1 outcomes** (Step 1) with the M1 outcome of the full sample.

*Co-reports (never decisive):* B at the unadjusted bound (−0.6), at −0.5, at the bound scaled by the label-free cross-session reliability (lower 95% limit of the Pearson r between returners' first two usable sessions, a lower bound; gate JSON `rL_cross`), and, after unlock, at the bound scaled by the split-echo reliability of covariate-residualised χ.

*Disconfirmation:* β̂_PS > 0 with two-sided p < .05 is reported as **disconfirming H1**; the disconfirming sentence then replaces any B wording (the B test is reported numerically).

- **Registered wording of B** (§4.6): *"An adult-magnitude SN–VTA susceptibility deficit is not a feature of clinically ascertained, help-seeking psychosis-spectrum youth as a group. This does not exclude a deficit confined to the minority who will develop psychosis, nor a deficit that emerges later."*
- **Prohibited phrasings:** "iron loss emerges with illness", "progression account supported", "not a vulnerability marker".
- **Dilution caveat (registered).** If a trait deficit of 0.6 SD existed only in future converters, a fraction p = 0.15–0.30 of the PS group, the group-level δ would be only ≈ 0.09–0.18 (Salazar de Pablo et al., 2021; Linscott & van Os, 2013). Such a deficit can neither be detected nor excluded here, and the paper says so.

**4.2 Claim qualifiers K1–K4 for H1 (intersection–union; no extra α).**

| Qualifier | Rule |
|---|---|
| **K1a "PS effect"** | Final outcome A, **and** both simple effects (PS-only − TD; PS+ADHD − ADHD-only) have point estimates < 0 |
| **K1b "PS-specific (not ADHD)"** | K1a, **plus** a TOST on β_ADHD (SN–VTA) at **±0.5·σ_B** that succeeds (both one-sided tests at .05), **plus** a 95% CI of the PS × ADHD interaction that includes 0. The last is a **non-rejection** only; it is not evidence of an absent interaction. The label is underpowered by construction (0.36 at δ = −0.6; Q21) |
| **K1c "Concentrated in comorbid youth"** | Final outcome A; est(PS+ADHD − ADHD-only) < 0 with p < .05; PS-only − TD p ≥ .05 **and** est(PS-only − TD) > 0.5 × est(PS+ADHD − ADHD-only) (with the comorbid estimate negative, this means the PS-only estimate is less negative than half of it). When K1c applies, "PS signature" wording is prohibited |
| Precedence | K1c > K1b > K1a: one K1 label is reported |
| **K2 "lower iron"** | Refit M1 with R2\*_primary as the outcome. "Iron" wording is allowed only if β_PS(R2\*) < 0 with one-sided p < .05; otherwise "lower susceptibility; the iron-versus-myelin source is unresolved" |
| **K3 reference agreement** | Refit M1 on Δref = χ_primary(whole-brain reference) − χ_primary(CSF reference), over the sessions with a valid CSF reference (§12.2; sessions without one are missing, counts reported by cell). **Trigger (inferential):** \|β̂_PS(Δref)\| > max(0.25·\|β̂_PS\|, 0.10·σ̂_resid) **and** two-sided HC3 p(Δref) < .10. Under the trigger the conjunctive rule of §4.1 Step 2 applies |
| **K4 "fragile"** | Any single leave-one-subject-out refit changes the M1 A/B/C category (this covers both the A and the B decisions) |
| **K4 "influence-dependent"** | The top-1 subject carries > 20% of the residual sum of squares or any \|externally studentised residual\| exceeds 5, **and** refitting without those subjects changes the M1 A/B/C outcome (for S1–S3: a Holm decision) |
| **K4 "robust"** | Final outcome A, **and** ≥ 80% of the evaluable specifications of the QSM specification curve have δ < 0, **and** the joint Freedman–Lane permutation p of the median specification-wise HC3 t is < .05 (one-sided, 10,000 permutations; §24.3). Not assessable if fewer than 400 of the 486 specifications are evaluable |
| **Consequence of K4 labels** | "Fragile" and "influence-dependent" must appear in the abstract and in the first sentence of the Discussion whenever H1 is mentioned; the outcome itself does not change. "Robust" is the only positive robustness claim allowed |

**Always co-reported with H1:**
- the PS-only − TD simple effect, which compares the two cells without an ADHD diagnosis. It is **not** described as "unexposed": stimulant exposure is unmeasured, and one prescribed-stimulant endorser sits in each of the two cells;
- the ses-1-only analysis (n = 95);
- the fact that all 6 threshold-psychosis youths sit in PS+ADHD.

**4.3 Key secondary hypotheses** (Holm FWER .05 at .0167 / .025 / .05; all tested whatever H1 shows; a Holm-significant test is "supported" only in the predicted direction, otherwise it is reported as "opposite-direction significant").

- **S1 (topographic contrast).**
  - **L = (β_PS,SN–VTA − β_PS,STR) − (β_ADHD,SN–VTA − β_ADHD,STR) < 0**, estimated as β_PS(D) − β_ADHD(D), where D = z(χ_primary) − z(χ_STR) with the pooled label-free constants of the gate JSON.
  - The label **"dissociation"** also requires β_PS(D) < 0 **and** β_ADHD(D) > 0, each one-sided at .05. If L is significant but only one component holds, the wording is "the PS effect is relatively more nigral" or "the ADHD effect is relatively more striatal"; if neither, the contrast only.
- **S2 (dimensional).**
  - On each subject's **first usable session only** (no substitution by a later session), in youths aged ≥ 11 at that session with PRIME, BDI and ARI present: higher PRIME total goes with lower SN–VTA χ (**β_PRIME < 0**), adjusted for BDI, ARI and ADHD diagnosis.
  - **Registered wording:** *"SN–VTA χ covaries with PRIME, a symptom dimension that partly indexes PS status."* "Beyond diagnosis" is allowed only if the + PS_e sensitivity model also gives β_PRIME < 0 with p < .05.
- **S3 (iron–microstructure coupling).**
  - Right dentate χ is positively associated with pSCP_R tissue-weighted ICVF (**β > 0**), adjusted for age, sex, race, diagnosis, global deep-grey-matter χ (χ_DGM), locally referenced b0, |pitch|, tract length, motion, streamline count, dentate volume and Q_phase.
  - S3 is confirmatory only if gates PC4-DN, G2-DN and PC5 all pass; otherwise it is exploratory and the Holm family becomes S1 and S2 (.025 / .05). **D12 (confirmed): S3 stays in the Holm family** when its gates pass.
  - **Registered reading:** *"S3 is adequately powered only for strong latent coupling; a null S3 does not exclude moderate coupling."* Quantitatively (Q21): the latent correlation ρ is defined with **total** right-dentate χ; if dentate χ shares no variance with the model covariates, power is 0.63 / 0.70 / 0.80 at ρ = 0.5 and 0.38 / 0.46 / 0.58 at ρ = 0.4 across the three Holm levels; if 30–50% of its variance is shared with χ_DGM, power at ρ = 0.5 falls to 0.28–0.62. The label-free r(χ_DN_R, χ_DGM) is posted in amendment 2.
  - **Specificity label S3b** ("specific to the dentate→SCP pair") requires S3a plus both negative-control differences > 0 (one-sided 95% cell-stratified case-bootstrap intervals, 5,000 resamples): NC1, dentate vs right-putamen χ as the predictor; NC2, pSCP_R vs distance-matched right MCP as the outcome.
  - A positive S3 means "coupled", never "causes" or "mediates".

**4.4 Temporal stability T1** (evaluated only if the final H1 outcome is A; otherwise descriptive).
- **"Did not decline by ≥ Δ":** the one-sided test of H0: β(PS_e:occ) ≥ Δ is rejected at .05 (Wald t at the verified REML maximum with nlme containment df), with **Δ = 0.4 total SD fixed a priori** (total SD = the pooled label-free SD of χ_primary over the frozen first-usable sessions), **and** M1 refitted on the follow-up sessions alone gives β̂_PS < 0 (only its sign is used).
- **"Declined":** β(PS_e:occ) > 0 with two-sided containment-df p < .05, and the follow-up-only β̂_PS ≥ 0.
- **"Indeterminate":** otherwise.
- Δ is reported in ppb next to the occasion-0 estimate in total-SD units. At the adult δ = −0.6 residual SD the occasion-0 deficit is about 0.45–0.55 total SD, so Δ = 0.4 excludes only near-complete disappearance; if the occasion-0 deficit is ≤ Δ, complete disappearance lies inside the alternative and this is stated. The words "stable" and "persistent" are prohibited.
- The two-occasion LMM assumes that missing follow-ups are missing at random given the model covariates.
- T1 is within-person stability in the same people. It is **never** called a replication.

**4.5 Tiers.**

| Tier | Content | Error control | Claim allowed |
|---|---|---|---|
| 0 Gates | Sign check, PC1, PC2w, PC3, PC4, PC4-DN, PC5, G2, G2-DN, G3, G4, G5 (Q24) | Registered criteria, computed label- and age-free before unlock by hashed scripts | Pipeline valid; branch taken |
| 0b Convergent validity | ΔR2\*_STR vs R2\*_STR and χ_STR (only if D10 applies, Annex B.8) | Descriptive r with 95% CI | Validity of the striatal measure |
| 1 Primary | H1 final outcome (M1 → K3) | α .05: two-sided for A, one-sided for B | Confirmatory |
| 1b Qualifiers | K1, K2, K4 | Intersection–union; K4 "robust" by the joint permutation test | Wording only |
| 1c Stability | T1 | One-sided .05 against Δ, conditional on final A | Label only |
| 2 Key secondary | S1, S2, S3a | Holm (.0167 / .025 / .05); S3 conditional on its gates | Confirmatory per hypothesis |
| 2b Labels | "Dissociation"; S3b | Intersection–union; bootstrap intervals | Label only |
| 3 Estimation | E-dev, E-ADHD, E-int, E-z, E-anchor (Q25) | None (CIs; a calibrated p for E-dev) | Estimates only |
| 4 Exploratory | X1–X7; the node-level profile; the DWI multiverse | BH-FDR q < .05 within block; node-level max-T | Hypothesis-generating only |

- No tier-3 or tier-4 result appears in the abstract as a finding. No test moves between tiers after unlock. No study-wide FWER is claimed across tiers.
- **Non-confirmatory by construction:** any statement about DWI or phenotype alone. Those data were analysed before this registration (Q18).

**4.6 Registered sentences.** The confirmatory report prints only these sentences (verbatim; `R/m1.R`, `R/report.R`), selected by the registered rules; the paper may add numbers and citations but not change their claims.

| Condition | Sentence |
|---|---|
| Final A | "Outcome A: a lower SN-VTA susceptibility is present in clinically ascertained, help-seeking psychosis-spectrum youth as a group (averaged over ADHD status)." |
| Final A+B | "Outcomes A and B: a lower SN-VTA susceptibility is present in help-seeking psychosis-spectrum youth as a group, and it is smaller than the adult-magnitude bound." |
| Final B | the registered wording of B (§4.1), followed by the dilution caveat |
| Final C | "Outcome C: the data neither establish nor exclude an SN-VTA susceptibility deficit of the tested size in help-seeking psychosis-spectrum youth; estimates are offered for meta-analysis." (followed by the dilution caveat) |
| Disconfirming | "Disconfirming H1: SN-VTA susceptibility is higher, not lower, in help-seeking psychosis-spectrum youth (beta_PS > 0, two-sided p < .05)." |
| K1a | "The lower susceptibility is seen in psychosis-spectrum youth both with and without ADHD (point estimates)." |
| K1b | "The lower susceptibility is attributed to psychosis-spectrum status and not to ADHD: the ADHD effect is equivalent to zero within +/-0.5 SD and the interaction is not significant (non-rejection only)." |
| K1c | "The lower susceptibility is concentrated in psychosis-spectrum youth with comorbid ADHD; no 'psychosis-spectrum signature' is claimed." |
| K2 met / not met | "Lower R2* accompanies the lower susceptibility; 'lower iron' wording is used." / "Lower susceptibility; the iron-versus-myelin source is unresolved." |
| B not evaluable | "The SN-VTA measure did not reach the registered precision (PC4 or r_L below the floor); outcome B was not evaluated and H1 is labelled 'low-precision measure'." |
| K4 | "H1 is labelled 'fragile': deleting a single subject changes the A/B/C outcome." / "H1 is labelled 'influence-dependent': removing the flagged dominant subject(s) changes the A/B/C outcome." / "H1 is labelled 'robust' across the registered QSM specification curve." |
| T1 | "Within the same youths, the psychosis-spectrum SN-VTA susceptibility difference did not decline by 0.4 total SD or more over about two years. This is within-person stability, not a replication." / "… declined over about two years. This is within-person change, not a failed replication." / "The temporal-stability test (T1) was indeterminate." / "T1 was not evaluated because outcome A does not hold; follow-up estimates are descriptive." |
| S1 | "Dissociation: the PS effect is relatively nigral and the ADHD effect relatively striatal." / "The PS effect is relatively more nigral than striatal." / "The ADHD effect is relatively more striatal than nigral." / "The topographic contrast L is significant; neither component is significant on its own." / "No topographic dissociation was established (S1 not significant under Holm); S1 was powered only for large contrasts." |
| S2 | "SN-VTA chi covaries with PRIME, a symptom dimension that partly indexes PS status." (+ "The PRIME association holds beyond PS diagnosis." if allowed) / "No association between SN-VTA chi and PRIME was established (S2 not significant under Holm)." |
| S3 | "Right dentate chi is coupled to right proximal superior cerebellar peduncle tissue-weighted neurite density (coupled; not causal)." (+ "The coupling is specific to the dentate-to-SCP pair (both negative controls excluded).") / the registered reading (§4.3) / "S3 is exploratory because a dentate gate or PC5 failed." |
| Opposite direction | "A Holm-significant result in the direction opposite to the prediction is reported as such and not claimed." |
| Branch caveats | "Pipeline not validated (positive-control or registration gate failure): H1 is reported descriptively only."; "developmental positive control not confirmed"; "inversion-algorithm agreement not established"; "method-dependent ROI measure"; "QC reliability not established" (Q24) |

---

## Part 2. Data description

### Q5. Datasets
**Primary dataset.**
- Penn LEAD (Penn Longitudinal Executive functioning in Adolescent Development), OpenNeuro **ds007116 v1.0.6**, CC0 (Sevchik et al., 2026, *Sci Data*; Sevchik et al., 2025, bioRxiv; dataset doi:10.18112/openneuro.ds007116.v1.0.6).
- 225 imaging sessions from 132 youths aged 8–16 at enrolment (TD 36, ADHD 27, psychosis-spectrum "PRO/CHR" 69), with up to three sessions about 2 years apart, all on one Siemens Prisma 3T.

**Subsets used.**
- **MEGRE:** 181 complete four-echo magnitude + phase sessions from 113 youths, the only QSM-capable data. None has ever been reconstructed into χ.
  - Sequence: 3D monopolar spoiled GRE; TE 7.5/15/22.5/30 ms; TR 35 ms; FA 15°; 1.5 mm isotropic; 120 × 160 × 96 matrix (96 slab partitions).
  - Coil combination: Adaptive Combine with Prescan Normalize. Hardware: Prisma_fit, syngo MR E11, Head_32.
  - Protocol epochs: `QSM_SWI_1.5mm` early epoch (16 sessions) vs `qsm_acq-1.5mm_GRE`.
  - Sidecar `MagneticFieldStrength` is 3 in every session. ImagingFrequency takes six values in three clusters (123.166–123.262 MHz; B0 = 2.8928–2.8950 T).
  - Pixel bandwidth, phase-encode direction, parallel imaging, partial Fourier, slab oversampling, flow compensation and acquisition time are read from the sidecars in the header-only audit and reported (consensus reporting; Bilgic et al., 2024) [FILL: values from the header audit].
- **T1w and T2w** of the same sessions.
- **Phenotype:** `participants.tsv`, `sessions.tsv`, and the `phenotype/` files PRIME, BDI-1A child, ARI, Tanner (boy/girl), substance, the short Penn CPT and RSAS.
- **Existing derivatives made by the first author from ds007116.**
  - QSIPrep 26.0.0 anatomical outputs: ACPC T1w, SynthSeg `desc-aseg_dseg`, and T1w → MNI152NLin2009cAsym transforms (used only on the DWI side, §12.3).
  - pyAFQ 2.1 bundles and profiles: 217 DWI sessions; SCP_R, MCP_R and VTAuncross_R used.
  - QSIRecon 26.0.0 AMICO NODDI maps (the container file is named `qsirecon-0.23.0.sif`, but its label and `qsirecon.toml` give version 26.0.0).
  - Paper 0's FW-DTI profiles.
- **Conditional upstream derivatives** (download approval A3; Annex B): ds007088 fMRIPrep 25.0.0 transforms, masks and confounds for ΔR2\* (D10); ds006744 ASLPrep CBF tables (X7); ds007402 XCP-D rest time series (X5); ds007091 AutoTrack scalarstats (optional S3 cross-pipeline check).

**Companion study (completed, DWI only).** Paper 0 (Q17) analysed the same cohort's diffusion data. No MEGRE value was used in Paper 0.

### Q6. Data availability
**Publicly available, no access restrictions** (CC0; OpenNeuro). The upstream derivatives above are also public on OpenNeuro.

### Q7. Data access
- https://openneuro.org/datasets/ds007116/versions/1.0.6 ; doi:10.18112/openneuro.ds007116.v1.0.6.
- The dataset asks to be acknowledged through its `HowToAcknowledge` text (the bioRxiv preprint) and its DOI. The *Sci Data* descriptor (doi:10.1038/s41597-026-08095-1) is also cited.

### Q8. Access dates

| Who | Item | Date |
|---|---|---|
| First author | Full ds007116 v1.0.6 fetched from `s3://openneuro.org/ds007116/`, verified byte-for-byte against OpenNeuro git tag 1.0.6 (CHANGES "1.0.6 2026-07-28"); this lab-server copy (read-only) is the analysed copy | 2026-08-25 |
| First author | First access to DWI data (own tractometry) | on or before 2026-09-06 |
| First author (AI agents working for him) | First access to MEGRE metadata (file listings, JSON sidecars, NIfTI headers; E4 rows M1–M4) | 2026-09-23 |
| First author (AI agents working for him) | First access to MEGRE voxel data (feasibility QC; E4 rows M5–M6) | 2026-09-23 |
| First author | Upstream derivatives: S3 object listings only, plus one ds007091 NODDI scalarstats table (DWI) for a format check | 2026-09-23 |
| Co-authors | [FILL: per co-author, date of first access to metadata and to data, or "none"] | |

Knowledge gained from the dataset's own documents (descriptor and preprint: recruitment, group definitions, cohort characteristics) is part of Q18.

### Q9. Data collection procedures
- Described in Sevchik et al. (2026, *Sci Data*) and the preprint. Youths were recruited from clinics of the Children's Hospital of Philadelphia and assessed with a modified K-SADS and the Structured Interview for Psychosis-Risk Syndromes.
- Group assignment follows the dataset: PS = `study_group == "PRO/CHR"` (a subthreshold risk syndrome or a threshold psychotic disorder); ADHD = `dx_adhd == 1`. Diagnoses were made at baseline only.
- The MEGRE was released "not curated or preprocessed".

**Ethics and data use.**
- **IRB (approval A12).** The first author requests a Northwell Health IRB determination for this secondary analysis of de-identified public data (expected: not human-subjects research, or exempt) [FILL: date requested; determination number, category and date if received before posting]. **No QSM reconstruction, the pilot included, starts before the determination is on file.** If it arrives after posting, it is recorded in amendment 1.
- The data are CC0; the dataset's acknowledgement requests are followed.
- **No attempt is made to identify any participant as a person.** "Re-identification" in this document means only linking a coded analysis copy back to its public pseudonymous dataset ID (Annex A).
- **Derived-data sharing.** No subject-level derived data are posted before publication. After acceptance, blinded-QC'd χ and R2\* maps are offered back to OpenNeuro, with the Penn LEAD authors' agreement.
- Medication, SIPS item and follow-up diagnosis records do not exist in the release. Obtaining them would need a data-use agreement and probably an IRB amendment (A10, optional). Any analysis of such records is labelled post hoc if run after unlock.

### Q10. Codebook and documentation
- **Documentation:** `participants.json`, the phenotype JSON sidecars, the MEGRE JSON sidecars, and the PennLINC QC tables `EF_desc-T1_qc.csv` and `EF_desc-diffusion_qc.csv`.
- **Variables used:**
  - `participants.tsv`: study_group, dx_adhd, dx_psychosis (threshold psychosis), dx_mdd, dx_moodnos, sex, race, ethnicity;
  - `sub-*_sessions.tsv`: age at session;
  - `prime.tsv`: `prime_total_score` and items `prime_1` … `prime_12` (duration items not used); BDI-1A child total; ARI total;
  - Tanner self-report items (§12.5); `substance_othr_010` and `_011` (§12.5);
  - `short_penn_continuous_performance_test.tsv:cnb_cptnl_total_fp` (S2 sensitivity); `rsas.tsv:rsas_total_score` (X6).
- **Data traps respected:** E-SWAN is never used (higher scores mean strength); the PRIME flag is miscoded and n/a below age 11, so the total score is used; Tanner code 5 ("don't know") → NA for the PDS items (for the picture items 5 is a valid stage; they are not used); ALES zeros where not administered → NA (ALES is not used); the stimulant item is not asked at age ≤ 11, is missing in 61% of ses-1 rows and does not ask about ADHD medication, so it enters sensitivity analyses only and is labelled uninformative about exposure; CNB valid_code "0" is not a failure flag.

---

## Part 3. Variables

### Q11. Manipulated variables
None; the study is observational. The PS × ADHD cells are measured groupings, coded PS_e = +0.5 if study_group == "PRO/CHR", else −0.5; ADHD_e = +0.5 if dx_adhd == 1, else −0.5; their product is the interaction. Cells: TD, ADHD-only, PS-only, PS+ADHD.

### Q12. Measured variables

**12.1 Primary outcome: χ_SNVTA (ppb).** The probability-weighted mean susceptibility of the bilateral SN–VTA complex in the subject's first usable MEGRE session:

χ_SNVTA = Σ_v w_v χ_v / Σ_v w_v, with w_v = min(1, p_SNc + p_SNr + p_VTA + p_PBP).

- **Atlas:** CIT168 v1.1.0 `CIT168toMNI152-2009c_prob.nii.gz` (MNI152NLin2009cAsym, 1 mm; 0-based volumes SNc 6, SNr 8, PBP 9, VTA 10; Pauli et al., 2018), warped to native MEGRE space by the composed inverse transforms (linear interpolation of the probabilities; one interpolation).
- The summed probability is **clipped at 1**, then thresholded at **w ≥ 0.5**. Hemispheres are pooled.
- **Support S** of each reconstruction = Mask 3 ∩ the reconstruction's own final (background-removal/inversion) mask ∩ {finite, non-zero χ}. S excludes the rim that V-SHARP erosion and MEDI leave as exact zeros inside Mask 3, whose share varies with brain size and background-removal method.
- **Voxels included:** those in S, more than 1 voxel (1.5 mm) from SynthSeg CSF (labels 4, 5, 14, 15, 43, 44, plus 24 where present), and **not implausible** (§12.4). No χ-based and no vesselness voxel selection is applied. χ is never resampled.
- **Reference (two-step, on S).** ref₀ = mean χ over S; implausible voxels are those in S with |χ − ref₀| > 1.0 ppm; the reference is the mean over S minus the implausible voxels, and χ is reported relative to it (whole-brain reference = 0; Bilgic et al., 2024). For G3 and the specification curve, the reference and the ROI means of all 9 background-removal × inversion nodes are computed on the per-session **intersection** of their supports, so the nodes differ only by algorithm.
- **Units:** ppb (ppm × 1,000).
- **B0.** The field-to-ppm conversion uses each session's centre frequency f0 (sidecar ImagingFrequency, MHz). Each engine is given the field strength that makes γ_engine·B0 = f0: QSMxT `MagneticFieldStrength` = round(f0/42.576, **4**) T in the coded sidecar; SEPIA CF = f0 in Hz. Using the nominal 3 T would underestimate χ by 3.5%.
- **Pipeline** (Annex B): phase-offset removal (MCPC-3D-S type) → ROMEO in template mode → weighted echo combination → V-SHARP → MEDI (CSF term off), **single pass on Mask 3**.
- Under the G3 fallback the primary ROI is the SN-only complex (CIT168 SNc + SNr) with the same rules; "χ_primary" denotes whichever ROI G3 selects.

**12.2 Other imaging variables.**

| Variable | Definition | Use |
|---|---|---|
| R2\*_ROI (s⁻¹) | ARLO fit on the 4 magnitude echoes (QSMxT, primary run), clipped voxelwise to [0, 200] s⁻¹; same voxel sets as the matching χ ROI, restricted to Mask 3 | K2, PC3, X2 |
| χ_STR | Mean over SynthSeg caudate (11/50) + putamen (12/51), each eroded 1 voxel, plus accumbens (26/58), not eroded | S1 |
| χ_DGM | Caudate + pallidum (13/52), eroded 1 voxel | S3 covariate |
| χ_DN_R | Probability-weighted mean over the SUIT right dentate (label 30, probseg p ≥ 0.5; Diedrichsen et al., 2011; NLin6AsymC mapped to 2009cAsym by a template-to-template SyN, Annex B.6) | S3; X1 |
| χ_DN+IP_R; χ_DN_L | SUIT right dentate + interposed (labels 30 + 32); left dentate (label 29), p ≥ 0.5 | S3 sensitivities |
| χ_PUT_R | SynthSeg right putamen, eroded | NC1 |
| χ of GP, thalamus (10/49), caudate, putamen, RN (CIT168 vol 7), STN (vol 15) | Eroded SynthSeg, or CIT168 p ≥ 0.5 | PC1, PC2w, X1 |
| SN-only; VTA+PBP | CIT168 SNc+SNr; VTA+PBP | G3 fallback; specification curve; descriptive |
| Δref | χ_primary(whole-brain reference) − χ_primary(CSF reference). CSF reference: lateral ventricles (4/43) eroded 2 voxels (1-voxel erosion if fewer than 200 voxels remain), within S, implausible voxels excluded; SynthSeg has no choroid-plexus label, so ventricular voxels with R2\* > ventricle median + 2 SD are dropped; ≥ 200 voxels and \|mean CSF − whole-brain mean\| ≤ 50 ppb are required. **Otherwise the session has no valid CSF reference** and is missing for Δref, K3 and the CSF level of the specification curve (counts by cell reported). The PLIC is not used as a reference | K3 |
| roi_vol (mm³) | Σ w_v × 3.375 mm³ over the atlas voxels with w ≥ 0.5 of the warped CIT168 SN–VTA atlas (atlas only: no mask, no CSF buffer; T1-derived, χ-independent) | M1 covariate |
| str_vol, DN_R_vol (mm³) | Native volumes of the matching ROI voxel sets | S1, S3 covariates |
| Q_phase | **Phase-noise/linearity index** (not a motion index): log of the median absolute multi-echo phase-prediction residual inside Mask 3 (echo 3 predicted as 2φ₂ − φ₁, echo 4 as 3φ₂ − 2φ₁, residuals wrapped to (−π, π]); free of φ₀ and of any static field | Covariate |
| motion_gsr | **Magnitude motion index:** ghost-to-signal ratio on the echo-1 magnitude: head = voxels above the Otsu threshold (holes filled); ghost band = the head shifted by half the FOV along the phase-encode axis, outside the head dilated 3 voxels; GSR = (mean ghost − mean remaining background)/mean head. Label-free | H1 sensitivity; tabulated by cell |
| abs_pitch (°) | Head pitch relative to B0 from the rotation (polar decomposition) of the composed MEGRE→T1w rigid and the affine stage of the T1w→MNI152NLin2009cAsym registration, **one frame for every session**; roll, yaw and total tilt also recorded | S3 covariate; H1 sensitivity; specification curve |
| ref_frac_wm / gm / csf; fov_brain_frac | SynthSeg tissue fractions of the reference domain (S minus implausible voxels); share of the T1w SynthSeg brain inside the MEGRE field of view | H1 sensitivity; tabulated by cell |
| χ_primary (harmonised reference) | χ_primary referenced to the mean over S restricted to the MNI-defined region inside the MEGRE field of view in ≥ 95% of sessions (computed label-free from the FOV masks) | H1 sensitivity |
| χ_primary (split-echo) | From the {1,3} and {2,4} echo-subset reconstructions | PC4, B co-report |
| f0 cluster | Rank (1–3) of the session's ImagingFrequency cluster (distinct values split at gaps > 0.01 MHz); kept in the encrypted key | H1 sensitivity |
| Protocol epoch | `QSM_SWI_1.5mm` ("early") vs later; kept in the encrypted key | H1 sensitivity |
| Code D | Rubric failure code D (dental or metal artefact), recorded for every session in which it is visible | H1 sensitivity; tabulated by cell |

**12.3 DWI variables for S3.** All come from existing DWI outputs, are keyed by real ID, are computed by the hashed `dwi_geometry.py` and are frozen before unlock (Annex A). DWI is non-blind by construction.
- **pSCP_R tract.** The pyAFQ 2.1 SCP_R bundle; each streamline is oriented so that node 0 lies at the dentate end and resampled to 100 nodes; node centroids are the mean resampled points.
- **Candidate nodes 0–30.** The range is **data-informed but label-free**: it was set from template geometry and the pooled label-free SCP_R length and ISOVF profiles (E2 §4.3); round-1 SCP_R nodes 0–16 had been compared by group (Q18, D2). Per subject (χ-free, label-free) nodes are removed by: (i) upper limit min(30, n_int − 1), n_int = the node nearest the subject's SCP_R intermediate-ROI centroid; (ii) nodes whose centroid lies < 4 mm from the nearest right-dentate point (SUIT p ≥ 0.5 voxel centres mapped to ACPC space; 3 mm plus the 1-voxel dilation at the 1 mm atlas grid); (iii) nodes with ISOVF > 0.30; (iv) for FW-DTI metrics only (DWI multiverse), nodes with tissue weight < 0.2.
- **n_min.** The largest value in {10, 8, 6} that keeps ≥ 90% of the 108 pre-QC S3 candidate sessions (attachment `s3_candidates_preqc.csv`, from metadata) after rules (i)–(iii); computed by `dwi_geometry.py` before posting. **n_min = [FILL: value from dwi/n_min.json]**. Sessions below n_min leave S3.
- **twICVF_pSCP_R** (primary) = the mean over retained nodes of profile(ICVF × TF)/profile(TF), TF = 1 − ISOVF (Parker et al., 2021), with profiles by Gaussian-weighted averaging along the bundle (dipy `afq_profile`, as pyAFQ).
- **Covariates:** `b0loc_pSCP_R` = node-wise N4-corrected b = 0 signal divided by the median N4 b0 of the right cerebellar WM (SynthSeg 46, eroded 1 voxel, excluding voxels within 2 voxels of the mapped dentate), averaged over the retained nodes; `tract_len` (mean SCP_R streamline length, mm); `nSL` (streamline count); `FD`, `NDC` (PennLINC DWI QC).
- **NC2 outcome.** Right-MCP tw-ICVF over MCP nodes 0–79 whose dentate distance lies within ±1 mm of a retained pSCP_R node's (terminal nodes 80–99 not eligible).
- **Sensitivity variables:** VTAuncross_R nodes 0–30 under rules (ii)–(iii); peak-aligned tw-ICVF (the session profile shifted by s ∈ [−5, 5] nodes to maximise its Pearson r with the pooled mean SCP_R profile over nodes 0–30, then the retained-node mean); calibrated T1w/T2w along pSCP_R (T2w rigidly registered to the ACPC T1w; ratio scaled so that the median in eroded lateral-ventricle CSF is 0 and in eroded cerebral WM is 1); split-half streamline reliability (20 random halvings, pooled Spearman–Brown r; reported).
- **Physically constrained FW-DTI rule** (DWI multiverse only): voxels with MD_t = 0 and FW < 0.999, or with MD_t > 3.0 × 10⁻³ mm²/s, are failed fits, excluded from the tissue weights and from FW.

**12.4 Physical-plausibility rules** (pre-specified and label-free; applied by `qsm_plausibility.py` on coded data before unlock; Paper 0 lessons i and iv).
- **Post-run mask assertion.** A reconstruction is invalid for a session ("mask_mismatch") unless ≥ 99% of its final mask lies inside Mask 3 and it covers ≥ 90% of Mask 3 eroded by 2 voxels (QSMxT falls back silently to its own mask if the custom mask is not found).
- **Voxel.** |χ − ref₀| > **1.0 ppm** within S (ref₀ as in §12.1) is a failed fit; excluded from ROI means and from the reference, retained in a sensitivity analysis. Rationale: venous blood lies about 0.45–0.7 ppm above tissue (Jain et al., 2012; approximate) and youth deep-grey-matter iron is below 0.15 ppm.
- **ROI failure** (per session and reconstruction): implausible voxels carry > 5% of the ROI's in-support weight; **or** < 95% of its weight lies in S (< 90% for the dentate); **or** its mean lies outside the physical envelope.
- **Envelopes (ppb):** GP [−50, 300]; putamen [−80, 200]; caudate [−80, 200]; SN–VTA [−50, 300]; RN [−50, 300]; dentate [−50, 300]; thalamus [−100, 100]; composites take their components' envelopes (STR [−80, 200]; DGM [−80, 300]; right putamen [−80, 200]; dentate + interposed [−50, 300]). Set from Li et al. (2014), Zhang et al. (2018) and Treit et al. (2021); far wider than any tested effect.
- **R2\* envelope.** Every deep-grey-matter ROI mean must lie in [5, 100] s⁻¹.
- **Consequences:** a primary-ROI failure fails the session under first-usable criterion (c′) (Q13); any other failed ROI is missing for that session; failure of ≥ 2 of {GP, SN–VTA, RN, DN} grades the session 2.
- **Session indices:** f_imp = the implausible fraction in S eroded by 3 voxels (> 1% is an automatic fail); each of f_imp, χ SD in eroded WM, χ SD in the eroded ventricles, the median of |local field| in the 5-mm inner shell of S, the posterior-fossa streaking ratio SD(χ, cerebellar WM)/SD(χ, cerebral WM), and Q_phase forces a visual review when it exceeds the upper quartile + 3·IQR of the pooled coded distribution.
- The same rules apply to the split-echo, engine and multiverse reconstructions; a failure makes the value missing for that variant only.

**12.5 Covariates and phenotypes.**
- **age:** age at the MEGRE session, entered as `splines::ns(age, df = 3)` (interior knots at the tertiles).
- **sex_e:** ±0.5 (male +0.5). **race3** ∈ {Black, White, Other}, deviation-coded (`contr.sum`): "Black or African American" → Black; "White" → White; "Asian", "More than one race" and "Unknown or not reported" → Other; treated as a proxy for unmeasured structural and nutritional factors.
- **S2 predictors** (same session as the MEGRE): PRIME total (Kobayashi et al., 2008), BDI-1A child total, ARI total.
- **Sensitivity covariates:**
  - years since the diagnostic assessment (age at MEGRE − age at ses-1); Hispanic ethnicity;
  - **PDS-like puberty score:** mean of the available same-session Tanner self-report items with "don't know" (5) → NA; boys items 1–5 (height, skin, facial hair, body hair, voice); girls items 1, 2, 5, 6 (height, skin, body hair, breasts) plus menarche (item 3: yes → 4, no → 1); at least 3 items required; **menarche_e** = ±0.5 in girls, 0 in boys;
  - mood-diagnosis flag (dx_mdd or dx_moodnos); intracranial volume (SynthSeg, coded side; in place of eTIV);
  - **prescribed-stimulant endorser:** same-session `substance_othr_010` ∈ {1, 2} (ever used stimulants) and `substance_othr_011` = 1 (prescribed);
  - DWI mean FD (same session); registration Dice (the minimum of the four G5 Dice values); |pitch|; the motion index; the reference-domain fractions; f0 cluster; protocol epoch.

**12.6 QC variables:** the blinded QC grade (0/1/2) and failure codes (including code D), SN–VTA placement pass/fail, G5 Dice pass/fail, the plausibility flags, the post-run mask assertion, and the S3 node-count flag.

### Q13. Units of analysis, inclusion and exclusion, expected N
**Unit:** the subject, analysed at one session (the first usable session) for H1, S1, S2 and S3. T1 uses one follow-up per returner (occasion 1 = the earliest usable session after the first usable one).

**First usable session:** the earliest session meeting all of the following, frozen by coded ID before unlock (`build_usable.py`):
- (a) a complete 4-echo magnitude + phase MEGRE;
- (b) a same-session T1w not rated "fail" by PennLINC QC (sub-21568 ses-1 has no T1-QC entry; it is kept, flagged and dropped in a sensitivity analysis);
- (c) a blinded MEGRE QC grade ≤ 1 (and no automatic grade 2);
- (c′) primary-ROI physical plausibility (§12.4) and a valid reconstruction (post-run mask assertion);
- (d) SN–VTA ROI placement rated "pass" on T1w and echo-1 magnitude overlays, **never on χ**;
- (e) both χ-independent registration checks passed (G5).

If a session fails any of (b)–(e) or (c′), the next session of that subject is used.

**S2 session rule:** the first usable session only; if the subject is < 11 years old or lacks PRIME, BDI or ARI at that session, the subject leaves S2 (no substitution).
**S3 session rule:** the H1 session (the first usable MEGRE session, whatever its wave) with a same-session DWI: PennLINC scalar pass, an SCP_R bundle, a coded bundle-QC pass (rubric in `qc_rubric.md` §6; single rater, the first author, **not blind in practice**) and ≥ n_min retained nodes. No substitution.

**Expected N (pre-MEGRE-QC; counts from public metadata).**

| Sample | TD | ADHD-only | PS-only | PS+ADHD | Total |
|---|---|---|---|---|---|
| MEGRE subjects (181 sessions) | 30 | 22 | 24 | 37 | 113 |
| **H1 / S1 first usable** (after the T1w rule: 175 sessions / 112 subjects) | 30 | 21 | 24 | 37 | **112** |
| of which from ses-1 / ses-2 / ses-3 | 25/4/1 | 17/4/0 | 22/2/0 | 31/5/1 | 95/15/2 |
| H1 projected after MEGRE QC and plausibility (5–10% loss) | 27–29 | 19–20 | 22–23 | 33–35 | ≈ 100–106 |
| **S2** (age ≥ 11; PRIME, BDI and ARI present) | 19 | 12 | 24 | 28 | **83** |
| **S3** pre-QC candidates (H1 session with a same-session DWI scalar pass and an SCP_R profile) | 30 | 20 | 23 | 35 | **108** |
| S3 conservative (dropping the 20-subject round-1 AFQ-QC list; sensitivity) | 26 | 16 | 19 | 33 | 94 |
| S3 projected after coded bundle QC and n_min | | | | | ≈ 80–97 |
| **T1** returners (first follow-up after the first usable session) | 17 | 12 | 8 | 16 | **53** |
| PC2w pairs (label-free; returners with 2 usable sessions) | | | | | 53 (≈ 48 after QC) |

- **Ages at the analysed session:** 8.8–17.8 y. r(age, PS) = 0.27: PS youths are about 1.2 y older.
- **Details:** 6 threshold-psychosis youths, all PS+ADHD, 2 of them analysed at ses-2; 13 early-epoch sessions (5/2/1/5); 5 prescribed-stimulant endorsers (1/2/1/1); 17 first-usable sessions are ses-2 or ses-3.
- **S3 exclusions:** sub-19861 and sub-20866 (no ses-1 DWI); sub-20916 and sub-21112 (ses-1 DWI scalar fail).

**Sensitivity exclusions** (each run separately, never replacing the primary): the 6 threshold-psychosis youths; 4 TD youths with a dx-coding gap or prescribed stimulants (sub-20139, sub-20252, sub-20303, sub-20953); the prescribed-stimulant endorsers; T1w "artifact" sessions; sub-21568; early-epoch sessions (also entered as a covariate); baseline age < 11; sessions with rubric code D; strict QC (grade 0); **automated exclusions only** (plausibility, G5 and T1w QC; the visual grade and placement rating ignored); without the (c′) plausibility exclusions; the **exposure-sensitivity set** (Q18).

**QC-by-group (after unlock, before any outcome), at the subject level.** Two tables by cell with Fisher's exact test: (1) subjects with ≥ 1 MEGRE session who contribute no first-usable session; (2) subjects whose earliest MEGRE session is not usable. If (1) gives p < .05, an inverse-probability-weighted M1 is co-reported with weights 1/P(in H1 | PS_e × ADHD_e, age at the earliest MEGRE session, sex_e) from a subject-level logistic model. Orientation angles, the motion index, ROI in-support fractions, slab coverage, reference-domain fractions, code D and CSF validity are tabulated by cell.

### Q14. Missing data
- **Outcomes:** complete case. A session missing χ_primary moves the first-usable rule to the next session; if no session remains, the subject leaves H1. Sessions without a valid CSF reference are missing for K3 only.
- **Differential loss** is examined by the subject-level QC-by-group tests (Q13).
- **Covariates:** roi_vol, Q_phase and the orientation angles exist for every reconstructed session. The one "Unknown or not reported" race youth is Other in the primary and set to missing in a sensitivity analysis.
- **Phenotype availability** at the first usable session (pre-QC, public metadata):

| Cell (n) | PRIME | age ≥ 11 | BDI | ARI | PDS computable | CPT |
|---|---|---|---|---|---|---|
| TD (30) | 27 | 22 | 28 | 30 | 18 | 27 |
| ADHD-only (21) | 20 | 12 | 19 | 21 | 13 | 19 |
| PS-only (24) | 24 | 24 | 24 | 24 | 20 | 23 |
| PS+ADHD (37) | 36 | 31 | 35 | 36 | 21 | 34 |

- **S2:** complete case (n = 83 of 86 with PRIME at age ≥ 11).
- **Puberty sensitivity:** `mice` is not installed and is not used; the registered branch is (i) complete case and (ii) one seed-fixed stochastic regression imputation of the PDS from age, sex_e, PS_e, ADHD_e and the subject's nearest other-session PDS with its age gap (when ≥ 20 observed pairs exist). The observed (pre-imputation) scores are reported by cell.
- **Attrition (E-dev, Tier 3):** IPW-weighted GEE (statsmodels); pattern-mixture returner × PS terms; a δ-adjusted tipping point (δ ∈ {−0.5, −0.25, 0, 0.25, 0.5} SD) by single stochastic imputation (the multilevel MI is dropped: `mice`/`jomo` are unavailable).
- **Missing ROIs:** an ROI failing §12.4 is missing for that session only. PC1 uses sessions with all its ROIs plausible.

### Q15. Outliers
- **No subject is excluded for an extreme outcome value.** Physically implausible values are handled only by the label-free rules of §12.4, which remove failed fits, not tissue values.
- **Robust flags (before unlock).** A label- and age-free robust z (median/MAD over plausibility-passing sessions) is computed per ROI; |z| > 5 is flagged and **retained**; only counts are visible.
- **Dominance diagnostics (at unlock, for M1, S1, S2, S3).** The share of the residual sum of squares carried by the top 1 and top 2 subjects, the maximum |externally studentised residual|, and the maximum Cook's D. If the top-1 share exceeds 0.20 or any |t_i| exceeds 5, the model is refitted without those subjects; a changed M1 outcome or Holm decision earns **"influence-dependent"** (K4). A Huber M-estimate is always co-reported.
- **Scale rules:** outcome B and the ADHD TOST use σ_B (§4.1). For S1 and T1, if SD/(1.4826·MAD) > 1.5 for either ROI, robust-scaled values are co-reported and the decision uses the registered scale. The DWI robust z (|z| > 5) is flagged and retained.
- **Expected final N:** as in Q13.

### Q16. Sampling weights
None. Inverse-probability weights enter only the conditional QC-by-group sensitivity and the attrition battery.

---

## Part 4. Knowledge of data

### Q17. Prior publications and dissemination
- **The first author's FEP paper** (Velioglu et al., 2026, *BP:CNNI*) used a different cohort. It is the source of the cerebello-VTA diffusion hypothesis that Paper 0 tested.
- **Paper 0**, a completed companion (Penn LEAD DWI only; manuscript draft v4, 2026-09-25, not yet public):
  - it asked whether PS youth show the lower cerebello-VTA MD reported in FEP;
  - all registered tests were null: DKI-MD d = −0.07 (90% CI −0.42 to 0.29); tissue MD_t d = −0.07 (−0.37 to 0.22); FW d = 0.03 (95% CI −0.37 to 0.44); every Holm-adjusted p = 1.0; no node survived family-wise correction;
  - its plan and amendments carry SHA-256 hashes (Annex E) [UNVERIFIED as a publication: no DOI yet; FILL: preprint/OSF DOI if posted first, D14].
- The first author's 2026-09 ses-1 tractometry of Penn LEAD is unpublished [FILL: any talk, poster or abstract, or "none"].
- **No MEGRE-derived result has been disseminated by anyone on this team.** The only MEGRE-derived quantities ever produced are the feasibility QC values in Q18. To our knowledge no published work has processed the Penn LEAD MEGRE (searches of 2026-09-23/24; re-run before submission).
- **Co-authors:** [FILL: per co-author, prior publications or presentations using these data, or "none"].

### Q18. Prior knowledge of the data
**Category.** "Registration prior to analysis of the data": the data have been accessed, but no analysis related to the research plan has been performed. Because the first author's viewing status is "unsure", the stricter categories cannot be claimed.

**Why the feasibility statistics are not "analysis related to the research plan."** Summary statistics were computed on MEGRE voxels before registration: the magnitude 99th percentile and phase ranges of all 181 sessions (M6), and whole-head voxelwise R2\* percentiles in 4 sessions (M8, M10). They were technical QC checks (phase scaling, monopolar readout, slab coverage), computed over the whole head or whole slab. None was ROI-based, none concerned the SN–VTA or the dentate, none was computed or displayed by group, and none is an input to any registered model. The template's category excludes "calculation of summary statistics" related to the research question; these are unrelated to it, and they are declared in full below.

**Exposure table** (audit of 2026-09-25 over all scripts, reports and 59 session transcripts, attachment E4; plus the P5 incident found in task F).

| # | Modality | What was accessed or computed | Sessions | Linked to group? |
|---|---|---|---|---|
| M1–M3 | MEGRE metadata | File listings; JSON sidecar fields (TE, TR, FA, software, coil combination, ProtocolName and others); NIfTI headers (dimensions, voxel size) | All 181 | **Availability and header metadata only**, cross-tabulated by group/cell. No image-derived value |
| M4 | MEGRE JSON | Full sidecars printed | 5 sessions | No |
| M5 | MEGRE voxels (in memory) | Echo-1/4 magnitude and phase min/max and unique values (scaling check) | sub-21432 ses-1, sub-20976 ses-2, sub-20188 ses-3 | No |
| M6 | MEGRE voxels (in memory) | Echo-1 magnitude and echo-1/4 phase: slab z-range, edge-slice tissue fraction, magnitude p99, phase ranges → `megre_coverage.tsv` (181 × 11). No script ever read this file; 2–4 rows were displayed three times | All 181 | **No** |
| M7 | MEGRE geometry + T1w | Slab coverage of FreeSurfer structures; for MEGRE, header affine only | 12 sessions | No |
| M8 | MEGRE voxels (in memory, server) | Odd/even-echo phase check; **whole-head voxelwise R2\*** in a crude scalp-inclusive mask (median/IQR only); a crude echo-difference frequency p95 | Same 3 as M5 | No join. sub-20188's group (ADHD) is named in an unrelated report |
| M9–M10 | MEGRE voxels (local copy) | 8 NIfTI copied to the Mac. Phase-range, φ₀-smoothness, echo-linearity and residue checks; **whole-head R2\*** percentiles | sub-20139 ses-1 | No join. Its group (TD) is named in an unrelated report |
| M11 | MEGRE images | Two figures (magnitude echoes 1/4, echo-4 phase, φ₀, coverage slices) were rendered and **viewed by an AI agent** and copied to the first author's report folder. File metadata show both copies were read on the Mac on 2026-09-24 at 17:29 UTC by a non-agent process, consistent with a Finder/Quick Look preview, but not proof | sub-20139 ses-1 | No |
| M12–M13 | MEGRE | Two design-phase re-loads of the M9 phase arrays (min/max; header fields for scrubbing); display of 2–4 rows of `megre_coverage.tsv` | As stated | No |
| M14 | MEGRE | **χ map, local field, background removal, dipole inversion, χ-separation, any ROI χ or R2\*, any MEGRE value by group** | **None, ever** | — |
| M15 | Audit | Imperfect numeric redaction during the audit showed a few whole-head QC numbers (magnitude p99 of 4 coverage rows, upper IQR bounds of the 3 M8 whole-head R2\* fits, M10 residue counts) | As stated | No |
| **P5** | Synthetic-code test (2026-09-25) | While the analysis code was tested on synthetic outcomes simulated on the real design, the console showed that 2 of the 4 exposure sessions (sub-20139 ses-1, sub-21432 ses-1) are in the synthetic first-usable sample; from the cell counts their public cells could be inferred (TD and PS+ADHD). This makes sub-21432's whole-head R2\* IQR bound (M15) group-linkable. No SN–VTA, dentate or ROI value was involved; synthetic runs now mask these details | 2 sessions | **Yes (public labels of 2 sessions)** |
| D1–D4 | DWI | **Fully non-blind at group level:** the first author's ses-1 tractometry (2026-09-06 to 09-13), round-1 forensics (PS × ADHD cell decomposition, ses-2 stability, brain–behaviour), and Paper 0 (8 bundles × 12 metrics incl. SCP_R/MCP_R). Round 1 reported SCP_R nodes 0–16 ICVF PS − TD t = −1.88 (p = .063). The first author has also seen the 20-subject round-1 AFQ-QC list | All DWI | **Yes** |
| P1–P4 | Phenotype | Group means and effect sizes of 77 summary variables (incl. PRIME, BDI, ARI, E-SWAN); stimulant items and PRIME by cell; Paper 0 Table 1; associations of PRIME/BDI/E-SWAN/2-back d′ with DWI | All | **Yes** |
| O1–O3 | ASL, BOLD, T1w/T2w | QC pass counts by group; FreeSurfer eTIV, 4th-ventricle, brainstem and Euler numbers by group (Paper 0). **No CBF, BOLD, ΔR2\* or T1w/T2w-ratio value was ever computed** | All | QC and FS only |
| E2/E3/F | Design and code work (2026-09-25/26) | Simulations on the public design matrix; label-free SCP_R cross-session ICCs and node-survival counts from Paper 0 profiles read **without** the group column; aggregate sidecar counts (f0 values, MagneticFieldStrength), with no IDs; phenotype availability counts by cell (Q14) | — | Availability by cell only |
| Descriptor | Dataset documents | The *Sci Data* descriptor and the preprint (recruitment, group definitions, cohort characteristics) | — | Published group-level description |

**The first author's own viewing: "UNSURE" (binding declaration).**
- The full dataset, all MEGRE included, has been on the lab server under the first author's project since 2026-08-25.
- **He is unsure whether he has viewed raw MEGRE images of any participant before registration, including the two figures (M11).**
- He has never seen a χ map or any group comparison of MEGRE data.
- The most any such viewing could contain is the qualitative magnitude and phase appearance of individual sessions. Neither route could yield a group-level estimate of SN–VTA or dentate χ, or of their association with PRIME or DWI.

**Consequences and handling (registered).**
1. **Wording.** The paper never calls the MEGRE "unseen". It states: *"No susceptibility map, local-field map or ROI value (χ or R2\*) had been computed before registration; whole-head voxelwise R2\* was fitted in memory for feasibility in 4 sessions and summarised as percentiles; the first author is unsure whether he had previously viewed raw MEGRE images."*
2. **Blinding claims.** "Blind" is always qualified as **"procedurally blinded"**. The labels remain public, and before the header-translation mitigation every coded session could have been re-identified from header geometry alone (Annex A.4); the blind rests on commitment plus an audit trail, not on cryptography. Analyst blinding remains valid prospectively, because χ has never existed and all QC, gates and ROI placement run on coded copies.
3. **Residual risk.** Pipeline choices could have been shaped by image appearance (φ₀ handling, the monopolar model, coverage handling). These choices are outcome-neutral, fixed here by consensus or atlas rules, and never tuned on χ. An SN–VTA difference of the tested size (about 0.5 SD, about 10 ppb) is not visually discernible on raw magnitude or phase. Hypothesis directions come from the adult literature.
4. **Exposure-sensitivity analysis.** H1, S3 and T1 are re-run without the four quantitatively inspected sessions **sub-20139 ses-1, sub-21432 ses-1, sub-20976 ses-2 and sub-20188 ses-3**. The set is defined by this fixed list; the rows are dropped, never substituted by a later session. Two of the four are follow-up sessions, which matters for T1. It is applied by real ID after unlock; nobody looks up these sessions' codes before unlock. In T1 only the listed rows are dropped: if a listed session is a subject's occasion 0, that subject's occasion-1 row stays and nothing is re-indexed. A changed conclusion is flagged, never used to override the primary.
5. **Commitments** (written, and repeated in the paper): no uncoded MEGRE, and no image derived from it, is opened after blind setup; the uncoded local copy of sub-20139 ses-1 and `megre_coverage.tsv`, both held in a session temporary folder on the Mac, were removed by the operating system's temporary-folder cleanup between 2026-09-27 and 2026-09-29 (folder observed empty and file absent on 2026-09-29; no copy elsewhere was found by a file-name search) [FILL: the first author confirms, or records his own deletion date]; `megre_coverage.tsv` was never joined to group; the two figures were removed from the shared report folder and deleted or sealed on [FILL: date]; SHA-256 values of these files are in Annex E; χ is computed only after the OSF timestamp.
6. **Seen versus never seen, for S3 (registered statement).** *"All Penn LEAD DWI group-level results, including those for the right SCP, were seen before this registration. No susceptibility value existed at registration, and the association between dentate χ and any DWI metric has never been computed. The pSCP_R segment values are computed label-free by a hashed script, frozen and hashed before unlock, and are not analysed by group. The DWI bundle-QC rater is not blind in practice."*
7. **Co-authors:** [FILL: per co-author prior knowledge, or "none"].

---

## Part 5. Analyses
All models run in **R 4.4.1 on the first author's Mac** from the tables exported at unlock, by the hashed `run_all.R`. **Environment pin (no renv):** sandwich 3.1-1, lmtest 0.9-40, lme4 1.1-38, nlme 3.1-164, MASS 7.3-61, boot 1.3-30, jsonlite 1.8.9, digest 0.6.37, Matrix 1.7-0; on real data `run_all.R` stops if any version differs. pbkrtest, clubSandwich, lmerTest, MBESS, mice, brms, irr, specr and geepack are not installed and are never called; the registered fallback of Q23 is therefore the only branch. The Python verification uses Python 3.11 with statsmodels 0.14.5.

### Q19. Statistical models

**Shared coding.** PS_e and ADHD_e as in Q11; sex_e = ±0.5; race3 with `contr.sum`. Continuous covariates are z-scored within the analysis sample. Cross-sectional models and T1 use `ns(age, 3)` of the session age (T1 adds occ; this is the model simulated for the registered power). The pooled all-session, E-dev and S2 all-session models split age into age_bar (subject mean) and age_w (age − age_bar).

**Primary M1 (H1)** on `first_usable` (n ≈ 100–112):
```r
fit <- lm(chi_SNVTA ~ PS_e * ADHD_e + splines::ns(age, df = 3) + sex_e + race3 +
            scale(roi_vol) + scale(Q_phase), data = first_usable)          # chi_SNVTA = chi_primary
V   <- sandwich::vcovHC(fit, type = "HC3"); dfr <- df.residual(fit)
b <- coef(fit)["PS_e"]; se <- sqrt(V["PS_e", "PS_e"]); p <- 2 * pt(-abs(b / se), dfr)
outA <- (b < 0) && (p < .05)
rob  <- MASS::rlm(formula(fit), data = first_usable, psi = MASS::psi.huber, k = 1.345,
                  scale.est = "proposal 2", maxit = 200)                  # non-convergence -> 1.4826 * MAD(OLS residuals)
sigB <- min(sigma(fit), rob$s)
B_evaluable <- g$PC4 && g$rL_primary >= 0.25                              # g = gates_preunlock.json (Annex A)
dB   <- 0.6 * sqrt(min(1, g$rL_primary))
outB <- B_evaluable && pt((b + dB * sigB) / se, dfr, lower.tail = FALSE) < .05
# simple effects (1, -/+0.5) on (PS_e, PS_e:ADHD_e); ADHD TOST at +/- 0.5 * sigB (K1b)
```
Implementation: `R/m1.R` (`fit_ols_hc3`, `huber_fit`, `h1_decide`, `k1_labels`, `run_h1`).

**Reported with M1** (none changes the decision):
- β in ppb with 95% and 90% CIs; δ = β/σ̂_resid with a cell-stratified bootstrap 95% CI (5,000 resamples); raw Hedges g; a noncentral-t CI for δ;
- the bound in ppb (δ_B·σ_B) and the 90% CI of β next to the adult −14.24 ppb (descriptive; reference region, resolution, field strength and pipeline differ);
- a **Freedman–Lane permutation p** (10,000 permutations; the reduced model omits PS_e and keeps PS_e:ADHD_e, so the null is the average PS effect; its residuals are permuted within ADHD strata; the statistic is the HC3 t; Winkler et al., 2014; exchangeability is approximate, so HC3 remains the decision test);
- analytic normal-likelihood Bayes factors on δ under N(−0.35, 0.25) and N(0, 0.5) priors (Dienes, 2014), and the shrunken estimate under the sceptical prior;
- the ses-1-only refit; the B co-reports (§4.1).

**K-refits.** K2: the same formula with R2\*_primary, one-sided lower-tail p. K3: the same formula with Δref on the sessions with a valid CSF reference; under the trigger, the CSF-referenced refit is run with the full decision rule.

**S1:**
```r
fitS1 <- lm(D ~ PS_e * ADHD_e + splines::ns(age, df = 3) + sex_e + race3 +
              scale(roi_vol) + scale(str_vol) + scale(Q_phase), data = first_usable)
# D = (chi_primary - m1)/s1 - (chi_STR - m2)/s2, pooled label-free constants from gates_preunlock.json (required)
# L = beta_PS(D) - beta_ADHD(D), HC3; components one-sided HC3
```
- Also reported: a stacked two-ROI `nlme::gls` model with ROI-specific intercepts and slopes, `varIdent(~1|roi)` and `corSymm(~rk|sub)` (REML), and per-ROI β_PS and β_ADHD; robust-scaled L if SD/MAD > 1.5.

**S2** on `s2_sample` (n = 83):
```r
fitS2 <- lm(chi_SNVTA ~ scale(PRIME_total) + scale(BDI_total) + scale(ARI_total) + ADHD_e +
              splines::ns(age, df = 3) + sex_e + race3 + scale(roi_vol) + scale(Q_phase), data = s2_sample)
```
- The test is β on scale(PRIME_total): HC3, two-sided.
- **Sensitivities:** + PS_e (the gate for "beyond diagnosis"); within PS only; rank-based inverse-normal PRIME; PRIME items 1–12 without the duration items; CPT commission errors (`cnb_cptnl_total_fp`) in place of dx_adhd; after T1, a between-person all-session PRIME model (PRIME split into subject mean and within-subject deviation; REML maximum). The age-normed PRIME-5 sensitivity of the draft is **dropped** (its normative scoring cannot be applied without tables not held by the team).

**S3** on `s3_sample` (n ≈ 80–108):
```r
fitS3 <- lm(scale(twICVF_pSCP_R) ~ scale(chi_DN_R) + scale(chi_DGM) + scale(b0loc_pSCP_R) + scale(abs_pitch) +
              scale(tract_len) + splines::ns(age, df = 3) + sex_e + race3 + PS_e * ADHD_e + scale(FD) + scale(NDC) +
              scale(log(nSL)) + scale(DN_R_vol) + scale(Q_phase), data = s3_sample)
```
- **S3a:** β(scale(chi_DN_R)), HC3, two-sided. A Freedman–Lane permutation p (reduced model drops chi_DN_R; 10,000 permutations within ADHD strata) is co-reported. HC3 is conservative here (type I 0.009 at .0167 and 0.035 at .05; Q21); it remains the decision test and the power figures of Q21 are for HC3.
- **S3b:** NC1 = β(chi_DN_R) − β(chi_PUT_R), with chi_PUT_R substituted in the same model; NC2 = β(chi_DN_R → pSCP_R) − β(chi_DN_R → MCP_R matched); each a one-sided 95% cell-stratified case-bootstrap interval (5,000 resamples; recomputed with 10,000 if its bootstrap p falls in [.03, .07]).
- **Sensitivities:** + node-wise R2\* along pSCP_R (if β(chi_DN_R) falls by > 50%, the label is "not separable from local relaxation (artefact or mediation)"); neither b0 nor R2\*; the VTAuncross_R 0–30 segment (its nodes 0–16 were the round-1 "resilient" cell, which has been seen); dentate + interposed χ; left dentate χ as the predictor (laterality control); peak alignment; calibrated T1w/T2w in place of tw-ICVF; the 94-subject conservative sample; the exposure-sensitivity set.
- **Node-level profile (Figure 5d only; Tier 4; amendment-2 module):** a Freedman–Lane max-T over the retained pSCP_R nodes. Nodes excluded for a subject are missing for that subject; the max is taken over nodes to which ≥ 80% of the S3 subjects contribute; one vector of uniform keys per permutation is shared across nodes and, within each node and ADHD stratum, the reduced-model residuals of the subjects present are reassigned by the keys' rank order (the same scheme as the specification curve).

**T1** on `two_occ` (occasion 0 = the first usable session, all subjects; occasion 1 = the earliest usable session after it, returners):
```r
two_occ$y_tsd <- two_occ$chi_SNVTA / g$pooled$sd_total_SNVTA   # SD over the frozen first-usable sessions (gate JSON)
fitT1 <- lmer(y_tsd ~ PS_e * ADHD_e + occ + PS_e:occ + splines::ns(age, df = 3) + sex_e + race3 +
                scale(roi_vol) + scale(Q_phase) + (1 | sub), data = two_occ, REML = TRUE)   # at the verified REML maximum
b <- fixef(fitT1)["PS_e:occ"]; se <- sqrt(vcov(fitT1)["PS_e:occ", "PS_e:occ"])
df_c <- containment df of the same term from nlme::lme                                 # registered decision df
notDeclined <- pt((b - 0.4) / se, df_c) < .05 && (b_followup_only < 0)
declined    <- b > 0 && 2 * pt(-abs(b / se), df_c) < .05 && (b_followup_only >= 0)
```
- **Co-reported (never decisive):** a parametric-bootstrap p (2,000 draws from the null-restricted REML fit); an **OLS working-independence fit with `sandwich::vcovCL(type = "HC2")` cluster-robust SE and a restricted wild cluster bootstrap** (9,999 Rademacher draws). This is not CR2 for the LMM estimate and is labelled "OLS cluster-robust (HC2)".
- Δ is also stated in ppb (0.4 × the pooled total SD), next to the occasion-0 estimate (§4.4).
- The follow-up-only M1 (≈ 45–53 subjects, ≈ 8 PS-only) is used for its sign only; an undefined HC3 SE (a leverage of 1) never stops the run.

**Mixed-model rule (Paper 0 lesson iii; applies to T1, the pooled sensitivity, E-dev and S2's all-session model).**
- Each model is fitted with `lme4::lmer` under bobyqa and Nelder–Mead and with `nlme::lme` under nlminb and optim, all scored on one lme4 REML deviance function; **the fit with the highest REML log-likelihood is reported** (ML for likelihood-ratio tests; ties → the listed order).
- A 1-D profile over θ = σ_subject/σ_resid (0–3, step 0.01) confirms the global maximum. If the profile finds a higher likelihood than every optimiser, lmer is refitted from the profile argmin and the best fit is kept; the event is reported.
- A boundary maximum (θ = 0) is stated, and the cluster-robust co-report is then emphasised.

**Pooled all-session H1 sensitivity (after T1; not a test):**
```r
lmer(chi_SNVTA ~ PS_e * ADHD_e + ns(age_bar, 3) + age_w + sex_e + race3 + scale(roi_vol) + scale(Q_phase) + (1 | sub), REML = TRUE)
```

**E-dev** (estimation; all usable sessions; amendment-2 module):
```r
fitDev <- lmer(chi ~ PS_e * ADHD_e + ns(age_bar_c, 3) + age_w + PS_e:age_w + ADHD_e:age_w + PS_e:age_bar_c +
                 sex_e + race3 + scale(roi_vol) + scale(Q_phase) + (1 + age_w || sub), REML = FALSE, data = all_sess)
# parametric-bootstrap LRT of PS_e:age_w (lme4::bootMer, 2,000 draws; pbkrtest is not installed)
```
- If the random-slope variance is singular, it is dropped from both models. A developmental reading is allowed only if the within- and between-person PS × age terms do not differ (Wald p ≥ .10); only the within-person estimate is then interpreted (Neuhaus & Kalbfleisch, 1998; Sørensen et al., 2021).
- A realised-design null calibration (1,000 simulations, saturation b ∈ {0, 0.2, 0.35, 0.5} y⁻¹) runs after unlock and before any E-dev fit.
- Age-confounding sensitivities: common support 10–18 y, overlap weights (Li, Morgan & Zaslavsky, 2018), age matching (caliper 0.5 y). A GAMM (mgcv) is used for display only (Wood, 2011); **never inferential**.

**E-z:** a TD-only model χ ~ ns(age, 3) + sex with leave-one-subject-out predictions; z = (χ − μ̂_TD)/σ̂_TD; centile charts; the proportion with z ≤ −1 (Wilson CIs); the 6 threshold-psychosis youths serve as the E-anchor. These are **not** called a normative model.

**Verification.** `verify_statsmodels.py` re-implements M1 (β_PS, HC3 SE, p, both simple effects), S1 (L and its SE), S2 and S3a with statsmodels OLS/HC3, and T1 (β and SE of PS_e:occ) by a direct numpy/scipy REML profile over θ (the statsmodels MixedLM optimiser failed numerically on one synthetic scenario, so it is not used), on the same table and with its own natural-spline basis (same column space as `ns`). Tolerance: 1e-6 relative for OLS, 1e-4 for the mixed model. On all three synthetic scenarios all 13 checks passed (largest relative difference 2 × 10⁻⁹). A difference beyond tolerance in a decision quantity halts the analysis until the coding error is found (Annex C.4). Declared: both implementations come from the same AI system, so they are not independent analysts, and resampling p values are not re-computed. Every manuscript number is traced to `results.json` through a numbers ledger.

### Q20. Effect sizes and SESOI
- **H1 benchmark:** δ = −0.6, the adult magnitude from established psychosis (Vano et al., 2025, 2026). Neither estimate is from youth or pre-onset samples.
  - The **B bound** δ_B = 0.6·√min(1, r̂_L) (about 0.50 at r̂_L = 0.70) is a heuristic, conservative shrinkage of the observed adult bound, not an attenuation correction (§4.1); it is co-reported in ppb against −14.24 ppb (SE 3.73; Vano 2025, Table 3). Covariate adjustment makes δ larger than a raw g, which makes B conservative for covariate variance.
- **ADHD equivalence:** ±0.5 σ_B.
- **S1:** L = −0.6 (plausible range −0.6 to −0.9).
- **S2:** partial r = 0.25–0.30 (the correlation of χ with the part of PRIME not explained by the other model terms).
- **S3:** latent ρ = 0.5 with **total** right-dentate χ (R_χ = 0.70, R_DWI = 0.55; observed r ≈ 0.31) is the smallest coupling with adequate power, and only if dentate χ shares little variance with χ_DGM (Q21).
- **T1:** Δ = 0.4 total SD.
- **PC2w:** the published youth rise of about 2.6–3.4 ppb (putamen) and 4.4–6.9 ppb (GP) over 1.95 y (Zhang et al., 2018; Li et al., 2014); the SD of the 2-y change is assumed to be 8–16 ppb (Naji et al., 2022, inflated for 1.5 mm and youth).

### Q21. Power
Design-based, on the real first-usable design matrix (n = 112; cells; age–PS correlation 0.27; ns(age, 3); sex; race3; the interaction; simulated standard-normal roi_vol and Q_phase redrawn in every draw). **Source: `sim_typeI_power.R`, which calls the registered R functions (HC3 decision tests, σ_B, the K1 rules, the containment-df T1 test); its output `sim/typeI_power.md` is authoritative.** Rows marked "E2" come from the v1.2 HC3 Python scripts (attachments), which use a `cr` spline basis and one covariate draw; they agree with the R functions within Monte Carlo error where both exist.

| Test | n | Assumption | Power / rate |
|---|---|---|---|
| M1 HC3 type I | 112 | δ = 0; normal / skewed + heteroscedastic errors | 0.038 / 0.043 |
| **H1 (A)** | 112 | δ = −0.3 / −0.4 / −0.502 / −0.6 | **0.26 / 0.44 / 0.63 / 0.80** (MDE80 ≈ 0.61) |
| H1 after QC (E2) | 104 / 100 / 96 | δ = −0.6 | 0.74 / 0.72 / 0.70 |
| **Outcome B** (bound −0.502, r̂_L 0.70) | 112 | true δ = 0 / −0.3; at the bound | **0.740** / 0.22; **false B 0.038** |
| B, true δ −0.1 / −0.2 (E2) | 112 | bound −0.50 | 0.58 / 0.39 |
| B, unadjusted −0.6 (E2) | 112 | true δ = 0 / −0.2 | 0.87 / 0.57 |
| B after QC, δ = 0 (E2) | 104–96 | bound −0.50; −0.6 | 0.71–0.67; 0.84–0.81 |
| B with 2 dominant scans (±6 SD) | 112 | at the bound; δ = 0 | false B **0.052** (MC SE 0.004); P(B) 0.604 |
| K1a "PS effect" / K1b "PS-specific" | 112 | −0.6 in both strata, ADHD 0 | 0.79 / 0.36 (null: 0.017 / 0.008) |
| ADHD TOST ±0.5 σ_B | 112 | ADHD δ = 0 | 0.47 |
| K1c "comorbid-concentrated" | 112 | PS-only 0, comorbid −0.8 | 0.38 |
| PS × ADHD interaction (E2) | 112 | — | SE 0.395; MDI ≈ 1.1, **not testable** |
| K3 trigger, no true drift | 112 | SD(Δref) = 0.2 / 0.5 residual SD | P(trigger) 0.004 / 0.07 at δ = −0.6 (0.02 / 0.08 at δ = 0); P(A) 0.785 → 0.784 / 0.777 → 0.756; P(B) 0.742 → 0.721 at 0.5 SD |
| Freedman–Lane (999 permutations) | 112 | δ = 0 / −0.6 | 0.054 / 0.84 |
| **S1** L type I | 112 | null, ρ = 0.5 | 0.015 at .0167; 0.036 at .05 |
| S1, L = −0.9 (sim) | 112 | ρ(SN–VTA, STR) = 0.5 | 0.56 at .0167 (0.74 at .05); "dissociation" 0.42 |
| S1, L = −0.6 (E2) | 112 | ρ = .3 / .5 / .7 | 0.16 / 0.24 / 0.43 at .0167 (0.29 / 0.40 / 0.61 at .05); "dissociation" 0.10–0.33 |
| S1, L = −0.8 (E2) | 112 | same | 0.32 / 0.45 / 0.71 (0.49 / 0.63 / 0.85); "dissociation" 0.21–0.60 |
| **S2** | 83 | null; partial r = 0.25 / 0.30 | type I 0.013 / 0.020 / 0.039; **0.38 / 0.45 / 0.56** and **0.56 / 0.62 / 0.73** at .0167 / .025 / .05 |
| **S3a** | 108 | null | type I 0.009 / 0.017 / 0.035 at .0167 / .025 / .05 |
| S3a, R²_x = 0 | 108 | latent ρ = 0.4 / 0.5 | **0.38 / 0.46 / 0.58** and **0.63 / 0.70 / 0.80** at .0167 / .025 / .05 |
| S3a, R²_x = 0.3 | 108 | latent ρ = 0.4 / 0.5 | 0.27 / 0.32 / 0.44 and 0.42 / 0.49 / 0.62 |
| S3a, R²_x = 0.5 | 108 | latent ρ = 0.4 / 0.5 | 0.16 / 0.21 / 0.30 and 0.28 / 0.34 / 0.47 |
| S3a with FA_t_pc (R 0.26) (E2, model-based) | 108 | latent ρ = 0.5 | 0.35 / 0.51 (why FW-DTI is not primary; HC3 lower) |
| S3b (E2) | 92–101 | — | 0.18–0.42 (label only) |
| **T1** (containment df) | 53 returners | true β(PS × occ) = Δ (boundary), ICC 0.7 | 0.056 (MC SE 0.007) |
| T1, deficit (0.5 total SD) persists | 53 | ICC 0.5 / 0.7 / 0.8 | 0.46 / **0.61** / 0.70; **given A: 0.36 / 0.56 / 0.69** |
| T1 false "did not decline" if a 0.5-SD deficit vanishes (E2; reviewer check 0.014–0.015) | 53 | — | 0.01–0.03 |
| PC2w (E2) | 48 pairs | Zhang-size change, σ_D 8 / 12 / 16 ppb | PASS 0.96 / 0.74 / 0.53; no change: PASS 0.05 / FAIL 0.50 |
| G2 [reviewer calculation] | 20 | true ICC = 0.90 | P(point ≥ 0.90) ≈ 0.5: the engine choice is close to random at that agreement level |
| K4 "fragile" [reviewer calculation, registered functions, 300 draws] | 112 | δ = −0.6; δ = 0; false-positive A under the null | P(fragile \| A) 0.09; P(fragile \| B) 0.10; P(fragile \| A) 0.36 |

- Converter-only dilution (δ 0.09–0.18) gives P(A) of only 0.06–0.14 (E2).
- **Registered honest reading:**
  - *"≈ 0.8 power only for an adult-sized deficit before QC loss (0.80); about 0.70–0.74 after the expected QC loss."*
  - The study cannot resolve δ ≈ 0.3, converter-restricted effects or interactions. The specificity labels will often be "not established".
  - H1 power assumes that roi_vol and Q_phase are independent of group; motion or phase noise that differs by group would inflate SE(β_PS). The K3 rule costs about 0.02 power for A and for B when SD(Δref) is half a residual SD, and nothing when it is small.
  - S2 has 0.38–0.56 power for partial r = 0.25 and 0.56–0.73 for 0.30 across the Holm steps.
  - S3 has about 0.6–0.8 power only for latent ρ ≈ 0.5 with total dentate χ that shares little variance with χ_DGM, depending on the Holm step.
  - T1 is informative only if the cross-session ICC is ≥ 0.7 (power 0.56 given A at ICC 0.7).
- **Realised values posted in amendment 2, before unlock (label- and age-free):** n per sample, ρ(SN–VTA, STR), r̂_L and R_SB, r_cross, the cross-session ICC(3,1) of χ_primary over the PC2w pairs (not age-adjusted; T1 is labelled "uninformative" if it is < 0.40), SD(Δref) in residual-SD units with the implied P(K3 trigger | no drift) and adjusted power from the table above, and r(χ_DN_R, χ_DGM) with the implied S3 power row. Age-residualised and TD-only ICCs are post-unlock descriptives.

### Q22. Inference criteria
- **H1:** A at two-sided α = .05; B at one-sided α = .05 against δ_B (different nulls; no correction between them). The final outcome follows the K3 rule (§4.1).
- **K-labels:** each component at its own α (intersection–union); one-sided where a direction is required (K2; the TOST halves; the "dissociation" components). K4 "robust" by the joint permutation test (§24.3).
- **Key secondary:** Holm across {S1, S2, S3a} at .0167 / .025 / .05, two-sided, with the direction requirement. S3a leaves the family if PC4-DN, G2-DN or PC5 fails, leaving S1 and S2 at .025 / .05. **D12:** S3 stays in the family (confirmed with the Q21 figures).
- **T1:** one-sided α = .05 against the fixed margin with containment df, conditional on the final outcome A.
- **Exploratory:** BH-FDR q < .05 within each named block (X1–X7).
- **Node-level inference:** only a Freedman–Lane max-T (≥ 5,000 permutations, 10,000 if decision-relevant; identity included; p = (#≥ + 1)/(N + 1)). Pointwise or simultaneous GAMM bands are never inferential (Paper 0 lesson ii).
- **Monte Carlo precision.** Every resampling p used for a decision or label (T1 PB co-report excepted) uses ≥ 10,000 draws or is recomputed with 10,000 if it falls in [.03, .07]; its Monte Carlo 95% interval is always reported.
- **Bayes factors** are reported, never decisive (Dienes, 2014). Every significant estimate is also shown as a sceptical-prior shrunken value (Button et al., 2013).

### Q23. Assumption violations and non-convergence
- **Heteroscedasticity and non-normality.** HC3 is the decision SE for all OLS tests (calibrated type I 0.038–0.043). Reported alongside: Freedman–Lane permutation p, a bootstrap CI, a rank-based ANCOVA, a noncentral-t CI for δ, and a Huber M-estimate.
- **Mixed models:** the REML-maximum rule (Q19); a singular random slope in E-dev is dropped from both models; **df by nlme containment** (the registered decision); a parametric bootstrap (lme4, 2,000 draws) and the OLS cluster-robust (HC2) fit with a wild cluster bootstrap (9,999 Rademacher draws) are co-reported. Kenward–Roger and CR2 are not used (pbkrtest and clubSandwich are not installed); this is fixed now, not at amendment 2.
- **Other fixed branches:** Bayes factors by the analytic normal likelihood; noncentral-t CI by root finding; weighted κ hand-coded; the specification curve hand-coded (§24.3); GEE by statsmodels; puberty imputation as in Q14.
- **Gate JSON:** missing or invalid decision fields stop the run (hard stop); gate *failures* follow the registered branches (Q24) and do not stop it.
- **QSM non-convergence or failure** in a session or variant makes that value missing (§12.4); it never triggers re-tuning for that session.

### Q24. Reliability and robustness

**24.1 Gates** (outcome-neutral; computed by hashed scripts on coded data, label- and age-free, through `audit_run.py`; posted in amendments 1 and 2).

| Gate | Registered criterion | If it fails |
|---|---|---|
| **Sign check** (pilot; `pilot_checks.py sign`) | Pooled median, over the 10 pilot sessions, of the whole-brain-referenced ROI means of GP, SN–VTA, RN and right dentate: all four > 0 | All four < 0: the phase sign is negated for every session (`run_qsmxt.sh --stage negate`), all runs are recomputed and the pilot χ is discarded (amendment 1). Mixed signs: no flip; one diagnostic revision from a closed menu (phase-scaling check from the header audit; Mask 3 coverage of the four ROIs; reference support share); if the signs are still mixed, the methods-report path applies |
| **PC1** ordering | In ≥ 95% of usable sessions with all six ROIs plausible, χ(GP) > 0 and χ(GP) > χ(thalamus), and the median χ is highest in GP among {GP, SN–VTA, RN, putamen, caudate, thalamus} (Treit et al., 2021) | The fallback pipeline **QSMxT PDF + iLSQR** (for either engine) with the gates re-evaluated (`fallback_PC_pass`). Rationale: the sign check and the mask assertion precede PC1, so a residual ordering failure is most likely a background-removal/inversion problem; iLSQR's known underestimation of deep grey matter is declared. If the fallback also fails, the paper becomes a methods report and H1 is reported only descriptively ("pipeline not validated") |
| **PC2w** within-person development | In returners' first two usable sessions, Δ̄ = mean of (Δ_PUT + Δ_GP)/2. **PASS:** Δ̄ > 0 (one-sided paired t p < .05), mean Δ_PUT ≥ 0 and mean Δ_GP ≥ 0. **FAIL:** Δ̄ ≤ 0. **INDETERMINATE:** otherwise. Δ(thalamus) and Δ(ROI) − Δ(thalamus) are reported; if the software-version class differs within > 20% of pairs, same-class pairs are reported separately | FAIL → the PC1 fallback path. INDETERMINATE → no switch; after unlock PC2x (age ↔ putamen/GP χ, one-sided p < .05 each) is computed first; if PC2x also fails, every H1 statement carries "developmental positive control not confirmed" |
| **PC3** convergence | r(χ_GP, R2\*_GP) ≥ 0.30 across first-usable sessions | As PC1 |
| **PC4** precision | **Split-echo (noise) reliability, an upper bound:** R_SB = 2r/(1 + r) ≥ 0.70, r = Pearson r between the {1,3} and {2,4} reconstructions of χ_primary (bilateral ROI) over the frozen first-usable sessions of the selected branch whose split values pass their own plausibility rules. **r̂_L** = Spearman–Brown of the lower limit of the 95% Fisher-z CI of r (z ± 1.96/√(n − 3)); if undefined, −1 | B is not evaluable; H1 labelled "low-precision measure"; attenuation-corrected δ reported |
| **PC4-DN** | The same for χ_DN_R, ≥ 0.70 | S3 exploratory |
| **PC5** DWI (real ID, label-free; after posting; `pc5_dwi.py`) | (a) pSCP_R tw-ICVF and whole-tract MCP_R tw-ICVF each rise with age (one-sided p < .05) in `lm(y ~ age + sex + FD + NDC)` on each subject's first scalar-pass DWI session; (b) \|r(ISOVF, tw-ICVF)\| < \|r(ISOVF, ICVF_NDI)\| over the retained pSCP_R node observations | S3 exploratory |
| **G2** engines (`gates_preunlock.py g2`) | **ICC(3,1)** (two-way mixed, consistency; McGraw & Wong, 1996; Shrout & Fleiss, 1979; F-based 95% CI) of SN–VTA χ between QSMxT and SEPIA on the 20 sessions of `gate20.txt`: **point ≥ 0.90 and lower 95% bound ≥ 0.75**; ICC(2,1), Pearson r and the Deming slope co-reported. λ per engine (Annex B.4). Rationale: Cornell MEDI is the reference implementation, and because the SEPIA arm shares unwrapping, echo combination and Mask 3, any disagreement is confined to background removal and inversion | SEPIA becomes the primary engine (decided blind). **If A2 is declined:** G2 is within-QSMxT algorithm agreement (V-SHARP+MEDI vs V-SHARP+RTS vs PDF+iLSQR, ICC(3,1) point ≥ 0.75); the engine stays QSMxT whatever the outcome, and a failure adds the caveat "inversion-algorithm agreement not established" |
| **G2-DN** | The same criteria for χ_DN_R | S3 exploratory |
| **G3** method agreement | ICC(3,1) of SN–VTA χ across the 9 background-removal × inversion nodes on their common support, over the first-usable sessions: point ≥ 0.75 (ICC(2,1) co-reported) | The primary ROI becomes SN-only (CIT168 SNc+SNr), and every gate constant is recomputed for it. If SN-only also has ICC(3,1) < 0.75, H1 carries the caveat "method-dependent ROI measure" |
| **G4** QC rater | Linearly weighted κ ≥ 0.60 between pass 1 and a random 20% re-rated ≥ 14 days later (**intra-rater**: one rater) | One blind rubric refinement, a full re-rating (pass 3) and a new 20% re-rating (pass 4); if κ still < 0.60, the pass-3 grades are final and every H1 statement carries "QC reliability not established". No further iteration |
| **G5** registration | Per session: T1w→MNI brainstem + ventral-DC Dice ≥ 0.80, **and** MEGRE↔T1w SynthSeg Dice ≥ 0.80 for brainstem (16) and cerebellar WM (7/46) and ≥ the VDC threshold for ventral DC (28/60). The VDC threshold is 0.80 unless > 20% of the pilot sessions fall below it on VDC alone, in which case it is the pilot's 10th percentile of VDC Dice (never below 0.70), fixed in amendment 1 | The session fails first-usable criterion (e). If > 10% of sessions fail, the chain is revised **once**, blind, from a closed menu (N4 on the T1w as well; the rigid metric MI → CC; the brain-mask dilation 2 → 4 voxels), applied to all sessions; if > 10% still fail, G5 is recorded failed and the methods-report path applies. At unlock, Dice is compared by group (Wilcoxon); if p < .05 it is added as a covariate in a co-reported sensitivity |

**Branch order.** 1. G2 selects the engine. 2. PC1, PC2w and PC3 are evaluated on that engine's primary pipeline (fallback: QSMxT PDF + iLSQR). 3. G3 selects the ROI. 4. PC4-DN, G2-DN and PC5 decide S3's status. 5. PC4, r̂_L and the pooled constants are computed for the selected branch and ROI.

**Stage II changes** (after χ exists, before unlock; Annex C.3) are allowed only when a registered rule requires them or when a pre-named label-free metric (the plausibility session indices, the mask assertion, G5 Dice, the post-run file checks) reveals a technical failure; a change applies to all sessions, cites only those metrics, and is posted in amendment 2 with before/after aggregates.

**24.2 Reliability reported in addition.**
- **Cross-session reliability (pre-unlock, label- and age-free):** Pearson r and ICC(3,1) of χ_primary over the returners' first two usable sessions (not age-adjusted; heterogeneous developmental change makes it a lower bound). ≥ 0.60 counts as adequate; < 0.40 makes T1 "uninformative". The realised ICC fixes which row of the T1 power table applies. **After unlock:** ICC(3,1) of age-residualised χ in all returners and in TD alone (Koo & Li, 2016), descriptive.
- **DWI reliability, already measured label-free** (pSCP_R nodes 0–30, ISOVF ≤ 0.30, ses-1 vs ses-2 over about 1.7 y, a lower bound): tw-ICVF ICC(3,1) 0.57 [0.35, 0.73] (52 pairs) and 0.54 [0.33, 0.70] (57 pairs); ICVF_mod 0.45/0.42; FA_t_pc 0.27/0.26; MD_t_pc 0.31/0.30. Split-half streamline reliability (20 random halvings) is computed by `dwi_geometry.py`.

**24.3 QSM specification curve for H1 and the K4 "robust" label** (Simonsohn et al., 2020; Salman et al., 2025; `R/spec_curve.R`).

| Dimension | Levels (* primary) |
|---|---|
| Background removal | V-SHARP*, PDF, LBV |
| Inversion | MEDI*, iLSQR, RTS |
| Reference | whole brain*, lateral-ventricle CSF (missing where no valid CSF reference) |
| ROI summary | primary-ROI weighted mean w ≥ 0.5*, primary-ROI median w ≥ 0.5, the other CIT168 complex's weighted mean (SN-only if the primary ROI is SN–VTA; SN–VTA under the G3 fallback) |
| Covariates | full M1*, minimal (PS_e × ADHD_e + ns(age, 3) + sex_e), full + \|pitch\| |
| Sample | all first-usable*, without threshold psychosis, visual QC grade 0 |

- 3 × 3 × 2 × 3 × 3 × 3 = **486 specifications**, all on the primary branch's frozen first-usable sessions and on the 9-node common support; a variant value that failed its own plausibility rules is missing in that specification only.
- **Evaluable:** n ≥ 30 and ≥ 80% of the first-usable subjects of its sample level. "Robust" is assessable only if ≥ 400 specifications are evaluable.
- **Joint test:** statistic = the median over evaluable specifications of the specification-wise HC3 t of β_PS. For each of **10,000 permutations** one vector of i.i.d. uniform keys is drawn over the first-usable subjects and shared by all specifications; within each specification and ADHD stratum, the reduced-model (no PS_e) residuals are reassigned by the rank order of the keys (Freedman–Lane) and the full model is refitted. One-sided p = (#{median t\* ≤ median t_obs} + 1)/(N + 1).
- **Co-reported (descriptive):** median δ; the share of evaluable specifications significant in the predicted direction; the Stouffer Z (Σ t-to-z/√k), each with its own joint-permutation p.
- **"Robust"** = final outcome A (A or A+B) and ≥ 80% of evaluable specifications with δ < 0 and joint p < .05. Its false-positive rate is bounded by P(A) under the null (0.017).

**24.4 H1 sensitivities** (each reported separately; `R/sensitivity.R`; a row whose input is absent is reported as "not run").
- **Sample:** the Q13 exclusions; the exposure-sensitivity set (Q18); without the (c′) plausibility exclusions; automated exclusions only; strict QC (grade 0); without code D sessions; implausible voxels retained.
- **Model and covariates:** Huber M-estimator; + puberty (complete case and the single imputation) and + menarche; + intracranial volume; + mood diagnosis; epoch as a covariate and early-epoch sessions excluded; + f0 cluster; + years since diagnosis; + |pitch|; + DWI FD; + registration Dice; + the motion index; + the reference-domain WM and CSF fractions and the brain-in-FOV fraction; the coverage-harmonised reference; without race; + Hispanic ethnicity; rank-based ANCOVA.
- **Always co-reported:** the ses-1-only analysis; β_PS for the SN-only and the VTA + PBP sub-regions (descriptive). Without flow compensation (to be confirmed in the header audit), the medial SN–VTA next to the basal vein, the posterior cerebral artery and the interpeduncular cistern is exposed to flow errors and venous partial volume (0.3–0.6 ppm, below the 1-ppm rule); the sub-region contrast is discussed.
- **Engine:** H1 refitted with the other engine on all sessions (SEPIA if QSMxT is primary, and the reverse), if A2 is approved.
- **Conditional:** a medication analysis if records arrive (A10); labelled post hoc if run after unlock.

**24.5 DWI multiverse (S3; Tier 4; amendment-2 module).** Tract: pSCP_R\*, VTAuncross_R 0–30, and (optionally, A3-4) the DSI Studio AutoTrack DRTT_rl from ds007091. Weighting: tissue-weighted\*, modulated, plain. Node rule: per-subject\*, fixed 0–30 with ISOVF only, peak-aligned. n_min: registered value\* or 6. Metric: tw-ICVF\*, ICVF_mod, FA_t_pc, MD_t_pc, FW_pc, MK, T1w/T2w. Relaxation covariates: b0\*, b0 + R2\*, neither. QC list: coded session-level\*, the Paper 0 subject-level list. Covariates: ± NDC, ± ICV, ± 4th ventricle. An AutoTrack result is "convergent" if the subject-level r ≥ 0.50 and the group-effect estimates share sign with overlapping 95% CIs; otherwise "pipeline-dependent".

### Q25. Exploratory analyses (never claimed; BH-FDR within block; amendment-2 modules)
- **Tier 3, estimation:** E-dev (within- vs between-person PS × age for SN–VTA, STR and dentate χ); E-ADHD (the ADHD main effect on χ_STR); E-int (the PS × ADHD interaction and the four cell means); E-z; E-anchor.
- **Tier 4, exploratory** (capped at 7):
  - **X1:** dentate and RN χ group effects;
  - **X2:** R2\*–QSM closed-form χ-separation (Dimov et al., 2022; no new install), with the relaxometric constant and the assumed spatially constant R2 taken from that paper and written into the amendment-2 module before any χ exists [UNVERIFIED: the numeric values]; with data from one field strength and no R2 map, the separation is model-driven (declared);
  - **X3:** sex × PS;
  - **X4:** right dentate χ and pSCP_R tw-ICVF vs in-scanner 2-back d′;
  - **X5:** cerebellum–VTA resting FC (XCP-D, 4S456 parcels) vs SN–VTA χ, only if A3-7 is approved;
  - **X6:** RSAS social anhedonia (cross-sectional) vs SN–VTA χ;
  - **X7:** midbrain CBF (ASLPrep, 3.75 mm, whole-midbrain signal) vs SN–VTA χ (D7 = yes; needs A3-6).
- Anything else is labelled post hoc.

---

## Part 6. Statement of integrity
- The authors confirm that this preregistration was written in good faith, and that it reports every piece of prior knowledge of the data we are aware of (Q18), including the first author's "unsure" viewing status.
- No other preregistration exists for these hypotheses and these data. Paper 0 is a separate, completed DWI study with its own plan, and it tests no LEAD-χ hypothesis.
- We will report every registered analysis, including null and disconfirming results, and every deviation (Annex C).

---

## Annex A. Procedural blinding (registered)

**A.1 Principle and disclosure.**
- There is no escrow holder (user decision 1). Blinding is procedural: group labels, ages and phenotypes are locked in a separate encrypted file inside the project, coded MEGRE copies have scrubbed headers and per-session translation offsets, QC and gates are group-blind, code SHAs are registered, and labels are unlocked once, at the confirmatory stage.
- **Registered disclosure:** *"Blinding was procedural. The group labels remained publicly available on OpenNeuro and in the first author's own earlier analysis files; the lock covered only the project's working copies. The first author held both the analysis and the lock key and could technically have unblinded at any time; before the header-translation mitigation, every coded session could also have been re-identified from header geometry. The blind therefore rests on commitment and an audit trail, not on cryptography. The lock key never entered an AI session. All pre-unlock outputs are aggregates computed by hashed scripts, checked by an audit wrapper that rejects per-code values, and time-stamped in OSF amendments before the single unlock event, which is itself logged."*
- Every pre-unlock gate is **free of age and labels**. In v1.1 an escrow holder was needed because per-code age plus sex re-identifies 160 of 175 sessions; v1.2 removed the need (PC2 replaced by PC2w, the E-dev calibration moved to after unlock, the pre-unlock ICC not age-adjusted).

**A.2 Roles.**

| Role | Who | Holds | Commits not to do before unlock |
|---|---|---|---|
| Lock custodian | First author | `lock_priv.pem` (RSA-3072, AES-256-encrypted) on the Mac; the passphrase kept offline; an offline backup | Decrypt anything; type the passphrase anywhere except OpenSSL's own prompt during `unlock_and_merge.py`; copy the private key to the server |
| Analyst and single QC rater | First author, assisted by an AI coding agent | The coded tree, hashed scripts, aggregate JSON, coded QC images, real-ID DWI tables | See A.4 |
| AI agent | Claude (Anthropic) | No key and no passphrase | Receive the passphrase; run any decrypt command; read `lock_priv.pem` or `*.enc`; open code-keyed continuous tables. A `CLAUDE.md` in `LEADchi/` states these rules for every session |

**A.3 Lock and coded copy.**
1. **Key pair.** Before posting, the user generates the key pair (`openssl genpkey -algorithm RSA -pkeyopt rsa_keygen_bits:3072 -aes-256-cbc`) and derives `lock_pub.pem`. Only the public key goes to the server. Its SHA-256 (the lock fingerprint) is registered.
2. **Label lock (`lock_labels.py`, run once by the user).** `participants.tsv`/`.json`, every `sub-*_sessions.tsv` and the phenotype TSVs used (PRIME, BDI-1A, ARI, Tanner boy/girl, substance, short CPT, RSAS) are packed in memory into one tar, encrypted to `lock_pub.pem` (random session secret wrapped with RSA-OAEP; AES-256-CBC with PBKDF2, 200,000 iterations) and written as `locked/labels_locked.tar.enc` + `labels_locked.sessionkey.enc`, with a manifest of per-file plaintext SHA-256 values. The three hashes are registered. The public originals remain on OpenNeuro; the pipeline scripts refuse to open label, age or phenotype files (`forbid_label_inputs`).
3. **Key table (`blind_setup.py`).** [code, ses_order, participant_id, session_id, epoch, f0 cluster, translation offset] is built **in memory only**, encrypted the same way and written as `blind/key.csv.enc` + `blind/key.sessionkey.enc` (hashes registered). The session secret travels by file descriptor, never on a command line.
4. **Codes.** `secrets.token_hex(3)` subject codes. Coded ses-k labels keep within-subject chronological order; processing order is randomised.
5. **Header scrubbing and geometry.** NIfTI `descrip` (which holds the acquisition time), `aux_file`, `intent_name`, `db_name` and extensions are blanked. JSON sidecars hold only EchoTime, RepetitionTime, FlipAngle and `MagneticFieldStrength` = round(f0/42.576, 4). **Every image of a session (MEGRE, T1w, T2w, QSIPrep ACPC T1w and SynthSeg) receives the same random world translation** (integer mm in [−40, 40] per axis, in both qform and sform): voxel data, rotations, relative geometry and the B0 direction are unchanged, but a coded header no longer matches its raw header. No QSIPrep transform is copied (route R). Every coded file is scanned for the real ID.
6. **Phase.** Stored as float32 radians = raw·π/4096; `blind_setup.py` asserts per echo that the range matches the registered mode and spans > 1.9π.
7. **Draws.** `blind/pilot10.txt`, `blind/lcurve2.txt` (the first 2 pilot sessions) and `blind/gate20.txt` (codes only; seed 20260924; among sessions with a usable T1w; the G2 set disjoint from the pilot, 20 distinct subjects). Pilot and gate sessions stay in the analysis and are declared.
8. **DWI-QC code set (`dwiqc.py`).** An independent code set for coded DWI bundle QC (renamed pyAFQ figures). Its code↔ID file is plaintext in `LEADchi/dwiqc/`, hashed, and not opened until the DWI-QC freeze. It links nothing to χ.

**A.4 What the analyst can and cannot see before unlock.**
- **Can see:** the coded subject and session structure; coded QC montages (gross failures) and T1w/echo-1 magnitude overlays (placement); per-code **pass/fail flags only** (QC grade, plausibility, mask assertion, G5, S3 node count); aggregate JSON and software logs; real-ID DWI tables (non-blind by construction).
- **Cannot see** (by commitment; code-keyed tables are read only by hashed scripts run through `audit_run.py`, which rejects any aggregate output containing a coded ID or a numeric vector longer than 12): any per-code continuous value (χ/R2\* ROI means, split-echo, engine or multiverse values, ROI/SynthSeg volumes, ICV, Dice values, orientation angles, Q_phase, the motion index, coverage and reference fractions); raw ds007116 MEGRE after setup; the round-1 real-ID MEGRE files (Q18); any comparison of coded with raw files (by hash, header geometry or voxels); group analyses of the real-ID pSCP_R tables.
- **No-link rule.** No quantity with a per-code counterpart is computed by real ID before unlock. T1w-derived scalars are computed once, on the coded side. The DWI geometry job outputs node distances, flags and DWI values only, never a dentate volume or χ. DWI is never recomputed on the coded copy. Code-keyed tables stay on the server until amendment 2 is posted.
- **Declared residual channels:** coded T1w and T2w voxels are identical to the raw ones (a voxel comparison, forbidden by commitment, would re-identify them); the number of sessions per subject is public; the f0 cluster is recoverable from the coded B0 value; the software-version class is visible per code; Paper 0 holds real-ID eTIV and 4th-ventricle volumes. **Uniqueness count** (header audit, aggregate): of 181 sessions, [FILL: n] are unique on (sessions per subject × f0 cluster × software class × epoch).

**A.5 Group-blind QC.**
- **Single-rater design (declared).** One rater, the first author, grades every session; he is also the person whose prior viewing is "unsure". G4 measures **intra-rater** consistency only.
- Coded sessions are graded under random QC IDs (the rater never sees codes or session order) as 0 (pass), 1 (minor) or 2 (fail) for gross failures only: motion ghosting, residual wraps, streaking across the whole midbrain or dentate level, a mask excluding the midbrain or dentate, R–L phase wrap-around, air artefact in ROIs, coverage, and **code D (dental or metal artefact), recorded for every session in which it is visible**.
- Grading uses montages of magnitude echoes 1 and 4, the total and local fields and χ (greyscale window ±0.2 ppm, no colour-bar values) at the dentate, midbrain and basal-ganglia levels. **ROI placement is rated on T1w and echo-1 magnitude, never on χ.**
- Automated plausibility flags force a review; automatic grade 2 cannot be overridden. A random 20% is re-rated ≥ 14 days later (G4; stop rule in Q24). The rubric is attachment `qc_rubric.md`.
- A registered sensitivity uses automated exclusions only (Q13), so that the visual grades are never the only route to a primary-sample decision.

**A.6 Sequence.**
1. **Posting** (end of week 1). All code that computes a gate, a branch, δ_B, a first-usable flag or any Tier 0–2b quantity, the K4 "robust" label included, is hashed at posting (Annex J).
2. Pilot, sign check, λ, G5 calibration and G2 → **amendment 1**, posted before the full run (and, if A2 is approved, the hashed SEPIA runner, before any SEPIA χ exists).
3. Full run → plausibility → ROI and multiverse tables → QC → DWI table → `build_usable.py --stage usable` → `gates_preunlock.py --stage final` (aggregate-only JSON with a fixed schema; `audit_run.py` rejects any output holding a per-code value). Contents: every gate and branch, r̂_L and R_SB (primary ROI and dentate), the pooled label-free mean, SD and total SD of χ_primary and χ_STR, the cross-session r and ICC, ρ(SN–VTA, STR), r(χ_DN_R, χ_DGM), SD(Δref), CSF-validity counts, the G5 threshold and fail share, n_min and PC5.
4. **Freeze manifest** (`build_usable.py --stage manifest`): the SHA-256 of the QC tables and κ, the usable lists, the plausibility flags, the gate JSON, the coded ROI and multiverse tables, the real-ID DWI S3 table and the DWI-QC freeze, the analysis code (`run_all.R`, `R/*.R`, schemas, `verify_statsmodels.py`), and the audit-log heads.
5. **Amendment 2** posts the manifest, the gate JSON and the realised values of Q21, with a diff of the analysis code against the version hashed at posting. **Permitted changes:** (i) bug fixes and technical fixes justified without per-code values; (ii) the amendment-2 modules of Annex J (Tier 3–4 only), which must implement the registered text verbatim; their code is posted with the diff. Nothing that affects a Tier 0–2b decision may change.
6. **Unlock.** `unlock_and_merge.py` runs once, on the Mac, after amendment 2: it verifies the ciphertext, manifest and freeze hashes and the order of events, asks for the confirmation phrase, lets OpenSSL prompt for the passphrase, writes the merged and multiverse tables to `unlocked/`, and logs its own hash and time.
7. **Fixed order afterwards:** QC-by-group, orientation, motion, coverage and plausibility by cell → PC2x and the E-dev realised-design calibration → H1 + K1–K4 → S1–S3 (χ and DWI merge here; node-wise R2\* along pSCP_R is computed now by a frozen script) → T1 → robustness (sensitivity battery, specification curve, exposure sets) → estimation and exploratory analyses → the statsmodels verification.

**A.7 Audit log.**
- `audit_run.py [--aggregate …] -- <cmd>` appends one JSON line per run to `LEADchi/logs/audit.jsonl`: UTC times, host, user, the command, the SHA-256 of the script, of every existing input and of every declared output, the git commit if any, the exit status and any rejection; each line carries the SHA-256 of the previous one. The pipeline scripts also keep a hash-chained `logs/pipeline_runs.jsonl`.
- The logs are pushed to a private GitHub repository (A8); their head hashes are quoted in each amendment. AI session transcripts are kept. A read-only check after every step confirms that nothing was written outside the project tree.

---

## Annex B. Processing pipeline (pinned; registered)

**B.1 Environment.** Everything runs on the lab server inside `/analysis/halilaziz/halil/Penn_LEAD/LEADchi/` under `nice -n 10`, with ≤ 16 concurrent jobs. `env.sh`, sourced by every job, redirects HOME, the XDG paths, MPLCONFIGDIR, MATLAB_PREFDIR, TMPDIR, the Apptainer cache and TemplateFlow into the tree, sets `QSMXT_NO_UPDATE_CHECK=1` and points the QSMxT caches into `LEADchi/tools/`. Containers run as `apptainer exec --containall --no-home -B $LX`. `software_versions.txt` records every version, tag, digest and binary SHA-256.

**B.2 QSM engine and primary chain.**
- **QSMxT v9.22.0** (commit `9bffc6e726b09670a68790c71643cca811c5dfa9`; QSM.rs v0.38.0; GPL-3.0; Stewart et al., 2022), installed from the musl tarball `qsmxt-v9.22.0-x86_64-unknown-linux-musl.tar.gz` (sha256 `6593438c029749b5e3fb5c78ce485f6f5e282545b1892832cf30539b0326ca49`) unpacked by hand into `LEADchi/tools/qsmxt/`; `install.sh`, `qsmxt update` and `qsmxt example` are never run.
- **Archived before posting:** `qsmxt --version`, `qsmxt run --help`, `qsmxt init` (the defaults of the pinned binary); per run, `pipeline_config.toml`, `methods.md`, the provenance JSON and the "B0 direction" log lines.
- Every flag below was checked against the v9.22.0 `cli.rs` (web source, 2026-09-26); `run_qsmxt.sh` aborts if any is absent from the captured `--help` (then Annex C.2 applies).

| Step | Registered setting |
|---|---|
| Input | Coded BIDS; phase in radians (A.3); `MagneticFieldStrength` = round(f0/42.576, 4); B0 direction from the affine; **`--obliquity-threshold -1`** (no resampling of oblique acquisitions; χ is never resampled; the slab obliquity is reported by the header audit) |
| Masks | **Mask 1** = SynthStrip (Hoopes et al., 2022; qsiprep-26.0.0 container; its default model and default CSF-inclusive border, i.e. without `--no-csf`; the model file hash is recorded) on the root-sum-of-squares magnitude. **Mask 2** = QSMxT phase-quality map (`qsmxt quality-map`) above its **Otsu threshold within Mask 1** (Otsu, 1979). **Mask 3** = fill-holes(Mask 1 ∩ Mask 2) **minus the first and last 2 slab partitions** (the voxel axis of length 96), supplied with `--use-custom-masks mask3`. A pre-run check requires a mask for every session; a **post-run assertion** (§12.4) invalidates any session whose reconstruction did not use it. **Single pass: `--two-pass` is not used** (it would replace Mask-3 values by a second reconstruction on a holey phase-quality mask, preferentially in iron-rich nuclei). Technical fallback if the custom-mask glob cannot match session directories: `--mask-preset bet-and-phase`, stated in amendment 1 |
| Phase offset | `--phase-offset-removal true --phase-offset-sigma 4 4 4` (voxels; MCPC-3D-S type; Eckstein et al., 2018) |
| Unwrapping | ROMEO (Dymerska et al., 2021), **template mode**: `--no-romeo-individual --romeo-template 2` (1-indexed) |
| Echo combination | `--b0-estimation weighted-avg --b0-weight-type phase-snr`, i.e. w_n = \|S_n\|·TE_n (not the inverse-variance weight \|S_n\|²·TE_n², which the CLI does not offer as such; declared) |
| Background removal | V-SHARP, `--vsharp-max-radius 12 --vsharp-min-radius 1.5 --vsharp-threshold 0.05` (Li, Wu & Liu, 2011; Wu et al., 2012) |
| Inversion | MEDI (Liu et al., 2011), `--medi-lambda <λ from B.4> --medi-merit false --medi-data-weighting 1 --medi-percentage 0.3 --medi-max-iter 30 --medi-cg-max-iter 10 --medi-tol 0.1`; **SMV off** through the TOML key read from `qsmxt init` (QSM.rs defaults SMV to on, 5 mm); no CSF term (MEDI+0 off; Liu Z et al., 2018). QSM.rs percentage 0.3 is the share of gradient components treated as edges, whereas Cornell/SEPIA's 0.9 is counted over the three components together, so the two settings correspond (QSM reviewer's reading of the source; confirmed at install from the QSM.rs `medi.rs` documentation and logged in amendment 1) |
| Reference | `--qsm-reference none`; the two-step whole-brain reference on S is applied by `qsm_plausibility.py` and `extract_rois.py` (§12.1) |
| R2\* | `--do-r2starmap --do-t2starmap` (ARLO; Pei et al., 2015), clipped downstream to [0, 200] s⁻¹ |
| Run | `qsmxt run <coded BIDS> <out> --config qsm_primary.toml --output-space acquired --n-procs 16` |

**B.3 Variants.**
- **Split-echo:** echo-subset BIDS copies {1,3} and {2,4} (ΔTE 15 ms; echoes renumbered); same parameters; no R2\*. **λ scaling (verified from source, QSM.rs v0.38.0 `pipeline/inversion.rs`):** the ppm field is converted to radians with the **first echo time of the input set** before MEDI, so {2,4} (TE1 = 15 ms) carries twice the radian scale of the primary; because the data term is quadratic and the regulariser L1 in the radian-scale χ, the equivalent regularisation needs **λ_{2,4} = 2λ**; {1,3} (TE1 = 7.5 ms) uses λ. They give PC4, PC4-DN and the split-echo co-report.
- **Multiverse:** background removal PDF and LBV (tolerances set explicitly to the values in the captured `qsmxt init` file; the config stage aborts if a key is ambiguous); inversion iLSQR (Li et al., 2015) and RTS (Kames et al., 2018) at the recorded defaults; all 9 combinations share unwrapping and Mask 3; they feed G3 and the specification curve on their common support.
- **G2 subset runs** on `gate20.txt`: the primary chain (and, if A2 is declined, V-SHARP+RTS and PDF+iLSQR).

**B.4 Pilot, sign check and λ rule** (registered as rules; outputs in amendment 1; `pilot_checks.py`).
- **Pilot:** the 10 sessions in `pilot10.txt` (sign check, mask and slab checks, G5 VDC calibration); the first 2 are the L-curve sessions.
- **λ grids** (strongest to weakest regularisation): QSMxT {0.004, 0.002, 0.001, 0.0005} (QSM.rs λ weights the regulariser on a radian-scale field, roughly the inverse of Cornell's scale); SEPIA/Cornell MEDI {250, 500, 1000, 2000}; engine defaults 0.001 and 1000.
- **Norms** per session, engine and λ: ρ = ‖M(D∗χ − f_local)‖₂ (ppm; D the dipole kernel in k-space with the B0 direction from the affine; f_local the reconstruction's local field in ppm; M its final mask) and η = ‖M∇χ‖₁ (forward differences, ppm per voxel).
- **Corner rule:** x = log₁₀ρ, y = log₁₀η over the 4 grid points, each min–max normalised to [0, 1]; the session's corner is the grid point nearest the origin; **flat curve** (largest perpendicular distance of the two interior points from the chord < 0.10) → the engine default; the 2 sessions are combined by the geometric mean of their choices snapped to the nearest grid value on the log scale, ties toward the engine default. The rule uses only residual and regularisation norms.
- **Sign check** and its branches: Q24.

**B.5 Gate G2 (engine comparison; only if A2 is approved).** SEPIA v1.3.0.1 (MIT; commit recorded; MATLAB 2024a; Chan & Marques, 2021) with the Cornell MEDI toolbox (archive SHA-256 recorded). SEPIA takes the **QSMxT total field** (after phase-offset removal, ROMEO and echo combination; in Hz with CF = MagneticFieldStrength × 42.576 × 10⁶) and the same Mask 3, and runs its native V-SHARP (radius 12:−1.5:1.5 mm, threshold 0.05) and MEDI with isSMV = false, merit = false, isLambdaCSF = false, wData = 1, percentage 0.9, zeropad 0 and λ from its own L-curve. `sepia_header.mat` is written explicitly (CF in Hz; B0 = CF/42.57747892 × 10⁶ T; TE; delta_TE; B0_dir from the affine, computed by the same code as QSMxT's). G2 therefore tests background removal and inversion, not unwrapping. Both engines use the same ROI warps and referencing. The SEPIA runner (`sepia_g2.m`, `sepia_config.m`) implements exactly these settings and is hashed in amendment 1 before any SEPIA χ exists.

**B.6 Registration and ROIs** (T1-only, χ-independent; `register_rois.sh`).
- **MEGRE → same-session T1w:** N4-corrected echo-1 magnitude, ANTs 2.5.3 (Avants et al., 2011) rigid 6 dof, Mattes MI, centre-of-mass initialisation, fixed-image mask = dilated SynthSeg brain.
- **T1w → MNI152NLin2009cAsym: route R for every session** — own `antsRegistrationSyN.sh -t s` of the SynthSeg-masked T1w (the QSIPrep ACPC T1w where available, otherwise the raw T1w) to the TemplateFlow 2009cAsym brain. The QSIPrep transforms are not used on the MEGRE side (no inverse exists; the coded copies are translated).
- **Atlases to native space** by the composed inverse in one `antsApplyTransforms` call each: linear interpolation for probabilities, GenericLabel for labels. **Atlas files:** CIT168 v1.1.0 (CC BY 4.0; sha256 `7a91c367…f66`); **SUIT** `atl-Anatom_space-MNI_probseg` (CC BY; D8: no public release of the He et al. (2017) dentate atlas [UNVERIFIED reference] was found on 2026-09-25, so SUIT is primary; that atlas enters as a dentate sensitivity only if its authors provide it before the freeze); JHU-ICBM-DTI-81 (FSL 6.0.7.11; PLIC for description only).
- **SUIT and JHU (NLin6Asym) to 2009cAsym:** one ANTs SyN between the TemplateFlow `tpl-MNI152NLin6Asym_res-01_T1w` and `tpl-MNI152NLin2009cAsym_res-01_T1w` brains (already on the server; copied read-only; computed once and hashed; no download).
- **SynthSeg** (Billot et al., 2023; `mri_synthseg --robust` in qsiprep-26.0.0) on the echo-1 magnitude for G5. The CSF buffer's label 24 exists only in SynthSeg 2.0; the SynthSeg version of the QSIPrep `desc-aseg` is checked at install; with 1.0 the interpeduncular cistern is not buffered (declared).
- **Head orientation:** from the rotation of the composed MEGRE→T1w rigid and the affine stage of the T1w→2009cAsym registration (one frame for all sessions).
- **Coverage:** the MEGRE field of view is warped to the T1w grid (brain-in-FOV fraction) and to 2009cAsym; the region inside the FOV in ≥ 95% of registered sessions (label-free) is warped back to each MEGRE grid for the coverage-harmonised reference.
- **Registration checks at SN scale.** G5 Dice on brainstem/VDC and the placement rating cannot detect 1–2 mm SN or dentate misplacement (the SN is 3–4 voxels thick; T1w shows no dentate), and the dentate carries heavy partial volume (declared). An aggregate, label-free midbrain boundary-distance metric (mean distance between the warped CIT168 SN–VTA boundary and the echo-1-magnitude SynthSeg ventral-DC boundary) is reported in amendment 2.

**B.7 DWI side** (existing outputs; read-only; `dwi_geometry.py`, `pc5_dwi.py`, `dwiqc.py`).
- `dwi_geometry.py --stage geometry` runs **before posting** over all DWI sessions with an SCP_R bundle (about 20–40 CPU-h): bundle orientation and resampling, node centroids, the dentate and SCP-intermediate-ROI points mapped from template space to ACPC with `antsApplyTransformsToPoints` (with a built-in transform-direction check: 7 template brainstem points must fall in the session's SynthSeg brainstem in ≥ 5 cases), the rule (i)–(iii) flags, the NODDI profiles, the N4 b0 profile, the segment values, the peak alignment and the split-half reliability; it fixes n_min on the pre-QC candidate list. It outputs real-ID node and session tables only, with no volumes and no χ.
- After posting: `dwiqc.py` (coded bundle QC; setup and freeze), `dwi_geometry.py --stage t1t2` and `--stage table` (the frozen S3 table with `s3_eligible`), and `pc5_dwi.py` (reads only participant/session, age and sex columns from the covariate file it is given).
- **Fallback:** if the geometry job cannot run before posting, rule (ii) becomes a 3-mm arc-length proxy, b0loc is unavailable, and **S3 is registered as exploratory** (Holm over S1 and S2); the posted version states which applies.

**B.8 Other derived measures.** Q_phase, the motion index, the reference-domain fractions (§12.2); T1w/T2w along pSCP_R only (real-ID DWI side; the SN–VTA T1w/T2w of the draft is dropped); χ-separation (X2): Dimov closed form (Q25); **ΔR2\* (tier 0b; D10)** only if download A3-5 (ds007088, about 1.7 GB) is approved before posting, otherwise dropped with the reason stated (BOLD motion-corrected with the ds007088 transforms, FD > 0.3 mm censored, ΔR2\*_STR within the SynthSeg striatum following Sullivan-Toole et al., 2026, preprint; details frozen in amendment 1 before any value exists [UNVERIFIED details]); **CBF (X7):** ds006744 `atlas-4S456Parcels_cbf.tsv` with PennLINC QEI QC; **FC (X5):** ds007402 4S456 time series (parcels 407/421 SN–VTA; 447–456 cerebellar regions [identity to verify]).

**B.9 Compute and approvals.** About 185–365 CPU-h in total and 25–110 GB of output (QSMxT writes uncompressed NIfTI; measured on the pilot), all inside the tree (A4). Required downloads: about 0.01 GB of atlases (A3-1/2) and the QSMxT tarball (A1); optional A3-3 to A3-7. A2 (SEPIA + Cornell MEDI; the user completes the Cornell registration form) is optional. No R or Python package is installed (A11 is not needed).

---

## Annex C. Deviations policy
1. **Registered branches are not deviations.** Gate outcomes, fallbacks, the λ rule, the G5 VDC calibration, the sign-check branches and the D10/A2/A3 conditionals are reported as "branch taken".
2. **Stage I (after posting, before the first χ map).** A change is allowed only for technical necessity (a missing flag or file, a binary that does not run, an install-time `VERIFY` item that does not hold). The nearest documented equivalent is used, or the next QSMxT tag with no change to any registered algorithm parameter; it is logged in a time-stamped `DEVIATIONS.md` and posted in amendment 1, before the full run.
3. **Stage II (χ exists, before unlock).** A change is allowed only when a registered rule requires it or when a pre-named label-free metric reveals a technical failure (Q24). Its rationale may cite only those metrics; it applies to all sessions; it is posted in amendment 2 with before/after aggregates. **No parameter is tuned on χ.**
4. **Stage III (after unlock) and coding errors.** Registered confirmatory analyses do not change. A coding error — any deviation of the hashed code from the registered text, found at any stage — is corrected; both the registered-code and the corrected results are shown, and the text-conforming result is decisive. Any new analysis is labelled post hoc.
5. **Blinding incidents.** An accidental pre-unlock exposure (for example, a per-code value seen next to an identifying detail, or an `audit_run.py` rejection) is logged immediately and described in amendment 2 and in the paper; the confirmatory status of the affected analyses is reassessed and stated.
6. **Reporting.** The paper includes a deviations table (stage, change, reason, impact). Schedule slippage is not a deviation.

## Annex D. Timeline
Week 0 starts at the user's approvals [FILL: date]. One researcher; 14 weeks.

| Week | Work | Deliverable |
|---|---|---|
| 0 | Approvals (A1, A3–A6, A8, A12; A2 optional); IRB determination request; lock-key generation; deletion or sealing of the round-1 MEGRE files (Q18) | Signed-off plan |
| 1 | Tree and environment; header-only audit; tool checks; QSMxT install and `--help`/`init` capture; atlas downloads; template-to-template transform; `dwi_geometry.py` and **n_min**; `lock_labels.py`; `blind_setup.py`; hashing; **OSF posting**, before any χ | Registration timestamp |
| 2 | Registrations; geometry; Mask 3; pilot (sign, λ per engine, G5 calibration); G2 runs and `gates_preunlock.py g2`; **amendment 1**; full QSM launch | Amendment 1 |
| 3 | Full QSM, R2\*, split-echo and multiverse; plausibility; ROI and multiverse tables; coded DWI bundle QC; ΔR2\* (D10) | Coded tables |
| 4 | QC pass 1; DWI-QC freeze; DWI table; PC5 | QC pass 1 |
| 5–6 | 20% re-rating (≥ 14 days later); κ; usable lists; `gates_preunlock.py final`; freeze; **amendment 2** | Freeze manifest |
| 7 | **Unlock** → the confirmatory sequence (A.6) and the verification | Confirmatory results |
| 8 | Amendment-2 modules: E-dev, E-z, X1–X7, node-level profile, DWI multiverse | Robustness report |
| 9–12 | Figures and manuscript | Draft |
| 13 | Novelty searches re-run; preprint; submission (*BP:CNNI* default) | Submission |
| 14 | Buffer and code audit | — |

---

## Annex E. Attachments and hashes (SHA-256)

All files below are uploaded to the embargoed OSF project with the registration (file list and upload order: `USER_ACTIONS_TR.md`). Hashes computed on 2026-09-29 (UTC) with `shasum -a 256`; `MANIFEST.sha256` lists every file of the `prereg/` folder. The SHA-256 of this document itself is quoted in the OSF description (a file cannot contain its own hash).

| Attachment | Path | SHA-256 |
|---|---|---|
| CHANGELOG (draft → FINAL) | `CHANGELOG.md` | `889d5073faa301742f076be590c115b78199e0fe5a036e2de963f3367c7ff696` |
| Review R_editor (R1–R3) | `review/R_editor.md` | `84e66a04f78b0fe752a9085e710ebeb94db31e65cc3ce83864c3ca0d18a767ba` |
| Review R_qsm (R1–R3) | `review/R_qsm.md` | `6e59c0ee5e4fb79d20a6e71eac44f846f43b691cfebc0618cef7cb2b7d0f6a94` |
| Review R_stats (R1–R3) | `review/R_stats.md` | `2539e04514768adffa97bc80ea3ad811118d961e19c5370044c6023751ef7834` |
| Draft v1.0 (superseded; for the record) | `LEADchi_prereg_v1.0.md` | `ae83df35638d86926112ae61c8618e7ac530676935ab830c8e59497dd8563a76` |
| Analysis code: SHA256SUMS (hashes every file of code/analysis/) | `code/analysis/SHA256SUMS` | `f1231bbc66689a74f759d048e535fcf25e32382ac35962b12830b2f93032bade` |
| Pipeline code: SHA256SUMS (hashes every file of code/pipeline/) | `code/pipeline/SHA256SUMS` | `f34e9f80f26665e3a5cd21fe98b60294e75ff5b2db92e7ec508d18c550e0a64f` |
|   analysis/run_all.R | `code/analysis/run_all.R` | `435b86b82d36dba653a1f9c440352843e1d539766501af1837e0e03b15bde184` |
|   analysis/R/m1.R | `code/analysis/R/m1.R` | `75173d15662564a8a8a12ebecd989cfd90100742f6ce5cb2ca2fa2046c78257d` |
|   analysis/R/spec_curve.R | `code/analysis/R/spec_curve.R` | `ee00f202bcbfa262058061c269757eda21061964a90a0ab113287a35f33923a7` |
|   analysis/verify_statsmodels.py | `code/analysis/verify_statsmodels.py` | `7b1b08a9d01f1a6a564f8ccdbd6180f64fdff62695e8f4aec9fbcc7e4f0f43c4` |
|   analysis/sim/typeI_power.md | `code/analysis/sim/typeI_power.md` | `fec73003bcfd0296947d1d379d12da50db1d89df8882df545964eed7c56b0ee6` |
|   pipeline/gates_preunlock.py | `code/pipeline/gates_preunlock.py` | `4fbf220cbdcce0009a4c9d37e9b956e37da74cb7a1ba7e5471700297a785beb3` |
|   pipeline/pilot_checks.py | `code/pipeline/pilot_checks.py` | `bd11f9112094829e156270cbed5cdf5a3d1d77863d925bd0eaf2a6c9b6db2b84` |
|   pipeline/audit_run.py | `code/pipeline/audit_run.py` | `528a4b28959cc3ec43a5211ab7668a69f145773c5c6ee08621e9294dd8e625fe` |
|   pipeline/dwi_geometry.py | `code/pipeline/dwi_geometry.py` | `d5ab6b43afc43a8f070fedc771f8ceeab5b4d37fa01d410eb4036a8d2a0eabc8` |
|   pipeline/blind_setup.py | `code/pipeline/blind_setup.py` | `76923c530adb338a57d01e79d73d22e39703df4f13134caac382e17dfbf1cd6a` |
|   pipeline/run_qsmxt.sh | `code/pipeline/run_qsmxt.sh` | `ef8144ed956e1aff15cc8bdd530d752775197b43373160776f81f064df34fa06` |
|   pipeline/s3_candidates_preqc.csv | `code/pipeline/s3_candidates_preqc.csv` | `053c31b5a4e2a88ffa7ba7cd22766211fb12f64e107e31c27bb1b4c83b273c6c` |
| E1/E4 exposure audit | `work/E1_exposure.md` | `2e9128222e7b8473fe2abb1a449ee3a061fd90521d14a1054b98527f6c8189b4` |
| E2 v1.2 change set | `work/E2_protocol_v1.2_delta.md` | `4bf0e9308754bbd10492c5202ebc169394a01bc3158e1ec1c48e77dbbabca48b` |
| E3 technical pinning | `work/E3_technical.md` | `255257174e6bff7020d67ae12db47d71a9918afff3482d4d4958130b6f512488` |
| E2 output e2_design_v12.txt | `work/e2_design_v12.txt` | `d8cc23f242ac113b72e133b3d7e27c17279a3137ee011aaac1fe2ce23620f7da` |
| E2 output e2_hc2_vs_hc3.txt | `work/e2_hc2_vs_hc3.txt` | `b0225ed4b7de3a367f1fa45573add13c37e9c083c20071f4a0fd3b9182ec2e9f` |
| E2 output e2_hc3_afterqc.txt | `work/e2_hc3_afterqc.txt` | `5893937f1f3b1b8e1d5908ea92cffb574e8355271920130acbfe2bee3aeccca1` |
| E2 output e2_hc3_power_check.txt | `work/e2_hc3_power_check.txt` | `8a3f3ce72678a24b0947ba2b25b16baa6a3b8c7b903c43ae7452cf2335869c83` |
| E2 output e2_labels_hc3.txt | `work/e2_labels_hc3.txt` | `7591a0514f4371da8cd14c3a15520951be5cbe2763cbcdff6b077eab4e5b7621` |
| E2 output e2_lit_ranges.txt | `work/e2_lit_ranges.txt` | `a67aa6bee653c8732f1b86efebb6334b2f1c462db5604e7db2e0cf660d9203cb` |
| E2 output e2_pc2w_states.txt | `work/e2_pc2w_states.txt` | `33d3a6c5ac6a0269267ec01c7c3e01e75a0122e7aff8168c37f06d32ec0b7a4f` |
| E2 output e2_power_v2_rerun.txt | `work/e2_power_v2_rerun.txt` | `d47e0212fdee902140acbd90bf4729d18ff1f4f754cd380d5dcf8c40c908db60` |
| E2 output e2_robust_bound_sim.txt | `work/e2_robust_bound_sim.txt` | `24d957914f27ba6c717461ba75950bd57c403a3b97bef6e4e7a314beb7cdcf7a` |
| E2 output e2_s1_hc3.txt | `work/e2_s1_hc3.txt` | `b535ea6fc382570369835207a8c9a2450ff15e76f296224139690dcc34dcd4cc` |
| E2 output e2_scpr_reliability.csv | `work/e2_scpr_reliability.csv` | `99cb2443f2855b5ee618c6f634a633fecfc33086055debd84a5eec09a3363243` |
| E2 output e2_scpr_reliability.txt | `work/e2_scpr_reliability.txt` | `13a70752335b5b095ae76b8ff4858e3e4ae751755afb37fd4a7d69474c645dd0` |
| E2 script e2_design_v12.py | `work/e2_scripts/e2_design_v12.py` | `16e45b05f0ba9e47d86feec64fff329183396095f19c5ed6027694ceb9dc52db` |
| E2 script e2_hc2_vs_hc3.py | `work/e2_scripts/e2_hc2_vs_hc3.py` | `0d0ccb956bd12992d996b87a84aab793135a5ec126d34a8c746b16eb4195b349` |
| E2 script e2_hc3_afterqc.py | `work/e2_scripts/e2_hc3_afterqc.py` | `a37f6728eecb58d59d9b3c840c1deed0e8a412e6e471677afc6ad3a328d774f9` |
| E2 script e2_hc3_power_check.py | `work/e2_scripts/e2_hc3_power_check.py` | `151c52835916a39557c6aaeaaa814b2d9c8bc972b56911e2432f8dbc11c461df` |
| E2 script e2_labels_hc3.py | `work/e2_scripts/e2_labels_hc3.py` | `978548240c6f95d80d19225b9b96a8ec1c9f3f877a360351a2ec0cd829790a5d` |
| E2 script e2_lit_ranges.py | `work/e2_scripts/e2_lit_ranges.py` | `c7e78190b9eec3a16320c26a0307567109fc8d012f15d5dc95a7d9a580c167ee` |
| E2 script e2_pc2w_states.py | `work/e2_scripts/e2_pc2w_states.py` | `9d7d3ff660586aaab4a494ac216f9f0afbe9ded32244de8166ed6b02f63397c3` |
| E2 script e2_power_v2_rerun.py | `work/e2_scripts/e2_power_v2_rerun.py` | `82cd3872cf204fe9713477b5ac1d25aa497ebfd38fd6edf2c368a699279ab8ad` |
| E2 script e2_robust_bound_sim.py | `work/e2_scripts/e2_robust_bound_sim.py` | `043dfebd44339c3a86a2be6a6e1e0c9c4ad352a83d262677111a55d90d505242` |
| E2 script e2_s1_hc3.py | `work/e2_scripts/e2_s1_hc3.py` | `7b1e258bd117f05b463f2e5320be647a01f290f472bfed29e661508b3ffa24b7` |
| E2 script e2_scpr_reliability.py | `work/e2_scripts/e2_scpr_reliability.py` | `1b5eaa73c38a8a738a34f29399f3af828533a64f2e94c0b1a52d2b3829bb7058` |
| Protocol v1.1 | `claude_rapor/1_tasarim/LEAD-chi_final_protokol_v1.1.md` | `1de49dd35aaf29de8f9af3b2253c0fef028f25219139d531d38b480dc48bcb1b` |
| Paper 0 Paper0_analysis_plan_v1.0.md | `paper0/plan/Paper0_analysis_plan_v1.0.md` | `f76e1051aedcf027954a5b15e39129437e9d5a311093d4183ec9def9eecafe13` |
| Paper 0 Paper0_amendment_v1.1.md | `paper0/plan/Paper0_amendment_v1.1.md` | `0c145c133aa823145dc015bf5f1f2179728ab8821e6ccbd8d44d4f4ea3f69a06` |
| Paper 0 Paper0_amendment_v1.2.md | `paper0/plan/Paper0_amendment_v1.2.md` | `7aee5fe006a1b0d2c63e0b95f42a63c4417322af95eea8720864b54f32dcc59d` |
| Paper 0 Paper0_amendment_v1.3_posthoc.md | `paper0/plan/Paper0_amendment_v1.3_posthoc.md` | `484b7c51fe230681ba77515c213f343517506becbcbccb97cb1784be6cffdc91` |
| Paper 0 draft v4 | `paper0/manuscript/Paper0_draft_v4.md` | `9402b92fedc83a02f0682e04f85f1501053c7b78f4bd0fe3e4a19529fbe0d0b8` |
| Power v1.1 script | `claude_rapor/4_analiz_kodlari/guc_analizleri/v2_work/power_v2.py` | `526517691240fdd05de9ed4ef8e0f009561492d66e56f0391afa1a30c2962767` |
| Power v1.1 output | `claude_rapor/4_analiz_kodlari/guc_analizleri/v2_work/power_v2.txt` | `36822603802597a7e3d34edf7be7e287d8ef5b70124d6af49955ae3310edd27a` |

**Pinned external items and exposure evidence** (not attachments; values fixed earlier):

| Item | SHA-256 |
|---|---|
| Exposure evidence (removed or to be removed): 8 sub-20139 ses-1 NIfTI | as listed in E4 §4 (`4683c78c…e3c9` … `ab55c47f…71f6c6`) |
| Exposure evidence: `megre_sub-20139_coverage.png` / `_phasecheck.png` / `megre_coverage.tsv` (values recorded in the draft; not re-read) | `7fb414002bee3e5737b518d728af9a6e1652305c0e15948ff2c7c3f740605fe0` / `090c82c4a473655b0cf33a5bbb4b790a21327e7306204405006d91688ae82ef2` / `6ee78f5a04c9c810d8f10e1b6ef06f2e46432580b908891c78c23e272add8ea8` |
| QSMxT v9.22.0 musl tarball (pinned; downloaded at A1) | `6593438c029749b5e3fb5c78ce485f6f5e282545b1892832cf30539b0326ca49` |
| CIT168 `CIT168toMNI152-2009c_prob.nii.gz` / `labels.txt` (pinned; downloaded at A3) | `7a91c36701f36d9d908479b862e65a879d15f54badcec6878e61d857ab570f66` / `5e00265a7bbbdf242356345e0f6615d84791ada54aa9739539096aab3c1e1b41` |

**To be filled by the finaliser before posting** (each item gets its SHA-256 or its value):
- **Blinding:** the `lock_pub.pem` fingerprint; `key.csv.enc`, `key.sessionkey.enc` (from `blind/setup_manifest.json`); `labels_locked.tar.enc`, `labels_locked.sessionkey.enc`, `labels_locked_manifest.json`; `pilot10.txt`, `lcurve2.txt`, `gate20.txt`; `coded_files.sha256`; the DWI-QC key file (`dwiqc/dwiqc_key.csv`, if set up before posting; otherwise amendment 1).
- **Install captures:** `qsmxt --version`; `qsmxt_run_help.txt`; `qsmxt_v9.22.0_defaults.toml`; `qsmxt_quality_map_help.txt`; `qsm_primary.toml`; `pipeline_choices.json`; `qsmxt_flags_verified.txt`; the qsiprep-26.0.0 container digest.
- **ROI and transform files:** the SUIT probseg and its label table; the JHU atlas; the NLin6Asym→2009cAsym transform; the two TemplateFlow T1w templates; the pyAFQ SCP_R intermediate ROI (A5).
- **DWI geometry:** `dwi/n_min.json` (and the value of n_min, §12.3); `dwi/dwi_sessions.csv`; `dwi/dwi_nodes.csv`.
- **Header audit:** the aggregate JSON of `blind_setup.py --header-audit` (scaling fields, obliquity, sequence parameters, visible-channel uniqueness).
- **Logs:** the audit-log and pipeline-log head hashes at posting.
- **Public-file hashes** (from `labels_locked_manifest.json`): `participants.tsv`, `sessions.tsv`, phenotype TSVs; the PennLINC QC CSVs.

## Annex F. References
All entries were verified via PubMed, Crossref or DataCite on 2026-09-23/26 (protocol v1.1 §28; E2 §3.4; E3; task F), unless marked.
- Avants BB et al. *NeuroImage* 2011;54:2033–2044. doi:10.1016/j.neuroimage.2010.09.025. PMID 20851191.
- Bilgic B et al. (QSM Consensus Organization Committee). *Magn Reson Med* 2024;91:1834–1862. doi:10.1002/mrm.30006. PMID 38247051.
- Billot B et al. *Med Image Anal* 2023;86:102789. doi:10.1016/j.media.2023.102789. PMID 36857946.
- Button KS et al. *Nat Rev Neurosci* 2013;14:365–376. doi:10.1038/nrn3475. PMID 23571845.
- Chan KS, Marques JP. *NeuroImage* 2021;227:117611. doi:10.1016/j.neuroimage.2020.117611. PMID 33309901.
- Diedrichsen J et al. *NeuroImage* 2011;54:1786–1794. doi:10.1016/j.neuroimage.2010.10.035. PMID 20965257.
- Dienes Z. *Front Psychol* 2014;5:781. doi:10.3389/fpsyg.2014.00781. PMID 25120503.
- Dimov AV et al. *Tomography* 2022;8:1544–1551. doi:10.3390/tomography8030127. PMID 35736875.
- Dymerska B et al. *Magn Reson Med* 2021;85:2294–2308. doi:10.1002/mrm.28563. PMID 33104278.
- Eckstein K et al. *Magn Reson Med* 2018;79:2996–3006. doi:10.1002/mrm.26963. PMID 29034511.
- Hoffman LJ et al. *J Neurosci* 2024;44:e1031242024. doi:10.1523/JNEUROSCI.1031-24.2024. PMID 39256045.
- Hoopes A et al. *NeuroImage* 2022;260:119474. doi:10.1016/j.neuroimage.2022.119474. PMID 35842095.
- Jain V et al. *Magn Reson Med* 2012;68:863–867. doi:10.1002/mrm.23282. PMID 22162033.
- Kames C et al. *NeuroImage* 2018;167:276–283. doi:10.1016/j.neuroimage.2017.11.018. PMID 29138089.
- Kobayashi H et al. *Schizophr Res* 2008;106:356–362. doi:10.1016/j.schres.2008.08.018. PMID 18809299.
- Koo TK, Li MY. *J Chiropr Med* 2016;15:155–163. doi:10.1016/j.jcm.2016.02.012. PMID 27330520.
- Lakens D. *Soc Psychol Personal Sci* 2017;8:355–362. doi:10.1177/1948550617697177. PMID 28736600.
- Li F, Morgan KL, Zaslavsky AM. *J Am Stat Assoc* 2018;113:390–400. doi:10.1080/01621459.2016.1260466.
- Li W, Wu B, Batrachenko A, et al. *Hum Brain Mapp* 2014;35:2698–2713. doi:10.1002/hbm.22360. PMID 24038837.
- Li W, Wu B, Liu C. *NeuroImage* 2011;55:1645–1656. doi:10.1016/j.neuroimage.2010.11.088. PMID 21224002.
- Li W et al. *NeuroImage* 2015;108:111–122. doi:10.1016/j.neuroimage.2014.12.043. PMID 25536496.
- Linscott RJ, van Os J. *Psychol Med* 2013;43:1133–1149. doi:10.1017/S0033291712001626. PMID 22850401.
- Liu T et al. *Magn Reson Med* 2011;66:777–783. doi:10.1002/mrm.22816. PMID 21465541.
- Liu Z et al. *Magn Reson Med* 2018;79:2795–2803. doi:10.1002/mrm.26946. PMID 29023982.
- McGraw KO, Wong SP. *Psychol Methods* 1996;1:30–46. doi:10.1037/1082-989X.1.1.30.
- Morandini HAE et al. *J Affect Disord* 2026;402:121329. doi:10.1016/j.jad.2026.121329. PMID 41653994.
- Naji N et al. *NMR Biomed* 2022;35:e4788. doi:10.1002/nbm.4788. PMID 35704837.
- Neuhaus JM, Kalbfleisch JD. *Biometrics* 1998;54:638–645. doi:10.2307/3109770.
- Otsu N. *IEEE Trans Syst Man Cybern* 1979;9:62–66. doi:10.1109/TSMC.1979.4310076.
- Pacchiano F et al. *Pediatr Radiol* 2026;56:1268–1284. doi:10.1007/s00247-026-06565-7. PMID 41801366.
- Parker CS et al. *NeuroImage* 2021;245:118749. doi:10.1016/j.neuroimage.2021.118749. PMID 34852276.
- Pauli WM, Nili AN, Tyszka JM. *Sci Data* 2018;5:180063. doi:10.1038/sdata.2018.63. PMID 29664465.
- Pei M et al. *Magn Reson Med* 2015;73:843–850. doi:10.1002/mrm.25137. PMID 24664497.
- Salazar de Pablo G et al. *JAMA Psychiatry* 2021;78:970–978. doi:10.1001/jamapsychiatry.2021.0830. PMID 34259821.
- Salman F et al. *Hum Brain Mapp* 2025;46:e70187. doi:10.1002/hbm.70187. PMID 40260740.
- Schulze M et al. *Biol Psychiatry Cogn Neurosci Neuroimaging* 2025;10:597–606. doi:10.1016/j.bpsc.2024.08.015. PMID 39218036.
- Sevchik BL et al. *bioRxiv* 2025. doi:10.1101/2025.11.10.687633 (preprint).
- Sevchik BL et al. *Sci Data* 2026. doi:10.1038/s41597-026-08095-1.
- Sevchik BL et al. Penn LEAD [dataset]. OpenNeuro ds007116 v1.0.6. doi:10.18112/openneuro.ds007116.v1.0.6.
- Shrout PE, Fleiss JL. *Psychol Bull* 1979;86:420–428. doi:10.1037/0033-2909.86.2.420.
- Simonsohn U, Simmons JP, Nelson LD. *Nat Hum Behav* 2020;4:1208–1214. doi:10.1038/s41562-020-0912-z. PMID 32719546.
- Sørensen Ø, Walhovd KB, Fjell AM. *NeuroImage* 2021;226:117596. doi:10.1016/j.neuroimage.2020.117596. PMID 33248257.
- Stewart AW et al. *Magn Reson Med* 2022;87:1289–1300. doi:10.1002/mrm.29048. PMID 34687073.
- Sullivan-Toole H et al. *bioRxiv* 2026. doi:10.64898/2026.04.10.717850. PMID 42039616 (preprint).
- Tang S et al. *Hum Brain Mapp* 2022;43:2495–2502. doi:10.1002/hbm.25798. PMID 35107194.
- Treit S et al. *Hum Brain Mapp* 2021;42:4597–4610. doi:10.1002/hbm.25569. PMID 34184808.
- Tubiolo PN et al. *Biol Psychiatry Cogn Neurosci Neuroimaging* 2026. doi:10.1016/j.bpsc.2026.09.002. PMID 42716338.
- Van den Akker OR et al. *Meta-Psychology* 2021;5:MP.2020.2625. doi:10.15626/MP.2020.2625.
- Vano LJ et al. *Am J Psychiatry* 2025;182:830–839. doi:10.1176/appi.ajp.20240512. PMID 40887951.
- Vano LJ et al. *Mol Psychiatry* 2026;31:941–952. doi:10.1038/s41380-025-03195-7. PMID 40913112.
- Vano LJ et al. *Am J Psychiatry* 2026. doi:10.1176/appi.ajp.20251257. PMID 42265595.
- Velioglu HA et al. *Biol Psychiatry Cogn Neurosci Neuroimaging* 2026. doi:10.1016/j.bpsc.2026.04.001. PMID 41966263.
- Velioglu HA et al. Paper 0 (Cerebello–ventral tegmental pathway diffusivity in psychosis-spectrum youth), manuscript draft v4, 2026 [UNVERIFIED as a publication; hashes in Annex E].
- Winkler AM et al. *NeuroImage* 2014;92:381–397. doi:10.1016/j.neuroimage.2014.01.060. PMID 24530839.
- Wood SN. *J R Stat Soc B* 2011;73:3–36. doi:10.1111/j.1467-9868.2010.00749.x.
- Wu B et al. *Magn Reson Med* 2012;67:137–147. doi:10.1002/mrm.23000. PMID 21671269.
- Zhang Y et al. *NeuroImage* 2018;171:176–189. doi:10.1016/j.neuroimage.2018.01.008. PMID 29325780.
- *Software sources consulted (documentation/source pages, read 2026-09-25/26):* QSMxT v9.22.0 `src/cli.rs` and `src/pipeline/runner.rs`; QSM.rs v0.38.0 `src/inversion/medi.rs` and `src/pipeline/inversion.rs` (github.com/QSMxT/QSMxT; github.com/astewartau/QSM.rs).

## Annex G. Finaliser checklist (complete before posting; nothing else changes)
1. Fill every `[FILL: …]`: authors, ORCIDs, co-author access and knowledge (Q8, Q17, Q18), embargo date, IRB status (Q9), the week-0 date, the deletion or sealing dates (Q18), the header-audit values (Q5, A.4), n_min (§12.3), and the Annex E install and blinding hashes.
2. **Install-time checks** (every `VERIFY` marker in the hashed code is one of these; an item that does not hold is a Stage I technical deviation, Annex C.2): every QSMxT flag of B.2 in the captured `--help`; the MEDI SMV-off TOML key and the PDF/LBV tolerance keys in `init`; the `qsmxt quality-map` command; the `--use-custom-masks` glob with `--dry`, and the post-run mask assertion on the pilot; the MEDI percentage correspondence (B.2); the QSMxT output names; the SUIT probseg layout and the JHU PLIC label values; the pyAFQ bundle and figure names; the phase-encode axis and the partition axis; the SynthSeg version of the QSIPrep aseg; the `dwi_geometry.py` transform check on one session.
3. Settle from the header-only audit: the NIfTI phase `scl_slope`/`scl_inter`, whether `SoftwareVersions` varies, the obliquity distribution, the sequence parameters and the visible-channel uniqueness count.
4. Record whether A2 (SEPIA + Cornell MEDI), A3-4 (AutoTrack), A3-5 (ΔR2\*), A3-6 (X7) and A3-7 (X5) were approved. Mark the text of any conditional analysis that will not run as "not run".
5. Re-verify the Li 2014, Zhang 2018 and Treit 2021 values behind §12.4 against the PDFs (they were extracted automatically).
6. Check each Q heading against the live OSF form and paste the answers into the matching fields. Title and authors go in the metadata; Part 6 and the embargo go in the description.
7. D14: if Paper 0 is posted first as a preprint or OSF deposit, replace its [UNVERIFIED] citation with the DOI.
8. Compute the SHA-256 of the final version of this file and quote it in the OSF description; re-run `shasum -a 256 -c MANIFEST.sha256` on every changed file.

## Annex H. Rule → code crosswalk (registered)

| Registered rule | File : function |
|---|---|
| Table validation; coding; samples (first usable, S2, S3, T1, ses-1-only, exposure sets) | `analysis/R/design.R` : `validate_table`, `code_vars`, `set_race3`, `build_samples` |
| Gate-JSON validation; branch; caveats; environment pin | `analysis/run_all.R`; `analysis/R/constants.R` (`gate_required`, `pinned`) |
| Plausibility re-verification; CSF validity; robust z | `analysis/R/plausibility.R` : `apply_plausibility`, `robust_z_counts` |
| M1, HC3, t reference | `analysis/R/m1.R` : `m1_formula`, `fit_ols_hc3`, `term_test`, `lin_test` |
| σ_B (Huber k 1.345, maxit 200, MAD fallback) | `analysis/R/plausibility.R` : `huber_fit`; `m1.R` : `h1_decide` |
| A/B/C, δ_B, B evaluability (PC4, r̂_L ≥ 0.25), disconfirmation | `analysis/R/m1.R` : `h1_decide`, `bound_test`, `run_m1_core` |
| Final outcome (K3 conjunctive rule); K2; K3 inferential trigger | `analysis/R/m1.R` : `run_h1` |
| B co-reports (−0.6, −0.5, cross-session, residualised split-echo) | `analysis/R/m1.R` : `run_h1`, `resid_split_reliability` |
| K1a/K1b/K1c and precedence | `analysis/R/m1.R` : `k1_labels` |
| K4 fragile; influence-dependent | `analysis/R/m1.R` : `run_h1`; `plausibility.R` : `dominance` |
| K4 robust; specification curve; joint permutation | `analysis/R/spec_curve.R` : `spec_grid`, `run_spec_curve` |
| Freedman–Lane (M1, S3); bootstrap δ CI; Bayes factors; noncentral-t CI | `analysis/R/m1.R` : `fl_perm`, `boot_delta`, `bayes_normal`, `nct_ci_delta` |
| Subject-level QC-by-group and IPW; PC2x | `analysis/R/m1.R` : `qc_by_group`, `pc2x` |
| S1; "dissociation" wording | `analysis/R/secondary.R` : `run_s1`, `s1_fit`, `s1_labels` |
| S2 and its sensitivities; S2 all-session model | `analysis/R/secondary.R` : `run_s2`, `run_s2_allsession` |
| S3a, S3b, S3 sensitivities incl. exposure set | `analysis/R/secondary.R` : `run_s3`, `s3_fit`, `s3_labels` |
| Holm with direction requirement; secondary influence | `analysis/R/secondary.R` : `holm_family`, `secondary_influence`; `utils.R` : `holm` |
| T1 (containment df decision; PB and cluster-robust co-reports); REML-maximum rule | `analysis/R/t1.R` : `t1_core`, `run_t1`, `reml_max`, `pb_test`, `wild_cluster`, `cr2_setup` |
| Pooled all-session H1 sensitivity | `analysis/R/t1.R` : `run_pooled` |
| H1 sensitivity battery; PDS; imputation fallback | `analysis/R/sensitivity.R` : `run_h1_sensitivities`, `pds_score`, `impute_pds_once`, `first_relaxed` |
| Monte Carlo precision; resampling p | `analysis/R/utils.R` : `mc_info` |
| Registered sentences; report | `analysis/R/m1.R` : `REGISTERED_SENTENCES`, `outcome_sentence`; `report.R` : `REG_SENT2`, `registered_reading`, `write_report` |
| Verification | `analysis/verify_statsmodels.py` |
| Label lock | `pipeline/lock_labels.py` |
| Coded copy, header scrubbing, translation offsets, key, draws, header audit | `pipeline/blind_setup.py` : `main`, `read_meta`, `header_audit`, `f0_clusters`, `obliquity_deg`; `lx_common.py` : `write_scrubbed`, `hybrid_encrypt` |
| QSMxT pins, variants, λ_{2,4} = 2λ, Mask 2/3, slab trim, sign negation | `pipeline/run_qsmxt.sh` : `build_args`, `variant_spec`, `mask3_one`, `stage_negate`; `lx_common.py` : `otsu_threshold`, `trim_partition_edges` |
| Registration route R; atlas warps; coverage region | `pipeline/register_rois.sh` : `stage_session`, `stage_harmref` |
| Support S; two-step reference; CSF-only reference; tissue fractions; mask assertion | `pipeline/lx_common.py` : `support_mask`, `two_step_reference`, `csf_reference`, `tissue_fractions`, `mask_assertion`, `variant_outputs` |
| ROI definitions and statistics | `pipeline/lx_common.py` : `build_rois`, `roi_stats`; `atlas_config.json` |
| Q_phase; motion index; orientation | `pipeline/lx_common.py` : `q_phase`, `ghost_to_signal`, `orientation_angles`; `extract_rois.py` : `geom_one` |
| Plausibility rules and flags | `pipeline/qsm_plausibility.py` : `one_session`, `stage_flags` |
| ROI tables; multiverse on the common support; G5 flags | `pipeline/extract_rois.py` : `roi_one`, `stage_multiverse`, `mv_one`, `stage_g5`, `g5_apply` |
| Sign check; λ rule; G5 VDC calibration | `pipeline/pilot_checks.py` : `stage_sign`, `stage_lcurve` (`norms`, `corner`, `snap`), `stage_g5cal` |
| G2, G2-DN, PC1–PC3, fallback, G3 (+ SN-only), PC4, r̂_L, PC4-DN, pooled constants, cross-session reliability, SD(Δref), G4, G5, PC5, n_min | `pipeline/gates_preunlock.py` : `stage_g2`, `pc_checks`, `stage_final`, `icc3_1`, `icc2_1`, `pearson_ci`, `deming_slope` |
| Visual QC, κ, stop rule, code D | `pipeline/qc_montage.py` : `cmd_sheet`, `cmd_kappa`, `cmd_finalize`; `qc_rubric.md` |
| First-usable lists; freeze manifest | `pipeline/build_usable.py` : `stage_usable`, `stage_manifest` |
| Audit wrapper | `pipeline/audit_run.py` |
| DWI node rules, segment values, n_min, S3 table; PC5; DWI bundle QC | `pipeline/dwi_geometry.py` : `session_job`, `retained_nodes`, `stage_geometry`, `stage_t1t2`, `stage_table`; `pc5_dwi.py`; `dwiqc.py` |
| Unlock and merge | `pipeline/unlock_and_merge.py` |

## Annex I. Registered implementation rules (formerly the "P5/P5b choices")
- **T1:** decision by the containment-df test (the max(containment, PB) rule of the draft is withdrawn); occasion 1 = the earliest usable session after the first usable one; T1 scale = the pooled SD over the frozen first-usable sessions.
- **Exposure sets:** the four listed rows are dropped with no substitution (H1, S3, T1).
- **QC-by-group and IPW:** subject level (Q13).
- **Bootstraps:** S3b and the δ CI are cell-stratified.
- **Huber scale:** k = 1.345, maxit 200, MAD fallback.
- **S1 stacked GLS:** ROI-specific intercepts and slopes for all covariates, `varIdent(~1|roi)`, `corSymm(~rk|sub)`, REML.
- **Holm:** "supported" only in the predicted direction.
- **Composite envelopes:** STR [−80, 200], DGM [−80, 300], right putamen [−80, 200], dentate + interposed [−50, 300] ppb.
- **Primary-ROI plausibility failure in a row marked usable:** hard stop (the frozen list contradicts (c′)); resolved only as a documented coding error.
- **Pipeline:** B0 4 dp; phase float32 radians; route R; one orientation frame; the 1-voxel SynthSeg-CSF buffer for every atlas ROI; SynthSeg ROIs eroded per label with 6-connectivity, NAcc not eroded; roi_vol atlas-only; implausible voxels excluded from the CSF means; failed non-primary ROIs written as NA; R2\* of every variant from the primary run; slab-coverage flag forces a review only; G5 VDC = 28 ∪ 60, MEGRE↔T1w Dice inside the MEGRE FOV, SynthSeg `--robust` on the echo-1 magnitude; pilot and G2 draws among sessions with a usable T1w, the G2 set disjoint from the pilot with 20 distinct subjects; the L-curve on the first 2 pilot sessions.
- **QC:** linear κ is the G4 statistic; the worse grade and the union of failure codes for re-rated items; random QC IDs; visual grades shared by the fallback and SEPIA branches.
- **Synthetic runs** mask which exposure sessions fall in which cell or sample.

## Annex J. Module status

| Status | Modules |
|---|---|
| **Written, tested on synthetic data and hashed at posting** (every Tier 0–2b quantity, the K4 "robust" label and the H1/S2/S3 sensitivities) | `analysis/`: `run_all.R`, `R/*.R` (constants, utils, design, plausibility, m1, secondary, t1, sensitivity, spec_curve, report), `schema/*.csv`, `verify_statsmodels.py`, the test, synthetic and simulation files. `pipeline/`: `env.sh`, `lx_common.py`, `atlas_config.json`, `lock_labels.py`, `blind_setup.py`, `run_qsmxt.sh`, `register_rois.sh`, `qsm_plausibility.py`, `extract_rois.py`, `qc_montage.py`, `qc_rubric.md`, `build_usable.py`, `pilot_checks.py`, `gates_preunlock.py`, `audit_run.py`, `dwi_geometry.py`, `pc5_dwi.py`, `dwiqc.py`, `unlock_and_merge.py`, `s3_candidates_preqc.csv`, `tests/test_reference.py` |
| **Amendment 1** (before any χ they touch exists; conditional) | `sepia_g2.m` + `sepia_config.m` (only if A2; implement B.5 verbatim); the D10 ΔR2\* job (only if A3-5) |
| **Amendment 2** (before unlock; Tier 3–4 only; must implement the registered text verbatim; code posted with the diff) | E-dev with its realised-design calibration and the attrition battery; E-z and E-anchor; E-ADHD, E-int; X1–X7; the node-level Freedman–Lane max-T and Figure 5d; the GAMM display; the DWI multiverse (§24.5); node-wise R2\* along pSCP_R (computed after unlock by a frozen script) |
