# LEAD-χ: amendment 0 (route and governance; internal)

**Frozen:** 2026-09-29, **before any χ, local-field or R2\* map of the LEAD-χ pipeline existed**, and before any server step of the pipeline had run.

## Decisions (user, 2026-09-29)

1. **No OSF submission or public preregistration.**
   - The registered document is `LEADchi_prereg_v1.0_FINAL.md` (SHA-256 `139a99fbc92e63771b48dfd1f5a36196db03c0058207d64b4b70100017baa1ac`). It is frozen **internally** by that hash plus this amendment's hash and time stamp, as was done for Paper 0.
   - Wherever the registration says "posting", "OSF time stamp" or "registration DOI", read "internal freeze of this amendment".
   - Amendments 1 and 2 are frozen the same way: a local SHA-256 file with a UTC time stamp and a copy on the server under `LEADchi/prereg/`.
   - The manuscript will say "internally time-stamped analysis plan", not "preregistered".
2. **IRB.** Handled by the user outside this workflow. The data are public and de-identified (OpenNeuro ds007116, CC0), so the user considers an IRB determination unnecessary. The IRB gate in Q9 and Annex D is removed; the `[FILL]` fields for IRB are not blocking.
3. **Server approvals granted: A1, A3, A4, A5, A6** (all server operations of the pipeline), under the absolute rule that nothing is deleted, modified or created outside `/analysis/halilaziz/halil/Penn_LEAD` on any server. Optional items:
   - **A2** (SEPIA/MEDI) is **not** used. The Cornell MEDI toolbox needs the user's registration form, so G2 uses the registered within-QSMxT branch.
   - **A3-4 to A3-7** (AutoTrack, the D10 ΔR2\* inputs, CBF, XCP-D FC) are approved in principle. They are fetched only after the primary pipeline is frozen (amendment 2), because they are Tier 3–4.
4. **Blinding** stays procedural as registered (Annex A). The lock key pair is generated on the Mac by the user; the private key and passphrase never enter an AI session.
5. **`[FILL]` fields** that come from server steps (header-audit values, n_min, tool and atlas hashes, λ, G2, G5) are recorded in amendment 1 **before any χ used for a gate or a test is computed**. The pilot χ maps are the registered exception (Annex B), as written.

Nothing else in the registered document changes.
