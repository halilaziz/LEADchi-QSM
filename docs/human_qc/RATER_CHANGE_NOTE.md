# LEAD-χ human visual-QC re-rating: change of rater (written 2026-10-28, before any rating)

**Status:** post hoc (amendment1_log #62). The kit, its key, its analysis code and the pre-specified interpretation rules (README_coordinator.md sections 6–7) are unchanged; their SHA-256 list is `human_qc_kit/SHA256SUMS` (sha256 `a6decd994644d45e44d5e7734fc4ea38989c2922e4cf9a4fd9ac4f9cd0003a76`).

## What changes
- **Rater:** the first author (Halil Aziz Velioglu) will rate the kit himself. The kit's coordinator note listed the first author as not suitable, because he is not independent. No independent rater was available, so the re-rating is done by the first author and will be reported as a **first-author re-rating**, not as an independent one.
- **What the rater has seen before rating:**
  - three example QC montages in `LEADchi/qc_ornekler/` (shown in the 2026-10-04 working session), one of which is in the kit;
  - the aggregate QC results (about one third of sessions graded 2 by the AI rater) and the confirmatory outcome;
  - the kit design in outline (all AI grade-2 sessions, a random sample of AI grade-0/1 sessions, and hidden repeats);
  - the rubric, which he co-designed for the registration.
- **What the rater has not seen and commits not to open until his ratings are e-mailed:** the per-item AI grades (`human_qc_kit/_key_DO_NOT_SEND/`, `mirror/qc/`, `prereg/qc_ai/`, `qc_montages_archive/` file names), the kit coordinator note, the simulated test runs (`human_qc_kit/analysis/test_run*`), and anything that maps kit IDs to codes, sessions or groups. The kit IDs are random, so the montages are blind to group and to the AI grade.

## Pre-specified handling (fixed now, before any rating)
1. The pre-specified example-montage rule (README_coordinator.md §6.2, review item L6) applies: the analysis is run with the example montage that is in the kit excluded from Part A, and this run is primary; the run without the exclusion is reported alongside.
2. All other rules and thresholds are unchanged (linear-weighted κ ≥ 0.60 counts as acceptable agreement; Part B scenarios as specified).
3. Reporting: "The first author, blind to group and to the item-level automated grades (randomised kit identifiers), re-graded 101 sessions (61 automated fails, 40 random automated passes, 10 hidden repeats); he was not independent of the study (he had seen example montages and aggregate QC results)."
4. Time evidence: this note's SHA-256 is e-mailed to a co-author together with the kit's SHA-256 list **before** rating starts; the finished CSV is e-mailed to the same co-author **before** any analysis is run.

---

**Türkçe özet:** Değerlendiriciyi ilk yazar olarak değiştiriyoruz. Bu post hoc bir sapmadır. Kit, anahtar, analiz kodu ve yorum kuralları değişmiyor. Değerlendirme makalede "bağımsız olmayan, ilk yazar tarafından yapılmış yeniden notlama" olarak raporlanacak. Örnek montaj kuralı birincil analizde uygulanacak. Not vermeye başlamadan önce bu notun ve kit listesinin hash'leri bir ortak yazara e-postayla gönderilecek; bitmiş CSV de analizden önce aynı kişiye gönderilecek.
