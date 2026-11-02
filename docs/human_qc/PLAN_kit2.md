# LEAD-χ human visual QC, kit 2: analysis plan (fixed before any kit-2 rating)

Written 2026-10-28 (about 19:20 UTC), before kit 2 was built and before any kit-2 rating exists. **POST HOC.** This plan changes no frozen file. It does not change the registered outcome (B; n = 80, β_PS = +4.9628 ppb).

---

## English

### 1. Purpose
Kit 1 (built 2026-10-07, rated by the first author on 2026-10-28) gave a human grade to 101 of the 180 AI-graded sessions. Kit 2 gives the same rater the remaining **79 AI-pass sessions** (the 119 AI grade-0/1 sessions minus the 40 in kit 1), so that **all 180 sessions have a human grade**. Kit 2 also holds **21 hidden repeats** of kit-1 items (11 AI-fail + 10 AI-pass, drawn at random from kit 1's 101 primary items, shown under new kit IDs). The repeats hide the composition of kit 2 and give a test-retest estimate. Total: **100 items**.

Same rubric, same montage pair (gross + placement), same image treatment, same `rate.html` and the same CSV export as kit 1. Sampling, kit IDs and order: `numpy.random.default_rng(20261028)` (recorded in `analysis/build_log_kit2.json`).

### 2. Complete-human-QC sample rule (fixed now)
- A session is **usable** if it meets every automated usability criterion of the primary analysis (the frozen `fail_reason` components `b_t1w`, `e_not_registered`, `no_qsm`, `c_auto_grade2`, `cprime_plausibility`, `e_g5`, kept unchanged) **and** its human grade is ≤ 1 **and** its human placement is `pass`. This is the kit-1 `H_replace` rule (the n = 89 run), now applied to all 180 sessions.
- Human grade and placement per session:
  - the **79 new sessions**: the kit-2 rating;
  - the **101 kit-1 sessions**: the kit-1 rating, combined exactly as in the n = 89 run (worse grade of the kit-1 original and its kit-1 repeat; placement fail if either fails). The kit-2 repeats are **not** used for the sample; they are used only for test-retest (section 4).
- A session that would become usable must also pass the frozen SN–VTA re-check of `apply_plausibility` (envelope [−50, 300] ppb, fimp ≤ 0.05, inmask ≥ 0.95), as in kit 1.
- **First usable session per youth** = the earliest usable `ses_order` (registered rule, `build_usable.py`).
- A kit-2 item left without a valid grade (or placement) keeps the AI grade (or placement), as in `H_replace`. The number of such sessions is reported; if it is not 0, the sample is not "complete human QC" and is reported as such.

### 3. Analysis (fixed now)
- `analysis/analyze_kit2.py` builds the per-session table and calls the **unchanged kit-1 R script** `human_qc_kit/analysis/refit_m1_human_qc.R` (SHA-256 `7c182ae2…e0e5`), which sources the frozen R helpers (`constants.R`, `utils.R`, `design.R`, `plausibility.R`, `m1.R`) read-only: `validate_table` → `code_vars` → `apply_plausibility` → `build_samples` → `run_m1_core`, as for the n = 89 run.
- Scenarios fitted in one call:
  1. `AI_frozen`: reproduction check; must give n = 80, β = 4.96280154409223 ppb (the run stops interpreting otherwise);
  2. `K1_H_replace`: kit-1 grades only (the n = 89 run). When the real kit-1 CSV (SHA-256 `1975118e…18c3`) is used, it must reproduce n = 89, β = 6.29013622956428 ppb;
  3. **`H_all180`: the complete-human-QC sample (primary result of this plan).**
- Reported for each: **n, β_PS, 95 % CI (HC3), two-sided p, and p_B** (one-sided bound test at the frozen registered bound −6.1318 ppb), plus the bound re-derived in the sample and the descriptive outcome label. Sample changes against the frozen sample (youths added or removed, changed first session) are listed.
- Secondary, descriptive: human-vs-AI agreement on the full census of 180 sessions (confusion matrix, linear-weighted κ, exact agreement; no design weights needed). Following the kit-1 rule, the primary census κ leaves out the example-montage session (`q15499e9a`, kit-1 item `hfaec73`) and the value with it is shown next to it. Also: the human fail rate among all 119 AI-pass sessions, with a Wilson 95 % CI.

### 4. Test-retest (fixed now)
- On the 21 repeats: kit-2 grade against the grade of the **kit-1 primary item** of the same session.
- Reported: **exact agreement** (with a Wilson 95 % CI) and **linear-weighted κ** (3 categories, same formula as kit 1; NA if either rating uses one category), plus the 3×3 table and placement agreement. Split by AI stratum (11 AI-fail, 10 AI-pass) as counts only.
- If the example-montage session is among the 21, the values without it are shown as well.
- Descriptive only; no threshold. The rater had seen the kit-1 images, so recall is possible and agreement may be optimistic.

### 5. Interpretation rule (fixed now)
- Post hoc. Reported **next to** the registered primary result (n = 80) and the kit-1 n = 89 run, never instead of them. The registered outcome stays B.
- A **nominally significant positive β_PS** (p < 0.05) is **not interpreted as a finding**: the registered hypothesis and its decision rule are the confirmatory test; this is a QC-sensitivity analysis.
- If `H_all180` keeps p_B < 0.05 and is not negative-significant (A), we write: "the result is robust to replacing the AI visual-QC grade with a first-author human grade for all 180 sessions." Otherwise: "the result is sensitive to the visual-QC rater."
- The rater (first author) is not independent (see `../rater_change_2026-10-28/RATER_CHANGE_NOTE.md`); this is stated wherever the result is reported.

### 6. Time evidence (before rating)
Before opening kit 2, the first author e-mails a co-author (`EMAIL_DRAFT_kit2.txt`): the SHA-256 of this plan, of the key `_key_DO_NOT_SEND/kit2_key.csv`, of the rater zip, and the `SHA256SUMS` list of every non-image kit-2 file with the list's own SHA-256. The finished CSV is e-mailed to the same co-author before `analyze_kit2.py` is run on it. The analysis code was tested only on simulated ratings (`analysis/test_run_SIMULATED/`).

---

## Türkçe

### 1. Amaç
Kit 1, AI'nin notladığı 180 oturumun 101'ine insan notu verdi. Kit 2, aynı değerlendiriciye kalan **79 AI-pass oturumu** (119 AI not-0/1 oturumdan kit 1'deki 40'ı çıkınca kalanlar) verir; böylece **180 oturumun hepsinin insan notu olur**. Kit 2'de ayrıca kit 1'in 101 birincil maddesinden rastgele çekilen **21 gizli tekrar** (11 AI-fail + 10 AI-pass, yeni kit ID'leriyle) var. Tekrarlar kitin bileşimini gizler ve test-tekrar test uyumunu verir. Toplam **100 madde**. Rubrik, montaj çifti, görüntü işlemi, `rate.html` ve CSV biçimi kit 1 ile aynı. Tohum: `numpy.random.default_rng(20261028)`.

### 2. Tam insan-QC örneklem kuralı
- Bir oturum, birincil analizin bütün otomatik kullanılabilirlik ölçütlerini karşılıyorsa (dondurulmuş `fail_reason` bileşenleri değişmeden) **ve** insan notu ≤ 1 **ve** insan placement'ı `pass` ise kullanılabilir. Bu, kit 1'deki `H_replace` kuralıdır (n = 89 koşusu), şimdi 180 oturumun hepsine uygulanıyor.
- 79 yeni oturumda kit-2 notu; kit-1'deki 101 oturumda kit-1 notu (n = 89 koşusundaki gibi: kit-1 orijinali ile kit-1 tekrarının kötüsü). Kit-2 tekrarları örneklem için kullanılmaz, yalnızca test-tekrar test içindir.
- Yeni kullanılabilir olan oturum dondurulmuş SN–VTA yeniden kontrolünden de geçmelidir.
- Her genç için ilk kullanılabilir oturum = en erken kullanılabilir `ses_order`.
- Notsuz kalan kit-2 maddesi AI notunu korur; sayısı raporlanır; 0 değilse örneklem "tam insan QC" sayılmaz.

### 3. Analiz
`analyze_kit2.py`, kit 1'in değişmemiş `refit_m1_human_qc.R` betiğini çağırır (dondurulmuş R yardımcıları salt okunur). Senaryolar: `AI_frozen` (n = 80, β = 4,9628 vermeli), `K1_H_replace` (gerçek kit-1 CSV'siyle n = 89, β = 6,2901 vermeli), **`H_all180` (bu planın birincil sonucu)**. Her biri için **n, β_PS, %95 GA, p ve p_B** (dondurulmuş sınır −6,1318 ppb) raporlanır. İkincil, betimsel: 180 oturumluk tam sayımda insan–AI uyumu (örnek montaj oturumu `hfaec73` birincil κ'dan çıkarılır, dahil hâli yanında) ve 119 AI-pass oturumda insan fail oranı (Wilson GA).

### 4. Test-tekrar test
21 tekrarda kit-2 notu, aynı oturumun kit-1 birincil maddesinin notuyla karşılaştırılır: **tam uyum** (Wilson GA) ve **doğrusal ağırlıklı κ**; 3×3 tablo ve placement uyumu. Yalnızca betimsel; hatırlama mümkün olduğundan iyimser olabilir.

### 5. Yorum kuralı
Post hoc; kayıtlı birincil sonucun (n = 80) ve n = 89 koşusunun **yanında** raporlanır. Kayıtlı sonuç B kalır. **Nominal olarak anlamlı pozitif bir β bulgu olarak yorumlanmaz.** `H_all180`'de p_B < 0,05 ve A yoksa: "sonuç, 180 oturumun hepsinde AI notunun ilk yazar insan notuyla değiştirilmesine karşı sağlam"; değilse "sonuç QC değerlendiricisine duyarlı". Değerlendiricinin bağımsız olmadığı her raporda belirtilir.

### 6. Zaman kanıtı
Kit 2 açılmadan önce bu planın, anahtarın ve zip'in SHA-256'sı ile `SHA256SUMS` listesi bir ortak yazara e-postayla gönderilir (`EMAIL_DRAFT_kit2.txt`). Bitmiş CSV, analizden önce aynı kişiye gönderilir. Analiz kodu yalnızca simüle notlarla test edildi.
