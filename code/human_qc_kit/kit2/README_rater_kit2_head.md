# LEAD-χ visual QC re-rating, kit 2: instructions for the rater

*English below; Türkçe yönergeler 7. bölümde. The rubric (section 8) is the registered English text.*

Thank you for rating these images. Please read this page once before you start. It takes about 10 minutes.

## 1. Purpose

The quality of every brain scan in this study was graded by an **automated rater**, using a written rubric and two image montages per scan. We need an experienced human rater to grade a set of the same montages **independently**, with the same rubric, so that the automated grades can be checked. This is the second kit of the re-rating: same rubric, same montages, same tool as the first kit; only the items are different.

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

There are **100 items**. Each item is named by a random kit ID (`h` followed by 6 characters). The kit ID tells you nothing about the participant, the session or the group. The montage titles show only the kit ID.

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
- **Repeats.** Some items may be scans you have already rated, in this kit or in the first kit, shown under new kit IDs. Rate every item as if you were seeing it for the first time. Do not search for repeats, and do not look up your earlier ratings.
- **No expected proportions.** The items are a selection of the study's scans in random order. Do not assume any particular share of passes or fails.
- **Grade only what you see.** Some automated exclusion rules (for example, the automatic grade 2 in the rubric) are applied separately by software. You do not need to apply them, and you cannot see them.
- **Do not try to identify participants, sessions or groups.** Do not open any other file of this project (QC sheets, grade tables, data tables, manuscript drafts) before your ratings are returned. Do not discuss the items with the first author, or with anyone who has seen the automated grades, until then.

## 5. How to use `rate.html`

1. **Unzip the whole kit** and keep its folder structure: `rate.html` must stay next to the `images/` folder.
2. Open `rate.html` in a desktop browser (Chrome, Edge, Firefox or Safari), on a normal-sized screen.
   - The header should read **"Kit check: all 200 images found"**.
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

**Expected time:** about **1.5 to 2 hours** in total, possibly more if you open many items at full size. Breaks are fine, and there is no time limit: please do not rush.

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

## 7. Türkçe yönergeler (1–6. bölümlerin özeti)

**Amaç.** Bu çalışmadaki her MR oturumunun kalitesi, yazılı bir rubrik ve oturum başına iki montajla **otomatik bir değerlendirici** tarafından notlandı. Otomatik notları sınamak için aynı montajların bir kısmını, aynı rubrikle, **bağımsız** olarak notlamanı istiyoruz. Bu, yeniden notlamanın ikinci kitidir: rubrik, montajlar ve araç ilk kitle aynı; yalnızca maddeler farklı. Katılımcının grubunu, tanısını, yaşını, cinsiyetini ve hiçbir maddenin otomatik notunu bilmiyorsun; notlarını gönderene kadar böyle kalsın.

**Kitte ne var.** `rate.html` (çevrimdışı notlama aracı; ağ bağlantısı kurmaz), `images/gross/` ve `images/placement/` (her madde için iki montaj), `rating_sheet_blank.csv` (yedek kâğıt form), bu sayfa. **100 madde** var. Her madde rastgele bir kit ID ile adlandırılmış (`h` + 6 karakter); ID katılımcı, oturum ya da grup hakkında hiçbir şey söylemez.

**Montajlar.** Gross montajı 3 satır × 5 sütun: satırlar dentat, orta beyin (SN–VTA) ve bazal ganglion düzeyleri; sütunlar eko-1 ve eko-4 magnitüd, toplam alan, yerel alan ve χ (−0,2 ile +0,2 ppm, sarı Mask-3 sınırıyla). Placement montajı 2 satır × 4 sütun: üstte T1w, altta eko-1 magnitüd; SN–VTA'dan aksiyel, koronal, sagital kesit ve dentat düzeyinde aksiyel kesit; kırmızı çizgiler ROI sınırları. **Zorunlu inceleme bayrağı** turuncu kutu olarak görünür: o maddede her paneli tam boyutta incele; bayrak notu kendiliğinden değiştirmez.

**Notlama kuralları.** Her madde için: (1) gross notu 0 (geçer), 1 (hafif) veya 2 (başarısız), rubrik 2. bölüme göre; (2) uygun bütün hata kodları (M, W, S, B, K, F, A, C, D; not 2 en az bir kod ister; D'yi her gördüğünde işaretle); (3) placement `pass` veya `fail`, rubrik 3. bölüme göre ve **yalnızca** placement montajının T1w ve eko-1 satırlarına bakarak (χ'ya ya da gross montajına asla bakma; fail ise nedenini nota yaz); (4) isteğe bağlı kısa not. Her maddeyi kendi başına, rubriğe göre notla; maddeleri yan yana karşılaştırarak kalibre etme. **Tekrarlar:** bazı maddeler bu kitte ya da ilk kitte notladığın taramalar olabilir, yeni kit ID'leriyle; her maddeyi ilk kez görüyormuş gibi notla, tekrar arama ve eski notlarına bakma. **Beklenen oran yok:** geçen ya da kalan maddelerin belirli bir payını varsayma. Yalnızca gördüğünü notla. Katılımcı, oturum ya da grup tanımaya çalışma; notların dönene kadar projenin başka hiçbir dosyasını (QC tabloları, not tabloları, veri tabloları, makale taslakları, ilk kitin anahtarı ve analiz çıktıları) açma.

**`rate.html` kullanımı.** Kiti klasör yapısını bozmadan tümüyle aç; `rate.html` `images/` klasörünün yanında kalmalı. Masaüstü tarayıcıda aç; başlıkta **"Kit check: all 200 images found"** yazmalı. Başlığa baş harflerini yaz. Yakınlaştırma düğmeleri (Fit, 100–300 %), tam ekran (`G`, `L`; `+`/`−`; `Esc`). Kısayollar: `0` `1` `2` not; `M` `W` `S` `B` `K` `F` `A` `C` `D` kod; `P` placement; `←` `→` önceki/sonraki; `N` sonraki notsuz madde; `R` rubrik; `?` kısayollar. Alttaki şerit: yeşil = not ve placement girildi, sarı = eksik. İlerleme bu tarayıcıda otomatik saklanır; her molada **Export CSV** ile yedek al; başka tarayıcıda **Import CSV** ile geri yükle. Süre yaklaşık **1,5–2 saat**; mola verebilirsin, acele etme.

**Notları gönderme.** Bütün maddeler yeşil olunca **Export CSV**'ye bas (`LEADchi_humanQC_<baş harfler>_<tarih-saat>.csv`, sütunlar `kit_id, grade, failure_codes, placement, note, timestamp, rater`). Yalnızca bu CSV'yi kiti gönderen koordinatöre gönder. Araç çalışmazsa `rating_sheet_blank.csv`'yi doldurup **CSV UTF-8 (virgülle ayrılmış)** olarak gönder (`order` ve `kit_id` sütunlarını değiştirme; `grade` 0/1/2, `failure_codes` `;` ile ayrılmış harfler, `placement` `pass`/`fail`). Gönderdikten sonra notları değiştirme; hata fark edersen koordinatöre söyle, düzeltme olarak kaydedilir.

**Rubrik.** Aşağıdaki 8. bölüm, kayıtlı rubriğin İngilizce metnidir ve bilerek çevrilmeden, değiştirilmeden bırakıldı; notlamada geçerli olan bu metindir.

---

## 8. The rubric (registered text, verbatim)

The two sections below are the operative rubric for this re-rating. They are copied unchanged from the study's registered QC rubric, and the tool's "Rubric" button shows the same text. Section 1 of the registered rubric describes the blinding of the original grading round; it is reproduced in the Appendix for context only. For this re-rating, sections 1 to 7 of this page replace its procedural rules (for example, where grades are entered and who grades).

<!-- BEGIN VERBATIM: qc_rubric.md sections 2-3 -->
