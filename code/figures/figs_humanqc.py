"""Manuscript Figures 1-3 for the complete-human-QC sample (H_all180, n = 85). POST HOC.

Copied from manuscript/code/figs_pcn.py (2026-10-28; the original is not edited). The plotting code (journal style,
Arial 7-8.5 pt, widths 112 / 168 mm, Okabe-Ito palette, cell_panel(), the Figure 2 forest and Figure 3 layout, the
600-dpi LZW TIFF + 300-dpi JPEG writer) is unchanged. What changed:
  * data: instead of figs_pcn_data.R (n = 80 assert) and the tier-4 CSVs of the AI-QC sample, every number comes from
      ../supplement_H_all180.json   (code/supplement_humanqc.R; frozen helpers + tier4 helpers on the H_all180 table)
      ../frozen_run_H_all180/results.json   (frozen run_all.R, unchanged, on the H_all180 table)
    and is asserted against the frozen-run JSON where both exist;
  * Figure 2 rows: primary M1, Huber, CSF reference (K3), session-1 only, without threshold-psychosis youths, without
    prescribed-stimulant endorsers, + registration Dice, + motion index, joint nuisance model (Dice + motion + CSF
    offset + ICV), SN only, VTA + PBP only. All SEs are the fitted HC3 SEs (no back-calculation from p), and each row's
    95% CI uses that row's own residual df. Bound line = the frozen bound (-6.1318 ppb), adult line = -14.24 ppb;
  * Supplementary Figure S1 (TD centile chart) is not drawn.
Writes ONLY to ../figures/.   Run:  python3 figs_humanqc.py
"""
import io, json, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
OUT = os.path.join(BASE, "figures")
os.makedirs(OUT, exist_ok=True)

MM = 1 / 25.4
FULL, MID = 168 * MM, 112 * MM
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.0, "pdf.fonttype": 42, "svg.fonttype": "none",
    "mathtext.default": "regular",
})
CELLS = ["TD", "ADHD-only", "PS-only", "PS+ADHD"]
COL = {"TD": "#0072B2", "ADHD-only": "#E69F00", "PS-only": "#009E73", "PS+ADHD": "#CC79A7"}
MRK = {"TD": "o", "ADHD-only": "s", "PS-only": "^", "PS+ADHD": "D"}
VERMILLION, BLUE, GREY = "#D55E00", "#0072B2", "#7f7f7f"
CHECKS = []


def check(name, a, b, tol=1e-9):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = a.shape == b.shape and np.allclose(a, b, rtol=0, atol=tol)
    CHECKS.append((name, ok, float(np.max(np.abs(a - b))) if a.shape == b.shape else np.nan))
    assert ok, f"CHECK FAILED: {name}"


def save(fig, stem, tiff=True, png=False):
    """600-dpi LZW RGB TIFF + 300-dpi JPEG (or a 300-dpi PNG for the supplement)."""
    def render(dpi):
        buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=dpi, facecolor="white"); buf.seek(0)
        return Image.open(buf).convert("RGB")
    if tiff:
        render(600).save(os.path.join(OUT, stem + ".tif"), compression="tiff_lzw", dpi=(600, 600))
        render(300).save(os.path.join(OUT, stem + ".jpg"), quality=95, dpi=(300, 300), subsampling=0)
    if png:
        render(300).save(os.path.join(OUT, stem + ".png"), dpi=(300, 300))
    plt.close(fig)


# --------------------------------------------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------------------------------------------
SUP = json.load(open(os.path.join(BASE, "supplement_H_all180.json")))
FR = json.load(open(os.path.join(BASE, "frozen_run_H_all180", "results.json")))
assert SUP["tag"] == "H_all180"
first = pd.DataFrame(SUP["fig_points"])
first["cell"] = pd.Categorical(first.cell, CELLS)
N = len(first)
check("n first-usable = frozen run n", N, FR["H1"]["n"])
check("H_all180 n = 85 and beta_PS = 6.2809 (kit2_m1_refit.json)", [N, round(FR["H1"]["decision"]["beta_ppb"], 4)], [85, 6.2809])
T1 = SUP["table1"]
A_chi = pd.DataFrame(T1["chi_SNVTA_adjusted"]).set_index("cell").loc[CELLS]
A_r2 = pd.DataFrame(T1["r2s_SNVTA_adjusted"]).set_index("cell").loc[CELLS]
rawcsv = pd.concat([pd.DataFrame(T1["chi_SNVTA_raw"]), pd.DataFrame(T1["r2s_SNVTA_raw"])])


def verify_raw(v):
    for cl in CELLS:
        x = first.loc[first.cell == cl, v].dropna().to_numpy()
        r = rawcsv[(rawcsv.variable == v) & (rawcsv.cell == cl)].iloc[0]
        check(f"{v} {cl} raw summary vs supplement raw_by_cell",
              [len(x), x.mean(), x.std(ddof=1), np.median(x), np.quantile(x, .25), np.quantile(x, .75), x.min(), x.max()],
              [r.n, r["mean"], r.sd, r["median"], r.q1, r.q3, r["min"], r["max"]], tol=1e-7)


verify_raw("chi_SNVTA"); verify_raw("r2s_SNVTA")
check("cell counts vs frozen run", [int((first.cell == c).sum()) for c in CELLS], [FR["H1"]["cells"][c] for c in CELLS])


# --------------------------------------------------------------------------------------------------------------
# cell plot (Fig 1, Fig 3B): box (Tukey 1.5 IQR, quartiles type 7 as ggplot) + jittered points + adjusted mean/CI
# --------------------------------------------------------------------------------------------------------------
def cell_panel(ax, v, A, ylab, seed):
    rng = np.random.default_rng(seed)
    for i, cl in enumerate(CELLS):
        x = first.loc[first.cell == cl, v].dropna().to_numpy()
        ax.boxplot([x], positions=[i], widths=0.5, whis=1.5, showfliers=False, patch_artist=False,
                   boxprops=dict(color=COL[cl], lw=0.8), whiskerprops=dict(color=COL[cl], lw=0.8),
                   capprops=dict(color=COL[cl], lw=0), medianprops=dict(color=COL[cl], lw=1.6))
        ax.scatter(i + rng.uniform(-0.14, 0.14, len(x)), x, s=11, marker=MRK[cl], facecolor=COL[cl],
                   edgecolor="white", linewidth=0.3, alpha=0.9, zorder=3)
        a = A.loc[cl]
        ax.errorbar(i + 0.36, a.adj_mean, yerr=[[a.adj_mean - a.lo], [a.hi - a.adj_mean]], fmt="D", color="black",
                    ms=3.6, lw=1.0, capsize=0, zorder=4)
    ns = [int(first.loc[first.cell == c, v].notna().sum()) for c in CELLS]
    ax.set_xticks(range(4)); ax.set_xticklabels([f"{c}\n(n = {n})" for c, n in zip(CELLS, ns)])
    ax.set_xlim(-0.55, 3.65); ax.set_ylabel(ylab)
    ax.tick_params(axis="x", length=0)
    ax.yaxis.grid(True, color="#e6e6e6", lw=0.5); ax.set_axisbelow(True)
    return ns


# ---------------- Figure 1 ----------------
fig, ax = plt.subplots(figsize=(MID, 78 * MM))
ns1 = cell_panel(ax, "chi_SNVTA", A_chi, "SN–VTA χ (ppb)", seed=1)
check("Fig1 n per cell vs adjusted means", ns1, A_chi.n.to_numpy())
fig.tight_layout(pad=0.4)
save(fig, "Figure1")

# ---------------- Figure 2 (rows of this sample) ----------------
H = FR["H1"]; P = 12
co = {c["term"]: c for c in H["coefficients"]}
sens = {s["sensitivity"]: s for s in SUP["sensitivity"]["battery_frozen_function"]}
fsens = {s["sensitivity"]: s for s in FR["sensitivity"]["table"]}
K3 = SUP["H1_extras"]["K3_csf"]; s1 = SUP["sensitivity"]["ses1_only"][0]; jt = SUP["sensitivity"]["joint_nuisance"][0]
check("Fig2 K3 CSF beta / p vs frozen run", [K3["beta"], K3["p_two"]], [H["K3"]["csf_beta"], H["K3"]["csf_p"]], 1e-10)
check("Fig2 ses-1-only beta / p / n vs frozen run", [s1["beta_ppb"], s1["p_two"], s1["n"]],
      [H["ses1_only"]["beta"], H["ses1_only"]["p"], H["ses1_only"]["n"]], 1e-10)
rows = [  # label, beta, se, n, back-calculated?, df
    ("Primary M1 (whole-brain reference)", co["PS_e"]["est"], co["PS_e"]["se_hc3"], H["n"], False, H["n"] - P),
    ("Huber M-estimate", H["huber"]["est"], H["huber"]["se"], H["n"], False, H["n"] - P),
    ("CSF reference (K3)", K3["beta"], K3["se"], K3["n"], False, K3["df"]),
    ("Session-1 only", s1["beta_ppb"], s1["se_hc3"], s1["n"], False, s1["df"]),
]
for key, lab in [("without the 6 threshold-psychosis youths", "Without threshold-psychosis youths"),
                 ("without prescribed-stimulant endorsers (same session)", "Without stimulant endorsers"),
                 ("+ registration Dice (minimum of the G5 Dice values)", "+ registration Dice"),
                 ("+ magnitude ghost-to-signal motion index", "+ motion index"),
                 ("JOINT", "+ Dice, motion, CSF offset, ICV (joint)"),
                 ("SN-only complex (descriptive)", "SN only (SNc + SNr)"),
                 ("VTA + PBP only (descriptive)", "VTA + PBP only")]:
    s = jt if key == "JOINT" else sens[key]
    if key != "JOINT":
        check(f"Fig2 row '{lab}' vs frozen run", [s["beta_ppb"], s["se_hc3"], s["n"]],
              [fsens[key]["beta_ppb"], fsens[key]["se_hc3"], fsens[key]["n"]], 1e-10)
    rows.append((lab, s["beta_ppb"], s["se_hc3"], int(s["n"]), False, int(s["df"])))
bound = -6.13180464348368          # frozen registered bound (results_confirmatory: delta_B x sigma_B of n = 80)
adult = FR["H1"]["ppb"]["adult_reference"]
check("Fig2 bound / adult reference", [round(bound, 2), adult], [-6.13, -14.24])
fig, ax = plt.subplots(figsize=(FULL, 92 * MM))
y = np.arange(len(rows))[::-1]
for yi, (lab, b, se, n, bc, df) in zip(y, rows):
    tcrit = stats.t.ppf(0.975, df)
    col = "black" if yi == y[0] else ("#555555" if not bc else "#8a8a8a")
    ax.errorbar(b, yi, xerr=tcrit * se, fmt="o", color=col, ms=4.2 if yi == y[0] else 3.4, capsize=1.8, lw=1.0,
                mfc=col if not bc else "white", mew=0.9)
    ax.text(18.5, yi, f"{b:+.1f}  (n {n})", va="center", fontsize=7.5)
ax.axvline(0, color="grey", lw=0.7)
ax.axvline(bound, color=VERMILLION, ls="--", lw=1)
ax.axvline(adult, color=BLUE, ls=":", lw=1.3)
ax.text(bound + 0.4, len(rows) - 0.3, f"outcome-B bound\n{bound:.2f} ppb".replace("-", "−").replace("outcome−B", "outcome-B"), color=VERMILLION,
        fontsize=7.5, ha="left", va="bottom", bbox=dict(fc="white", ec="none", pad=0.6))
ax.text(adult - 0.4, len(rows) - 0.3, f"adult early-course\nschizophrenia\n{adult:.2f} ppb".replace("-1", "−1"), color=BLUE, fontsize=7.5, ha="right", va="bottom")
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows]); ax.tick_params(axis="y", length=0)
ax.set_xlim(-24, 25); ax.set_ylim(-0.8, len(rows) + 1.6)
ax.set_xlabel("PS effect on SN–VTA χ, ppb (β$_{PS}$, 95% CI)    ← lower in PS  |  higher in PS →")
ax.spines["left"].set_visible(False)
fig.tight_layout(pad=0.4)
save(fig, "Figure2")
FIG2_ROWS = rows

# ---------------- Figure 3 ----------------
ROIS = [("chi_SNVTA", "SN–VTA (primary)"), ("chi_SNonly", "SN only (SNc+SNr)"), ("chi_VTAPBP", "VTA + PBP"),
        ("chi_GP", "globus pallidus"), ("chi_PUT", "putamen"), ("chi_CAU", "caudate"), ("chi_STR", "striatum (CAU+PUT+NAcc)"),
        ("chi_THAL", "thalamus"), ("chi_RN", "red nucleus"), ("chi_STN", "subthalamic nucleus"), ("chi_DN_R", "right dentate"),
        ("chi_DN_L", "left dentate"), ("chi_DGM", "deep GM composite (CAU+GP)")]
fr = pd.DataFrame(SUP["fig3A_regional"]).set_index("outcome").loc[[r for r, _ in ROIS]]
assert len(fr) == 13
fr = fr.assign(lo=fr.ci_lo / fr.sd_raw, hi=fr.ci_hi / fr.sd_raw)
check("Fig3A point = beta_PS / sd_raw", fr.beta_per_sd_raw, fr.beta_PS / fr.sd_raw, 1e-12)
check("Fig3A SN-VTA row = frozen-run M1", fr.loc["chi_SNVTA", "beta_PS"], FR["H1"]["decision"]["beta_ppb"], 1e-10)
check("Fig3A SN-only / VTA+PBP rows = frozen-run sensitivity rows",
      [fr.loc["chi_SNonly", "beta_PS"], fr.loc["chi_VTAPBP", "beta_PS"]],
      [fsens["SN-only complex (descriptive)"]["beta_ppb"], fsens["VTA + PBP only (descriptive)"]["beta_ppb"]], 1e-10)
fig, (axA, axB) = plt.subplots(1, 2, figsize=(FULL, 82 * MM), gridspec_kw=dict(width_ratios=[1.0, 1.0]))
yy = np.arange(len(fr))[::-1]
axA.axvline(0, color="grey", lw=0.6)
axA.errorbar(fr.beta_per_sd_raw, yy, xerr=[fr.beta_per_sd_raw - fr.lo, fr.hi - fr.beta_per_sd_raw], fmt="o", color=BLUE,
             ms=3.4, lw=1.0, capsize=0)
xt = fr.hi.max() + 0.08
for yi, n in zip(yy, fr.n):
    axA.text(xt, yi, f"n = {n}", va="center", fontsize=6.8, color="#4d4d4d")
axA.set_yticks(yy); axA.set_yticklabels([lab for _, lab in ROIS]); axA.tick_params(axis="y", length=0)
axA.set_xlim(fr.lo.min() - 0.08, xt + 0.32); axA.set_ylim(-0.6, len(fr) - 0.4)
axA.set_xlabel("β$_{PS}$ in raw-SD units of the ROI\n(M1 form, HC3 95% CI)")
axA.spines["left"].set_visible(False)
axA.xaxis.grid(True, color="#ececec", lw=0.5); axA.set_axisbelow(True)
nsB = cell_panel(axB, "r2s_SNVTA", A_r2, "SN–VTA R2* (s$^{-1}$)", seed=2)
axB.tick_params(axis="x", labelsize=7)
check("Fig3B n per cell vs adjusted means", nsB, A_r2.n.to_numpy())
check("Fig3B R2* beta_PS = frozen-run K2", SUP["fig3B_r2s"]["beta_PS"], FR["H1"]["K2"]["est"], 1e-10)
fig.tight_layout(pad=0.4, w_pad=1.5, rect=(0, 0, 1, 0.95))
fig.text(0.005, 0.995, "A", fontsize=10, fontweight="bold", va="top", ha="left")
bbB = axB.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
fig.text(bbB.x0, 0.995, "B", fontsize=10, fontweight="bold", va="top", ha="left")
save(fig, "Figure3")

# --------------------------------------------------------------------------------------------------------------
print("sample H_all180, n =", N, "; merged sha256", SUP["merged_sha256"])
print("\nFigure 1 adjusted means [HC3 95% CI]:")
for c in CELLS:
    a = A_chi.loc[c]; print(f"  {c:10s} n={int(a.n):2d}  {a.adj_mean:.1f} [{a.lo:.1f}, {a.hi:.1f}]")
print("Figure 3B adjusted means R2*:")
for c in CELLS:
    a = A_r2.loc[c]; print(f"  {c:10s} n={int(a.n):2d}  {a.adj_mean:.2f} [{a.lo:.2f}, {a.hi:.2f}]")
print("Figure 2 rows (beta, 95% CI, n, df):")
for lab, b, se, n, bc, df in FIG2_ROWS:
    t = stats.t.ppf(0.975, df); print(f"  {lab:42s} {b:+6.2f} [{b - t * se:+6.2f}, {b + t * se:+6.2f}] n={n} df={df}")
print(f"  bound {bound:.4f}  adult {adult}")
print("Figure 3A rows (beta/SD [lo, hi]):")
for (v, lab) in ROIS:
    r = fr.loc[v]; print(f"  {lab:28s} {r.beta_per_sd_raw:+.3f} [{r.lo:+.3f}, {r.hi:+.3f}] n={r.n}")
print(f"\n{sum(ok for _, ok, _ in CHECKS)}/{len(CHECKS)} checks passed")
for nm, ok, d in CHECKS: print(f"  {'OK ' if ok else 'FAIL'} {nm} (max |diff| {d:.2e})")
print("\nwrote:", *sorted(os.listdir(OUT)), sep="\n  ")
