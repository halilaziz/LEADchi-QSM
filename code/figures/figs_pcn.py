"""Submission figures for Psychiatry and Clinical Neurosciences (PCN).

Redraws manuscript Figures 1-3 and Supplementary Figure S1 from the frozen outputs, without in-figure titles,
in journal style (Arial, 7-8.5 pt at final size; widths 112 / 168 mm; Okabe-Ito colour-blind-safe palette).
Writes ONLY to ../figures_pcn/. Reads (never writes) the frozen trees:

  Figure 1  raw SN-VTA chi of the first-usable sessions (box + points) + M1-adjusted means, HC3 95% CI
            raw points  : first-usable sample rebuilt by figs_pcn_data.R (frozen helpers, read-only)
            adj. means  : results_exploratory/tier4_exploratory/adjusted_cell_means.csv (outcome chi_SNVTA)
            (same content as tier4_exploratory.R fig_D1 -> manuscript/figures/Fig1_SNVTA_chi_by_cell.png)
  Figure 2  forest of beta_PS: row logic executed verbatim from fig2_forest.py (results_confirmatory/results.json)
  Figure 3  (A) PS effect by ROI in raw-SD units: group_effects_by_outcome.csv, D2 rows (fig_D2_forest)
            (B) SN-VTA R2* by cell: raw points (figs_pcn_data.R) + adjusted means (adjusted_cell_means.csv, r2s_SNVTA)
            (fig_D4_r2s_SNVTA_by_cell)
  Fig. S1   TD-referenced centile chart (t3a_tier3_estimation.R ez_centile_chart): points from
            results_exploratory/tier3_estimation/ez_subject_z.csv; centile curves from the TD-only model refitted by
            figs_pcn_data.R and checked against the frozen mu_TD values.

Every plotted number is asserted against the source CSV/JSON before drawing.   Run:  python3 figs_pcn.py
"""
import io, json, os, subprocess
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "../.."))
OUT = os.path.normpath(os.path.join(HERE, "../figures_pcn"))
T4 = os.path.join(ROOT, "results_exploratory/tier4_exploratory")
T3 = os.path.join(ROOT, "results_exploratory/tier3_estimation")
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
# Okabe-Ito (colour-blind safe) + redundant marker shapes
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
rjs = subprocess.run(["Rscript", os.path.join(HERE, "figs_pcn_data.R")], capture_output=True, text=True, check=True, cwd=HERE)
RD = json.loads(rjs.stdout)
first = pd.DataFrame(RD["first"])
assert len(first) == 80
adj = pd.read_csv(os.path.join(T4, "adjusted_cell_means.csv"))
rawcsv = pd.read_csv(os.path.join(T4, "raw_by_cell.csv")).drop_duplicates(["variable", "cell"])
ge = pd.read_csv(os.path.join(T4, "group_effects_by_outcome.csv"))
ez = pd.read_csv(os.path.join(T3, "ez_subject_z.csv"))
anc = pd.read_csv(os.path.join(T3, "eanchor_z.csv"))


def verify_raw(v):
    """raw points vs the frozen raw_by_cell.csv summary (n, mean, sd, median, quartiles type 7, min, max)."""
    for cl in CELLS:
        x = first.loc[first.cell == cl, v].dropna().to_numpy()
        r = rawcsv[(rawcsv.variable == v) & (rawcsv.cell == cl)].iloc[0]
        check(f"{v} {cl} raw summary vs raw_by_cell.csv",
              [len(x), x.mean(), x.std(ddof=1), np.median(x), np.quantile(x, .25), np.quantile(x, .75), x.min(), x.max()],
              [r.n, r["mean"], r.sd, r["median"], r.q1, r.q3, r["min"], r["max"]], tol=1e-7)


verify_raw("chi_SNVTA"); verify_raw("r2s_SNVTA")
# chi values of the R-rebuilt sample equal the frozen E-z subject table (independent frozen output, same 80 sessions)
m = first.merge(ez[["sub", "session_id", "cell", "age", "chi_SNVTA"]], on=["sub", "session_id"], suffixes=("", "_ez"))
assert len(m) == 80 and (m.cell == m.cell_ez).all()
check("chi_SNVTA first-usable vs ez_subject_z.csv", m.chi_SNVTA, m.chi_SNVTA_ez, 1e-7)
check("age first-usable vs ez_subject_z.csv", m.age, m.age_ez, 1e-9)

# Table 1 adjusted means (rounded to 0.1) must be reproduced by the plotted values
A_chi = adj[adj.outcome == "chi_SNVTA"].set_index("cell").loc[CELLS]
A_r2 = adj[adj.outcome == "r2s_SNVTA"].set_index("cell").loc[CELLS]
TABLE1 = {"TD": (49.7, 44.1, 55.3), "ADHD-only": (51.6, 43.5, 59.8), "PS-only": (57.9, 52.3, 63.5), "PS+ADHD": (53.4, 48.9, 57.9)}
check("Fig1 adjusted means/CI rounded vs Table 1", np.round(A_chi[["adj_mean", "lo", "hi"]].to_numpy(), 1),
      np.array([TABLE1[c] for c in CELLS]), 1e-9)


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
check("Fig1 n per cell vs adjusted_cell_means.csv", ns1, A_chi.n.to_numpy())
fig.tight_layout(pad=0.4)
save(fig, "Figure1")

# ---------------- Figure 2 (row logic of fig2_forest.py executed verbatim, up to the plotting block) ----------------
src = open(os.path.join(HERE, "fig2_forest.py")).read()
ns2 = {"__file__": os.path.join(HERE, "fig2_forest.py")}
exec(compile(src.split("\nfig, ax = plt.subplots")[0], "fig2_forest.py[data part]", "exec"), ns2)
rows, bound, adult, P = ns2["rows"], ns2["bound"], ns2["adult"], ns2["P"]
check("Fig2 bound / adult reference", [round(bound, 2), adult], [-6.13, -14.24])
fig, ax = plt.subplots(figsize=(FULL, 92 * MM))
y = np.arange(len(rows))[::-1]
for yi, (lab, b, se, n, bc) in zip(y, rows):
    tcrit = stats.t.ppf(0.975, n - P)
    col = "black" if yi == y[0] else ("#555555" if not bc else "#8a8a8a")
    ax.errorbar(b, yi, xerr=tcrit * se, fmt="o", color=col, ms=4.2 if yi == y[0] else 3.4, capsize=1.8, lw=1.0,
                mfc=col if not bc else "white", mew=0.9)
    ax.text(18.5, yi, f"{b:+.1f}  (n {n})", va="center", fontsize=7.5)
ax.axvline(0, color="grey", lw=0.7)
ax.axvline(bound, color=VERMILLION, ls="--", lw=1)
ax.axvline(adult, color=BLUE, ls=":", lw=1.3)
ax.text(bound + 0.4, len(rows) - 0.3, f"outcome-B bound\n{bound:.2f} ppb".replace("-", "\u2212").replace("outcome\u2212B", "outcome-B"), color=VERMILLION,
        fontsize=7.5, ha="left", va="bottom", bbox=dict(fc="white", ec="none", pad=0.6))
ax.text(adult - 0.4, len(rows) - 0.3, f"adult early-course\nschizophrenia\n{adult:.2f} ppb".replace("-1", "\u22121"), color=BLUE, fontsize=7.5, ha="right", va="bottom")
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows]); ax.tick_params(axis="y", length=0)
ax.set_xlim(-24, 25); ax.set_ylim(-0.8, len(rows) + 1.6)
ax.set_xlabel("PS effect on SN–VTA χ, ppb (β$_{PS}$, 95% CI)    ← lower in PS  |  higher in PS →")
ax.spines["left"].set_visible(False)
fig.tight_layout(pad=0.4)
save(fig, "Figure2")
FIG2_ROWS = [(lab, b, se, n, bc) for lab, b, se, n, bc in rows]

# ---------------- Figure 3 ----------------
ROIS = [("chi_SNVTA", "SN–VTA (primary)"), ("chi_SNonly", "SN only (SNc+SNr)"), ("chi_VTAPBP", "VTA + PBP"),
        ("chi_GP", "globus pallidus"), ("chi_PUT", "putamen"), ("chi_CAU", "caudate"), ("chi_STR", "striatum (CAU+PUT+NAcc)"),
        ("chi_THAL", "thalamus"), ("chi_RN", "red nucleus"), ("chi_STN", "subthalamic nucleus"), ("chi_DN_R", "right dentate"),
        ("chi_DN_L", "left dentate"), ("chi_DGM", "deep GM composite (CAU+GP)")]
fr = ge[(ge.module == "D2") & ge.test.str.startswith("regional PS effect: ")].set_index("outcome").loc[[r for r, _ in ROIS]]
assert len(fr) == 13
fr = fr.assign(lo=fr.ci_lo / fr.sd_raw, hi=fr.ci_hi / fr.sd_raw)
check("Fig3A point = beta_PS / sd_raw", fr.beta_per_sd_raw, fr.beta_PS / fr.sd_raw, 1e-12)
check("Fig3A SN-VTA row = frozen M1 (4.963 ppb)", round(fr.loc["chi_SNVTA", "beta_PS"], 3), 4.963)
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
check("Fig3B n per cell vs adjusted_cell_means.csv", nsB, A_r2.n.to_numpy())
fig.tight_layout(pad=0.4, w_pad=1.5, rect=(0, 0, 1, 0.95))
# panel letters at the top-left of each panel's full extent (panel A's ROI labels extend left of its axes)
fig.text(0.005, 0.995, "A", fontsize=10, fontweight="bold", va="top", ha="left")
bbB = axB.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
fig.text(bbB.x0, 0.995, "B", fontsize=10, fontweight="bold", va="top", ha="left")
save(fig, "Figure3")

# ---------------- Supplementary Figure S1 ----------------
E = RD["ez"]
ag, mu, sig, Bk = np.array(E["age_grid"]), np.array(E["mu_grid"]), E["sigma_td"], E["boundary"]
chk = pd.DataFrame(E["check_nonTD"]).merge(ez[ez.cell != "TD"][["sub", "mu_TD", "sigma_TD"]], on="sub", suffixes=("", "_frozen"))
assert len(chk) == 58
check("FigS1 TD model mu vs frozen ez_subject_z mu_TD (non-TD rows)", chk.mu_TD, chk.mu_TD_frozen, 1e-8)
check("FigS1 sigma_TD vs frozen", [sig], [chk.sigma_TD.iloc[0]], 1e-10)
zq = [-1.96, -1, 0, 1, 1.96]
anc_us = anc[anc.usable == 1]
yl = (min(ez.chi_SNVTA.min(), (mu - 1.96 * sig).min(), anc_us.chi_SNVTA.min()),
      max(ez.chi_SNVTA.max(), (mu + 1.96 * sig).max(), anc_us.chi_SNVTA.max()))
yl = (yl[0] - 0.04 * (yl[1] - yl[0]), yl[1] + 0.04 * (yl[1] - yl[0]))
xl = (min(ez.age.min(), ag.min()), max(ez.age.max(), ag.max()))
fig, axs = plt.subplots(2, 2, figsize=(FULL, 128 * MM), sharex=True, sharey=True)
for ax, cl in zip(axs.flat, CELLS):
    ax.fill_between(ag, mu - 1.96 * sig, mu + 1.96 * sig, color=BLUE, alpha=0.07, lw=0)
    ax.fill_between(ag, mu - sig, mu + sig, color=BLUE, alpha=0.10, lw=0)
    for q in zq:
        ax.plot(ag, mu + q * sig, color=BLUE, lw=1.4 if q == 0 else 0.8, ls="-" if q == 0 else (0, (4, 2)), alpha=0.9 if q == 0 else 0.6)
    for q, t in zip(zq, ["2.5", "16", "50", "84", "97.5"]):
        ax.text(Bk[1] + 0.12, mu[-1] + q * sig, t, fontsize=6.3, color="#4d4d4d", va="center", ha="left")
    dd = ez[ez.cell == cl]
    ax.scatter(dd.age, dd.chi_SNVTA, s=14, marker=MRK[cl], facecolor=COL[cl], edgecolor="white", linewidth=0.3, zorder=3)
    if cl == "PS+ADHD":
        aa = dd[dd.threshold_psychosis]
        ax.scatter(aa.age, aa.chi_SNVTA, s=34, marker="^", facecolor="black", edgecolor="black", zorder=4)
        ax.legend(handles=[Line2D([], [], ls="", marker=MRK[cl], mfc=COL[cl], mec="white", ms=4.5, label="PS+ADHD"),
                           Line2D([], [], ls="", marker="^", mfc="black", mec="black", ms=5.5, label="threshold-psychosis youth (first usable)")],
                  loc="lower right", frameon=False, handletextpad=0.3, borderaxespad=0.2)
    if cl == "TD":
        ax.legend(handles=[Line2D([], [], ls="", marker=MRK[cl], mfc=COL[cl], mec="white", ms=4.5, label="TD (scored leave-one-out)"),
                           Line2D([], [], color=BLUE, lw=0.8, ls=(0, (4, 2)), label="TD centiles (sex-averaged)")],
                  loc="lower right", frameon=False, handletextpad=0.4, borderaxespad=0.2)
    ax.text(0.02, 0.97, f"{cl} (n = {len(dd)})", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top", ha="left")
    ax.set_xlim(xl[0] - 0.2, xl[1] + 1.0); ax.set_ylim(*yl)
    ax.grid(True, color="#ececec", lw=0.5); ax.set_axisbelow(True)
    check(f"FigS1 n {cl}", len(dd), {"TD": 22, "ADHD-only": 16, "PS-only": 17, "PS+ADHD": 25}[cl])
for ax in axs[1]: ax.set_xlabel("Age at first usable session (years)")
for ax in axs[:, 0]: ax.set_ylabel("SN–VTA χ (ppb)")
fig.tight_layout(pad=0.4, h_pad=0.8, w_pad=0.8)
save(fig, "FigureS1", tiff=False, png=True)

# --------------------------------------------------------------------------------------------------------------
print("input sha256", RD["input_sha256"])
print("\nFigure 1 adjusted means [HC3 95% CI] (adjusted_cell_means.csv):")
for c in CELLS:
    a = A_chi.loc[c]; print(f"  {c:10s} n={a.n:2d}  {a.adj_mean:.1f} [{a.lo:.1f}, {a.hi:.1f}]")
print("Figure 3B adjusted means R2*:")
for c in CELLS:
    a = A_r2.loc[c]; print(f"  {c:10s} n={a.n:2d}  {a.adj_mean:.2f} [{a.lo:.2f}, {a.hi:.2f}]")
print("Figure 2 rows (beta, 95% CI, n, back-calculated):")
for lab, b, se, n, bc in FIG2_ROWS:
    t = stats.t.ppf(0.975, n - P); print(f"  {lab:40s} {b:+6.2f} [{b - t * se:+6.2f}, {b + t * se:+6.2f}] n={n} {'(bc)' if bc else ''}")
print(f"  bound {bound:.4f}  adult {adult}")
print("Figure 3A rows (beta/SD [lo, hi]):")
for (v, lab) in ROIS:
    r = fr.loc[v]; print(f"  {lab:28s} {r.beta_per_sd_raw:+.3f} [{r.lo:+.3f}, {r.hi:+.3f}] n={r.n}")
print(f"\n{sum(ok for _, ok, _ in CHECKS)}/{len(CHECKS)} checks passed")
for nm, ok, d in CHECKS: print(f"  {'OK ' if ok else 'FAIL'} {nm} (max |diff| {d:.2e})")
print("\nwrote:", *sorted(os.listdir(OUT)), sep="\n  ")
