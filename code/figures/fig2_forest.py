"""Figure 2: PS effect on SN-VTA chi (ppb) against the outcome-B bound and the adult reference.
Reads only results_confirmatory/results.json (frozen run). SEs for the CSF-referenced and ses-1-only refits are not
stored in the JSON; they are back-calculated from the two-sided p with the model df (n - 12) and marked in the legend."""
import json, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "../../results_confirmatory/results.json")))
H = R["H1"]; P = 12  # parameters in M1
co = {c["term"]: c for c in H["coefficients"]}
sens = {s["sensitivity"]: s for s in R["sensitivity"]["table"]}

def from_p(beta, p, n):
    t = stats.t.ppf(1 - p / 2, n - P)
    return abs(beta) / t

rows = [  # label, beta, se, n, back-calculated?
    ("Primary M1 (whole-brain reference)", co["PS_e"]["est"], co["PS_e"]["se_hc3"], H["n"], False),
    ("Huber M-estimate", H["huber"]["est"], H["huber"]["se"], H["n"], False),
    ("CSF reference (K3)", H["K3"]["csf_beta"], from_p(H["K3"]["csf_beta"], H["K3"]["csf_p"], H["n"]), H["n"], True),
    ("Session-1 only", H["ses1_only"]["beta"], from_p(H["ses1_only"]["beta"], H["ses1_only"]["p"], H["ses1_only"]["n"]), H["ses1_only"]["n"], True),
]
for key, lab in [("automated exclusions only (no visual grade or placement rating)", "Automated exclusions only"),
                 ("without the 6 threshold-psychosis youths", "Without threshold-psychosis youths"),
                 ("without prescribed-stimulant endorsers (same session)", "Without stimulant endorsers"),
                 ("+ registration Dice (minimum of the G5 Dice values)", "+ registration Dice"),
                 ("+ magnitude ghost-to-signal motion index", "+ motion index"),
                 ("+ puberty (PDS), complete case", "+ puberty (complete case)"),
                 ("SN-only complex (descriptive)", "SN only (SNc + SNr)"),
                 ("VTA + PBP only (descriptive)", "VTA + PBP only")]:
    s = sens[key]; rows.append((lab, s["beta_ppb"], s["se_hc3"], s["n"], False))

bound = H["ppb"]["bound_ppb"]
adult = H["ppb"]["adult_reference"]
fig, ax = plt.subplots(figsize=(7.2, 5.4))
y = np.arange(len(rows))[::-1]
for yi, (lab, b, se, n, bc) in zip(y, rows):
    tcrit = stats.t.ppf(0.975, n - P)
    col = "black" if yi == y[0] else ("#555555" if not bc else "#8a8a8a")
    ax.errorbar(b, yi, xerr=tcrit * se, fmt="o", color=col, ms=5 if yi == y[0] else 4, capsize=2, lw=1.2)
    ax.text(18.5, yi, f"{b:+.1f}  (n {n})", va="center", fontsize=8)
ax.axvline(0, color="grey", lw=0.8)
ax.axvline(bound, color="#c0392b", ls="--", lw=1)
ax.axvline(adult, color="#2c3e50", ls=":", lw=1.2)
ax.text(bound + 0.4, len(rows) - 0.3, f"outcome-B bound\n{bound:.2f} ppb", color="#c0392b", fontsize=8, ha="left", va="bottom")
ax.text(adult - 0.4, len(rows) - 0.3, f"adult early-course\nschizophrenia\n{adult:.2f} ppb", color="#2c3e50", fontsize=8, ha="right", va="bottom")
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
ax.set_xlim(-24, 25); ax.set_ylim(-0.8, len(rows) + 1.2)
ax.set_xlabel("PS effect on SN–VTA χ, ppb (β_PS, 95% CI)  ← lower in PS   |   higher in PS →", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
out = os.path.join(HERE, "../figures/Fig2_PS_effect_forest.png")
fig.savefig(out, dpi=300); print("wrote", out, "bound", bound, "adult", adult)
