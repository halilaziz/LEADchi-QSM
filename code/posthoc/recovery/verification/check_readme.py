#!/usr/bin/env python3
"""Check each number printed in ../README.md against the independent recomputation, at the printed rounding."""
import json
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

HERE = Path(__file__).resolve().parent
me = json.loads((HERE / "recomputed.json").read_text())
fu, al = me["k_first_usable"], me["k_all"]
obs, cf, ts = me["observed"], me["counterfactual"], me["true_scale"]
sp, kg = me["k_spearman_L14"], me["k_by_group_L14"]
L = {"-5": "-5.0", "-10": "-10.0", "-14.24": "-14.24", "-20": "-20.0", "+10": "10.0"}


def rnd(x, printed):
    """Round x to the number of decimals shown in the printed string (sig. digits for 'e' notation)."""
    s = printed.replace("−", "-").replace("+", "")
    if "e" in s:
        mant, ex = s.split("e")
        nd = len(mant.split(".")[1]) if "." in mant else 0
        return f"{x:.{nd}e}".replace("e-0", "e-").replace("e+0", "e+")
    if s.startswith(".") or s.startswith("-."):
        nd = len(s.split(".")[1])
        r = str(Decimal(repr(x)).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))
        return r.replace("0.", ".", 1) if not r.startswith("-") else r.replace("-0.", "-.", 1)
    nd = len(s.split(".")[1]) if "." in s else 0
    return str(Decimal(repr(x)).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))


items = []


def chk(label, printed, x):
    s = printed.replace("−", "-").replace("+", "")
    r = rnd(x, printed)
    items.append((label, printed, x, r, "OK" if r == s else "MISMATCH"))


for lv in ["-5", "-10", "-14.24", "-20", "+10"]:
    v = fu[L[lv]]
    p = {"-5": ("0.470", "0.088", "0.474", "0.16", "0.96"), "-10": ("0.467", "0.059", "0.472", "0.30", "0.74"),
         "-14.24": ("0.473", "0.052", "0.473", "0.37", "0.68"), "-20": ("0.469", "0.047", "0.470", "0.37", "0.63"),
         "+10": ("0.477", "0.059", "0.480", "0.33", "0.70")}[lv]
    for nm, pr, x in zip(["mean", "SD", "median", "min", "max"], p, [v[1], v[2], v[3], v[6], v[7]]):
        chk(f"table k (FU) L={lv} {nm}", pr, x)
chk("n per level (all sessions)", "180", me["n_sessions_per_level"]["-14.24"])
chk("pooled slope", "0.474", me["pooled_linear"]["slope"])
chk("pooled intercept", "−0.05", me["pooled_linear"]["intercept"])
chk("pooled R2", "0.96", me["pooled_linear"]["r2"])
chk("pooled n", "900", me["pooled_linear"]["n"])
chk("k~group PS est", "+0.0006", kg["PS"]["est"])
chk("k~group PS CI lo", "−0.021", kg["PS"]["ci95"][0])
chk("k~group PS CI hi", "0.022", kg["PS"]["ci95"][1])
chk("k~group PS p", ".95", kg["PS"]["p_two"])
chk("k~group ADHD est", "+0.003", kg["ADHD"]["est"])
chk("k~group ADHD p", ".79", kg["ADHD"]["p_two"])
chk("cell means min", "0.46", min(kg["cell_means"]))
chk("cell means max", "0.49", max(kg["cell_means"]))
for v, pr in [("dice_min", ".02"), ("motion_gsr", "−.09"), ("roi_vol", "−.01"), ("Q_phase", ".04"), ("age", ".26")]:
    chk(f"Spearman rho {v}", pr, sp[v]["rho"])
chk("Spearman p age", ".02", sp["age"]["p"])
chk("ref shift (L=-14.24 mean)", "+0.006", me["d_ref_ppb"]["-14.24"][1])
chk("THAL shift L=-14.24", "−0.17", me["d_THAL_ppb"]["-14.24"][1])
chk("GP shift (L=-14.24 mean)", "+0.18", me["d_GP_ppb"]["-14.24"][1])
chk("GP SD (L=-14.24)", "1.6", me["d_GP_ppb"]["-14.24"][2])
chk("k multiplier", "0.47", ts["k_mean_L14"])
chk("beta_PS observed", "+4.96", obs["beta"])
chk("SE observed", "3.125", obs["se"])
chk("df observed", "68", obs["df"])
for row, (sh, pp) in zip(me["D_table"], [("−2.36", ".011"), ("−4.73", ".0014"), ("−6.73", ".0002"), ("−9.46", "9e-6")]):
    chk(f"D={row['D']} expected shift", sh, row["shift"])
    chk(f"D={row['D']} one-sided p", pp, row["p_full"])
chk("true-scale k", "0.473", ts["k_mean_L14"])
chk("B bound observed", "−6.13", obs["bound"])
chk("B bound true", "−13.0", ts["bound_true_ppb"])
chk("one-sided 95% lower observed", "−0.25", obs["ci90"][0])
chk("one-sided 95% lower true", "−0.53", ts["lower95_one_sided_true_ppb"])
chk("two-sided 95% lower observed", "−1.27", obs["ci95"][0])
chk("two-sided 95% lower true", "−2.7", ts["lower95_two_sided_true_ppb"])
for c, (b, o) in zip(cf, [("+2.61", "B"), ("+0.22", "B"), ("−1.80", "C"), ("−4.46", "C"), ("+9.78", "B")]):
    chk(f"counterfactual L={c['L_true_ppb']} beta", b, c["beta"])
    items.append((f"counterfactual L={c['L_true_ppb']} outcome", o, c["outcome"], c["outcome"], "OK" if o == c["outcome"] else "MISMATCH"))
chk("counterfactual +10 p_two", ".003", cf[4]["p_two"])
items.append(("counterfactual +10 disconfirming", "yes", cf[4]["disconfirming"], cf[4]["disconfirming"], "OK" if cf[4]["disconfirming"] else "MISMATCH"))
chk("MEDI k (caveat)", "0.47", fu["-14.24"][1])

out = ["label\tprinted\trecomputed\trounded\tstatus"] + ["\t".join(map(str, i)) for i in items]
(HERE / "readme_check.tsv").write_text("\n".join(out) + "\n")
print(sum(i[-1] == "OK" for i in items), "/", len(items), "OK")
for i in items:
    if i[-1] != "OK":
        print(i)
