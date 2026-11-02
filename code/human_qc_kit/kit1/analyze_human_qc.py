#!/usr/bin/env python3
"""analyze_human_qc.py - LEAD-chi: validate the AI visual-QC rater (amendment 1c) against a blinded human rater.

POST HOC / EXPLORATORY. Nothing here changes a frozen file or re-decides H1; every output is labelled post hoc.

Usage
  python3 analyze_human_qc.py --ratings <returned CSV> --out <new output folder> [--simulated]
         [--key ../_key_DO_NOT_SEND/kit_key.csv] [--nboot 2000] [--seed 20261008] [--rscript Rscript]
         [--nextrap 1000] [--extrap-seed 20261009] [--exclude-from-a hfaec73]

The returned CSV is the export of rate.html (kit_id, grade, failure_codes, placement, note, timestamp, rater) or the
filled rating_sheet_blank.csv (extra columns are ignored). Items without a valid grade are reported and left out.

Part A - agreement (kit sample: every AI grade-2 session + a simple random sample of 40 AI grade-0/1 sessions)
  * 3x3 confusion matrix human x AI; linear-weighted kappa (the registered G4 statistic, same formula as
    qc_montage.py wkappa) and quadratic-weighted kappa; exact agreement.
  * Two estimands: (i) the kit sample as rated (unweighted); (ii) design-weighted to the 180 graded sessions
    (AI grade-2 stratum is a census, weight 1; AI grade-0/1 stratum weight 119/40). 95% CIs by stratified bootstrap.
  * Within strata: AI-fail stratum (AI grade constant, kappa undefined) -> human confirmation rate P(human = 2),
    Wilson CI; AI-pass stratum -> linear-weighted kappa, P(human = 2) = AI "missed fail" rate, Wilson CI, and the
    implied number of the 119 AI-pass sessions a human would fail.
  * Binary fail (grade 2) vs pass (0/1): Cohen kappa, PPV, NPV, sensitivity and specificity of the AI with the human
    as reference (design-weighted).
  * Agreement by forced-review flag (the flag is shown to both raters, as in the registered rubric).
  * Failure codes: per-code presence agreement (Cohen kappa), placement pass/fail agreement.
  * Intra-rater reliability on the 10 hidden duplicates: linear-weighted kappa, exact and placement agreement.
Part B - POST HOC sensitivity of H1 to the visual-QC rater
  * Per session, the human rating is the conservative combination of the original and its duplicate (worse grade;
    placement fail if either fails), the same rule qc_montage.py finalize applies to re-rated sessions.
  * The registered first-usable rule (build_usable.py) is re-applied: a session is usable iff it has none of
    b_t1w, e_not_registered, no_qsm, c_qc_grade (grade not 0/1), c_auto_grade2, cprime_plausibility,
    d_placement (placement not pass), e_g5; first usable = earliest usable ses_order per subject. All components
    except c_qc_grade and d_placement are taken unchanged from the frozen `fail_reason` column of the unlocked merged
    table; c_qc_grade and d_placement are recomputed from the scenario's grade/placement (the AI placement comes from
    grades_final.csv, which the merged table does not carry). The reconstruction is first checked to reproduce the
    frozen `usable`, `first_usable` and `fail_reason` exactly (the run stops otherwise).
  * Scenarios: AI_frozen (reproduction check, must give beta_PS = +4.9628 ppb, n = 80), H_replace (human grade and
    placement replace the AI for rated sessions; primary post hoc scenario), H_grade_only (placement kept from the
    AI), Strict_union (fail if either rater fails), Lenient (pass if either rater passes). Unrated sessions keep the
    AI rating.
  * A newly usable session must also pass the frozen SN-VTA re-verification of apply_plausibility (envelope
    [-50, 300] ppb, fimp_SNVTA <= 0.05, inmask_SNVTA >= 0.95); otherwise it is kept not usable and listed.
  * M1 is refitted by refit_m1_human_qc.R, which sources the frozen R helpers (identical design matrix), and the
    bound test is reported at the frozen registered bound (-6.1318 ppb) and at a bound re-derived in the sample.

Changes of 2026-10-07 after the adversarial review (REVIEW_by_checker.md; all made BEFORE any human rating existed;
see ../CHANGES_2026-10-07.md):
  * M1  wkappa returns NA (not 0) when either rater uses a single category; the forced-review flag = 1 subset is
        reported as a confirmation rate P(human = 2) with a Wilson CI, never as a kappa.
  * M2  the linear-weighted kappa on flag = 0 items (sample and design-weighted, bootstrap CI) is co-reported with
        the same prominence as the primary kappa.
  * M3  coverage of the H1 sensitivity is reported, and an extrapolation imputes human fails among the AI-pass
        sessions that are not in the kit (--nextrap draws; rate ~ Jeffreys Beta(m + 1/2, n - m + 1/2) from the
        observed session-level missed-fail rate), refitting M1 for every draw.
  * M4  the ratings file may be comma-, semicolon- or tab-separated (decided from the header line), UTF-8 with or
        without BOM, UTF-16, or Windows-1254; grades such as "1.0" or "1,0" are accepted; the run stops with a clear
        message and writes nothing when no valid grade is left.
  * L4  identical code sets are counted only over duplicate pairs in which both items have a grade.
  * L5  secondary comparator: agreement with the AI pass-1 grade (grade column only; AI notes are never read) and the
        AI's own pass-1 vs pass-2 repeatability on its re-rated sessions.
  * L6  --exclude-from-a <kit_id,...> drops the listed items (and any repeat of the same session) from Part A only.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import io
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

LX = "/path/to/Penn_LEAD/LEADchi"
HERE = os.path.dirname(os.path.abspath(__file__))
DEF = {
    "key": os.path.join(os.path.dirname(HERE), "_key_DO_NOT_SEND", "kit_key.csv"),
    "merged": os.path.join(LX, "mirror/unlocked/merged_qsmxt_primary_SNVTA.csv"),
    "grades": os.path.join(LX, "mirror/qc/grades_final.csv"),
    "gates": os.path.join(LX, "mirror/gates/gates_preunlock.json"),
    "rdir": os.path.join(LX, "mirror/scripts/analysis"),
    "confirmatory": os.path.join(LX, "results_confirmatory/results.json"),
    "pass1": os.path.join(LX, "mirror/qc/rating_pass1.csv"),       # only columns qc_id and grade are read (L5)
    "pass2": os.path.join(LX, "mirror/qc/rating_pass2.csv"),       # only columns qc_id and grade are read (L5)
}
CODES = list("MWSBKFACD")
REASONS = ["b_t1w", "e_not_registered", "no_qsm", "c_qc_grade", "c_auto_grade2", "cprime_plausibility",
           "d_placement", "e_g5"]
FROZEN_BETA = 4.96280154409223
FROZEN_BOUND = -6.13180464348368
N_AI_PASS_POOL = 119


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------ statistics
def wkappa(x, y, k=3, weights="linear", w=None):
    """Weighted Cohen kappa (same formula as qc_montage.py wkappa), optional case weights. Returns (kappa, O).

    Returns NaN (reported as NA) when either rater uses a single category: kappa is then undefined (P = E gives
    exactly 0, which would wrongly read as 'chance agreement'). Review item M1, 2026-10-07."""
    x = np.asarray(x, int)
    y = np.asarray(y, int)
    w = np.ones(len(x)) if w is None else np.asarray(w, float)
    O = np.zeros((k, k))
    np.add.at(O, (x, y), w)
    n = O.sum()
    if n <= 0:
        return float("nan"), O
    if (O.sum(1) > 0).sum() < 2 or (O.sum(0) > 0).sum() < 2:
        return float("nan"), O
    P = O / n
    i, j = np.indices((k, k))
    W = (np.abs(i - j) / (k - 1) if weights == "linear" else ((i - j) / (k - 1)) ** 2 if weights == "quadratic"
         else (i != j).astype(float))
    E = np.outer(P.sum(1), P.sum(0))
    den = (W * E).sum()
    return (float(1 - (W * P).sum() / den) if den > 1e-12 else float("nan")), O


def wilson(k, n, z=1.959964):
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [float(c - h), float(c + h)]


def binmetrics(h2, a2, w):
    """AI vs human (reference) for fail = grade 2. Design-weighted."""
    h2, a2, w = np.asarray(h2, bool), np.asarray(a2, bool), np.asarray(w, float)
    tp, fp = (w * (a2 & h2)).sum(), (w * (a2 & ~h2)).sum()
    fn, tn = (w * (~a2 & h2)).sum(), (w * (~a2 & ~h2)).sum()
    f = lambda a, b: float(a / (a + b)) if (a + b) > 0 else float("nan")
    return {"PPV_P(human2|AI2)": f(tp, fp), "NPV_P(human<2|AI<2)": f(tn, fn), "sensitivity_P(AI2|human2)": f(tp, fn),
            "specificity_P(AI<2|human<2)": f(tn, fp), "est_human_fail_n_of_180": float(tp + fn)}


def boot_ci(df, fn, nboot, rng):
    """Stratified (by AI stratum) percentile bootstrap of a scalar statistic fn(df)."""
    groups = [g.index.to_numpy() for _, g in df.groupby("stratum")]
    vals = []
    for _ in range(nboot):
        idx = np.concatenate([rng.choice(g, size=len(g), replace=True) for g in groups])
        v = fn(df.loc[idx])
        if v == v:
            vals.append(v)
    if len(vals) < 0.5 * nboot:
        return {"ci95": [float("nan"), float("nan")], "n_valid": len(vals)}
    return {"ci95": [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))], "n_valid": len(vals)}


def codeset(s):
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return set()
    return {c.strip().upper() for c in str(s).replace(",", ";").replace(" ", ";").split(";") if c.strip()}


# ------------------------------------------------------------------ inputs
def read_ratings(path, key):
    """Read the returned ratings robustly (review item M4): delimiter from the header line, BOM / UTF-16 / cp1254."""
    raw = open(path, "rb").read()
    enc = "utf-16" if raw[:2] in (b"\xff\xfe", b"\xfe\xff") else "utf-8-sig"
    try:
        text = raw.decode(enc)
    except UnicodeDecodeError:
        enc = "cp1254"          # Turkish Windows code page (spreadsheet 'CSV' that is not UTF-8); affects notes only
        text = raw.decode(enc, errors="replace")
    lines = text.lstrip("\ufeff").splitlines()
    if not lines or not lines[0].strip():
        sys.exit("ratings file is empty; nothing was written")
    head = lines[0]
    sep = "," if "," in head else ";" if ";" in head else "\t" if "\t" in head else ","
    r = pd.read_csv(io.StringIO(text.lstrip("\ufeff")), sep=sep, dtype=str, keep_default_na=False)
    r.columns = [str(c).strip().lstrip("\ufeff").strip().lower() for c in r.columns]
    for c in ("kit_id", "grade", "placement"):
        if c not in r.columns:
            sys.exit(f"ratings file lacks column {c} (delimiter read as {sep!r}; columns {list(r.columns)}); "
                     "nothing was written")
    for c in ("failure_codes", "note", "timestamp"):
        if c not in r.columns:
            r[c] = ""
    r["kit_id"] = r["kit_id"].str.strip().str.lower()
    issues = {"file_encoding": enc, "file_delimiter": {",": "comma", ";": "semicolon", "\t": "tab"}[sep],
              "unknown_kit_ids": sorted(set(r.kit_id) - set(key.kit_id) - {""}),
              "duplicated_rows": sorted(r.kit_id[r.kit_id.duplicated() & (r.kit_id != "")].unique().tolist())}
    r = r[r.kit_id.isin(key.kit_id)].drop_duplicates("kit_id", keep="last")
    r["grade"] = r["grade"].str.strip()
    r["placement"] = r["placement"].str.strip().str.lower()
    bad_codes = {k: sorted(codeset(v) - set(CODES)) for k, v in zip(r.kit_id, r.failure_codes) if codeset(v) - set(CODES)}
    issues["invalid_failure_codes"] = bad_codes
    r["h_codes"] = [";".join(sorted(codeset(v) & set(CODES))) for v in r.failure_codes]
    g = pd.to_numeric(r["grade"].str.replace(",", ".", regex=False), errors="coerce")     # "1", "1.0", "1,0" -> 1
    r["h_grade"] = g.where(g.isin([0, 1, 2]))
    r["h_place"] = r["placement"].where(r["placement"].isin(["pass", "fail"]))
    issues["missing_or_invalid_grade"] = sorted(set(key.kit_id) - set(r.kit_id[r.h_grade.notna()]))
    issues["missing_or_invalid_placement"] = sorted(set(key.kit_id) - set(r.kit_id[r.h_place.notna()]))
    return r[["kit_id", "h_grade", "h_place", "h_codes", "note", "timestamp"]], issues


# ------------------------------------------------------------------ Part A
def part_a(it, dup, nboot, rng, ai_rep=None):
    out = {}
    a = it[it.h_grade.notna()].copy()
    a["h"] = a.h_grade.astype(int)
    a["a"] = a.ai_grade.astype(int)
    out["n_rated_primary_items"] = int(len(a))
    out["n_by_stratum"] = a.stratum.value_counts().to_dict()
    T = np.zeros((3, 3), int)                      # robust to empty categories (review item M4)
    np.add.at(T, (a.h.to_numpy(int), a.a.to_numpy(int)), 1)
    tab = pd.DataFrame(T, index=[f"human_{i}" for i in range(3)], columns=[f"AI_{i}" for i in range(3)])
    out["confusion_rows_human_cols_AI"] = tab.astype(int).to_dict(orient="index")
    out["_tab"] = tab
    for lab, wcol in (("sample_unweighted", None), ("design_weighted_to_180", "sampling_weight")):
        w = None if wcol is None else a[wcol].to_numpy(float)
        kl, _ = wkappa(a.h, a.a, 3, "linear", w)
        kq, _ = wkappa(a.h, a.a, 3, "quadratic", w)
        kb, _ = wkappa((a.h == 2).astype(int), (a.a == 2).astype(int), 2, "nominal", w)
        ww = np.ones(len(a)) if w is None else w
        exact = float((ww * (a.h == a.a)).sum() / ww.sum())
        res = {"kappa_linear": kl, "kappa_quadratic": kq, "kappa_binary_fail_vs_pass": kb, "exact_agreement": exact}
        fl = lambda d, wc=wcol: wkappa(d.h, d.a, 3, "linear", None if wc is None else d[wc].to_numpy(float))[0]
        fb = lambda d, wc=wcol: wkappa((d.h == 2).astype(int), (d.a == 2).astype(int), 2, "nominal",
                                       None if wc is None else d[wc].to_numpy(float))[0]
        res["kappa_linear_boot"] = boot_ci(a, fl, nboot, rng)
        res["kappa_binary_boot"] = boot_ci(a, fb, nboot, rng)
        if wcol:
            res.update(binmetrics(a.h == 2, a.a == 2, a[wcol]))
            for m in ("sensitivity_P(AI2|human2)", "specificity_P(AI<2|human<2)"):
                res[m + "_boot"] = boot_ci(a, lambda d, m=m: binmetrics(d.h == 2, d.a == 2, d[wcol])[m], nboot, rng)
        out[lab] = res
    # within strata
    f = a[a.stratum == "AI_fail_all"]
    p = a[a.stratum == "AI_pass_sample"]
    k2 = int((f.h == 2).sum())
    out["AI_fail_stratum"] = {"n": int(len(f)), "kappa": "undefined (AI grade constant = 2 in this stratum)",
                              "human_grade_counts": {str(g): int((f.h == g).sum()) for g in (0, 1, 2)},
                              "confirmation_rate_P(human=2)": k2 / len(f) if len(f) else float("nan"),
                              "confirmation_rate_wilson95": wilson(k2, len(f)),
                              "linear_weighted_agreement": float(1 - (np.abs(f.h - f.a) / 2).mean()) if len(f) else float("nan")}
    kp, _ = wkappa(p.h, p.a, 3, "linear") if len(p) else (float("nan"), None)
    m2 = int((p.h == 2).sum())
    out["AI_pass_stratum"] = {"n": int(len(p)), "kappa_linear": kp,
                              "kappa_linear_boot": boot_ci(p, lambda d: wkappa(d.h, d.a, 3, "linear")[0], nboot, rng)
                              if len(p) else None,
                              "exact_agreement": float((p.h == p.a).mean()) if len(p) else float("nan"),
                              "human_grade_counts": {str(g): int((p.h == g).sum()) for g in (0, 1, 2)},
                              "missed_fail_rate_P(human=2)": m2 / len(p) if len(p) else float("nan"),
                              "missed_fail_wilson95": wilson(m2, len(p)),
                              "implied_human_fails_among_119_AI_pass": [N_AI_PASS_POOL * x for x in
                                                                        [m2 / max(1, len(p))] + wilson(m2, len(p))]}
    # by forced-review flag (review items M1, M2). The flag is a strong cue for the AI grade (21 of the 22 flagged
    # primary items are AI grade 2), so: flag = 0 -> kappa with bootstrap CI, co-reported with the primary kappa;
    # flag = 1 -> confirmation rate only (kappa is undefined or meaningless when one rater is (almost) constant).
    out["by_forced_review_flag"] = {}
    for fv, g in a.groupby("forced_review"):
        d = {"n": int(len(g)), "n_AI_grade2": int((g.a == 2).sum()),
             "human_grade_counts": {str(x): int((g.h == x).sum()) for x in (0, 1, 2)},
             "exact_agreement": float((g.h == g.a).mean())}
        if int(fv) == 0:
            ws = g.sampling_weight.to_numpy(float)
            d["kappa_linear_sample"] = wkappa(g.h, g.a, 3, "linear")[0]
            d["kappa_linear_design_weighted"] = wkappa(g.h, g.a, 3, "linear", ws)[0]
            d["kappa_binary_design_weighted"] = wkappa((g.h == 2).astype(int), (g.a == 2).astype(int), 2, "nominal", ws)[0]
            d["exact_agreement_design_weighted"] = float((ws * (g.h == g.a)).sum() / ws.sum())
            d["kappa_linear_sample_boot"] = boot_ci(g, lambda x: wkappa(x.h, x.a, 3, "linear")[0], nboot, rng)
            d["kappa_linear_design_weighted_boot"] = boot_ci(
                g, lambda x: wkappa(x.h, x.a, 3, "linear", x.sampling_weight.to_numpy(float))[0], nboot, rng)
        else:
            k_all = int((g.h == 2).sum())
            g2 = g[g.a == 2]
            k2 = int((g2.h == 2).sum())
            d["kappa"] = "not reported (pre-specified 2026-10-07): the flag is a cue for the AI grade"
            d["P(human=2)"] = k_all / len(g)
            d["P(human=2)_wilson95"] = wilson(k_all, len(g))
            d["confirmation_rate_P(human=2|AI=2)"] = k2 / len(g2) if len(g2) else float("nan")
            d["confirmation_rate_wilson95"] = wilson(k2, len(g2))
            d["AI_pass_items_human_grades"] = {str(x): int((g[g.a < 2].h == x).sum()) for x in (0, 1, 2)}
        out["by_forced_review_flag"][str(int(fv))] = d
    # failure codes
    pc = {}
    hc = [codeset(x) for x in a.h_codes]
    ac = [codeset(x) for x in a.ai_failure_codes]
    for c in CODES:
        hx = np.array([c in s for s in hc], int)
        ax = np.array([c in s for s in ac], int)
        pc[c] = {"both": int((hx & ax).sum()), "human_only": int((hx & (1 - ax)).sum()),
                 "AI_only": int(((1 - hx) & ax).sum()),
                 "kappa": wkappa(hx, ax, 2, "nominal")[0] if len(set(hx)) > 1 and len(set(ax)) > 1 else float("nan")}
    out["failure_codes_presence"] = pc
    # placement
    b = it[it.h_place.notna()]
    hp, ap = (b.h_place == "fail").astype(int), (b.ai_placement == "fail").astype(int)
    out["placement"] = {"n": int(len(b)), "human_fail": int(hp.sum()), "AI_fail": int(ap.sum()),
                        "agreement": float((hp == ap).mean()) if len(b) else float("nan"),
                        "kappa": (wkappa(hp, ap, 2, "nominal")[0] if hp.nunique() > 1 and ap.nunique() > 1
                                  else float("nan")),
                        "note": "kappa is reported as NA when either rater is constant (the AI failed placement in none of 180)",
                        "human_fail_kit_ids": b.kit_id[hp == 1].tolist()}
    # intra-rater
    pairs = []
    for _, d in dup.iterrows():
        o = it[it.kit_id == d.is_duplicate_of].iloc[0]
        pairs.append({"original": o.kit_id, "duplicate": d.kit_id, "g1": o.h_grade, "g2": d.h_grade,
                      "p1": o.h_place, "p2": d.h_place, "c1": o.h_codes, "c2": d.h_codes})
    P = pd.DataFrame(pairs)
    q = P[P.g1.notna() & P.g2.notna()]
    qp = P[P.p1.notna() & P.p2.notna()]
    out["intra_rater_duplicates"] = {
        "n_pairs": int(len(P)), "n_pairs_with_both_grades": int(len(q)),
        "kappa_linear": wkappa(q.g1.astype(int), q.g2.astype(int), 3, "linear")[0] if len(q) else float("nan"),
        "exact_agreement": float((q.g1 == q.g2).mean()) if len(q) else float("nan"),
        "placement_agreement": float((qp.p1 == qp.p2).mean()) if len(qp) else float("nan"),
        # L4: only pairs in which both items carry a grade (an unrated item is not 'identical')
        "codes_identical": float((q.c1.fillna("") == q.c2.fillna("")).mean()) if len(q) else float("nan"),
        "pairs": P.to_dict(orient="records"),
        "note": "10 pairs: kappa is imprecise; exact agreement is the main intra-rater figure"}
    # L5: secondary comparator = AI pass-1 grade (the final AI grade is the worse of pass 1 and pass 2 for the 36
    # sessions the AI rated twice; the human looks once); and the AI's own repeatability next to the human's.
    if "ai_grade_pass1" in a.columns and a.ai_grade_pass1.notna().all():
        a1 = a.ai_grade_pass1.astype(int)
        wv = a.sampling_weight.to_numpy(float)
        out["vs_AI_pass1_grade"] = {
            "n": int(len(a)), "AI_pass1_grade_counts": {str(x): int((a1 == x).sum()) for x in (0, 1, 2)},
            "kappa_linear_sample": wkappa(a.h, a1, 3, "linear")[0],
            "kappa_linear_design_weighted": wkappa(a.h, a1, 3, "linear", wv)[0],
            "exact_agreement": float((a.h == a1).mean()),
            "kappa_linear_design_weighted_boot": boot_ci(
                a, lambda x: wkappa(x.h, x.ai_grade_pass1.astype(int), 3, "linear",
                                    x.sampling_weight.to_numpy(float))[0], nboot, rng),
            "note": "secondary comparator (review L5); design weights still refer to the final-grade strata"}
    if ai_rep is not None:
        out["AI_repeatability_pass1_vs_pass2"] = ai_rep
    return out


def ai_repeatability(p1path, p2path):
    """AI pass-1 vs pass-2 grade on the sessions the AI rated twice (grade columns only; notes are never read)."""
    try:
        p1 = pd.read_csv(p1path, usecols=["qc_id", "grade"])
        p2 = pd.read_csv(p2path, usecols=["qc_id", "grade"])
    except Exception as e:      # pragma: no cover
        return {"error": f"not available: {e}"}
    m = p2.merge(p1, on="qc_id", suffixes=("_pass2", "_pass1")).dropna()
    g1, g2 = m.grade_pass1.astype(int), m.grade_pass2.astype(int)
    return {"n_sessions_rated_twice": int(len(m)), "identical": int((g1 == g2).sum()),
            "exact_agreement": float((g1 == g2).mean()) if len(m) else float("nan"),
            "kappa_linear": wkappa(g1, g2, 3, "linear")[0] if len(m) else float("nan"),
            "crosstab_pass1_rows_pass2_cols": {f"pass1_{i}": {f"pass2_{j}": int(((g1 == i) & (g2 == j)).sum())
                                                              for j in (0, 1, 2)} for i in (0, 1, 2)}}


# ------------------------------------------------------------------ Part B
def reconstruct(mg, grade, place):
    """Registered first-usable rule (build_usable.py) with the scenario's grade/placement."""
    rows = []
    for i, r in mg.iterrows():
        comp = set(r.reason_set) - {"c_qc_grade", "d_placement"}
        g, pl = grade[i], place[i]
        if not (g == g and int(g) in (0, 1)):
            comp.add("c_qc_grade")
        if pl != "pass":
            comp.add("d_placement")
        if not comp and bool(r["_posthoc_plaus_fail"]):     # would become usable but fails the frozen re-check
            comp.add("posthoc_snvta_plausibility_recheck")
        rows.append(comp)
    usable = np.array([int(not c) for c in rows])
    first = np.zeros(len(mg), int)
    ok = pd.DataFrame({"code": mg.code.values, "ses": mg.ses_order.values, "u": usable})
    fmin = ok[ok.u == 1].groupby("code").ses.min()
    for j, (c, s) in enumerate(zip(mg.code, mg.ses_order)):
        first[j] = int(c in fmin.index and fmin[c] == s)
    reason = [";".join([x for x in REASONS + ["posthoc_snvta_plausibility_recheck"] if x in c]) for c in rows]
    return usable, first, reason


def part_b(it, dup, args, outdir, label):
    res = {"label": "POST HOC sensitivity: not part of the registered confirmatory analysis"}
    mg = pd.read_csv(args.merged, dtype={"code": str}, keep_default_na=True)
    gf = pd.read_csv(args.grades, dtype={"code": str})
    mg = mg.merge(gf[["code", "ses_order", "qc_grade", "placement"]].rename(
        columns={"qc_grade": "ai_grade_gf", "placement": "ai_place"}), on=["code", "ses_order"], how="left")
    same = ((mg.qc_grade == mg.ai_grade_gf) | (mg.qc_grade.isna() & mg.ai_grade_gf.isna())).all()
    if not same:
        sys.exit("merged-table qc_grade differs from grades_final.csv: stop")
    mg["reason_set"] = [set(x.split(";")) if isinstance(x, str) and x else set() for x in mg.fail_reason]
    unknown = set().union(*mg.reason_set) - set(REASONS)
    if unknown:
        sys.exit(f"unknown fail_reason components {unknown}: stop")
    # frozen SN-VTA re-verification (apply_plausibility) for any session that could become usable
    mg["_posthoc_plaus_fail"] = ((mg.chi_SNVTA < -50) | (mg.chi_SNVTA > 300) | (mg.fimp_SNVTA > 0.05) |
                                 (mg.inmask_SNVTA < 0.95) | mg.chi_SNVTA.isna() | mg.roi_vol.isna() |
                                 mg.Q_phase.isna()).fillna(False)
    ai_g = mg.ai_grade_gf.to_numpy(float)
    ai_p = mg.ai_place.fillna("missing").to_numpy(object)
    u0, f0, r0 = reconstruct(mg, ai_g, ai_p)
    repro = {"usable_equal": bool((u0 == mg.usable.to_numpy()).all()),
             "first_usable_equal": bool((f0 == mg.first_usable.to_numpy()).all()),
             "fail_reason_equal": bool(all(set(a.split(";")) - {""} == b for a, b in zip(r0, mg.reason_set)))}
    res["reconstruction_check_AI_frozen"] = repro
    if not all(repro.values()):
        sys.exit(f"the first-usable reconstruction does not reproduce the frozen table: {repro}")
    # human per-session rating: conservative over original + duplicate
    allr = pd.concat([it, dup], ignore_index=True)
    hs = {}
    for (c, s), g in allr.groupby(["code", "ses_order"]):
        gg = g[g.h_grade.notna()]
        pp = g[g.h_place.notna()]
        hs[(c, int(s))] = {"h_grade": float(gg.h_grade.max()) if len(gg) else float("nan"),
                           "h_place": ("fail" if (pp.h_place == "fail").any() else "pass") if len(pp) else None,
                           "n_ratings": int(len(g)), "kit_ids": ";".join(g.kit_id)}
    keys = list(zip(mg.code, mg.ses_order.astype(int)))
    h_g = np.array([hs[k]["h_grade"] if k in hs else np.nan for k in keys], float)
    h_p = np.array([hs[k]["h_place"] if k in hs else None for k in keys], object)
    rated_g = ~np.isnan(h_g)
    rated_p = np.array([x is not None for x in h_p])
    res["n_sessions_in_kit"] = int(sum(k in hs for k in keys))
    res["n_sessions_with_human_grade"] = int(rated_g.sum())
    res["n_sessions_with_human_placement"] = int(rated_p.sum())
    res["unrated_kit_sessions_keep_AI"] = int(sum(k in hs for k in keys) - rated_g.sum())
    # M3: coverage of the sensitivity analysis (sessions outside the kit always keep the AI grade)
    inkit = np.array([k in hs for k in keys])
    res["coverage"] = {
        "frozen_first_usable_sessions": int(mg.first_usable.sum()),
        "frozen_first_usable_not_in_kit": int(((mg.first_usable == 1).to_numpy() & ~inkit).sum()),
        "AI_pass_sessions": int((ai_g < 2).sum()),
        "AI_pass_sessions_not_in_kit": int(((ai_g < 2) & ~inkit).sum()),
        "note": "sessions outside the kit keep the AI grade in every scenario; the claim supported is robustness to "
                "replacing the AI grade by the human grade for the re-rated sessions only (review item M3)"}
    sc = {}
    sc["AI_frozen"] = (ai_g, ai_p)
    sc["H_replace"] = (np.where(rated_g, h_g, ai_g), np.where(rated_p, h_p, ai_p))
    sc["H_grade_only"] = (np.where(rated_g, h_g, ai_g), ai_p)
    sc["Strict_union"] = (np.where(rated_g, np.fmax(h_g, ai_g), ai_g),
                          np.where(rated_p & (h_p == "fail"), "fail", ai_p))
    sc["Lenient"] = (np.where(rated_g, np.fmin(h_g, ai_g), ai_g),
                     np.where(rated_p & (h_p == "pass") & (ai_p != "missing"), "pass", ai_p))
    tab = pd.DataFrame({"sess_key": mg.participant_id + "_" + mg.session_id, "code": mg.code,
                        "ses_order": mg.ses_order, "ai_grade": ai_g, "human_grade": h_g,
                        "ai_placement": ai_p, "human_placement": h_p,
                        "kit_ids": [hs[k]["kit_ids"] if k in hs else "" for k in keys]})
    summ = {}
    base_first = set(tab.sess_key[f0 == 1])
    for name, (g, p) in sc.items():
        u, f, r = reconstruct(mg, g, p)
        tab[f"usable__{name}"] = u
        tab[f"first_usable__{name}"] = f
        tab[f"fail_reason__{name}"] = r
        fs = set(tab.sess_key[f == 1])
        subj0 = {k.rsplit("_ses-", 1)[0] for k in base_first}
        subj1 = {k.rsplit("_ses-", 1)[0] for k in fs}
        summ[name] = {"n_usable_sessions": int(u.sum()), "n_first_usable": int(f.sum()),
                      "subjects_added": sorted(subj1 - subj0), "subjects_removed": sorted(subj0 - subj1),
                      "subjects_with_changed_first_session": sorted(
                          {k.rsplit("_ses-", 1)[0] for k in (fs ^ base_first)} & subj0 & subj1),
                      "posthoc_plausibility_recheck_excluded": tab.sess_key[[
                          "posthoc_snvta_plausibility_recheck" in x for x in r]].tolist()}
    res["sample_changes_vs_frozen"] = summ
    scen_csv = os.path.join(outdir, f"{label}scenario_usable.csv")
    tab.to_csv(scen_csv, index=False)
    # R refit with the frozen helpers
    rj = os.path.join(outdir, f"{label}m1_refit.json")
    cmd = [args.rscript, os.path.join(HERE, "refit_m1_human_qc.R"), "--merged", args.merged, "--scen", scen_csv,
           "--gates", args.gates, "--rdir", args.rdir, "--out", rj, "--bound", repr(FROZEN_BOUND)]
    pr = subprocess.run(cmd, capture_output=True, text=True)
    res["rscript_returncode"] = pr.returncode
    if pr.returncode != 0:
        res["rscript_stderr"] = pr.stderr[-4000:]
        return res, None
    fit = json.load(open(rj))
    ai = fit["scenarios"].get("AI_frozen", {})
    res["refit_reproduces_frozen"] = bool(abs(ai.get("beta_PS_ppb", np.nan) - FROZEN_BETA) < 1e-8 and ai.get("n") == 80)
    if not res["refit_reproduces_frozen"]:
        res["WARNING"] = "AI_frozen refit does not reproduce the confirmatory beta/n: do not interpret the scenarios"
    if args.nextrap > 0:
        res["extrapolation_unrated_AI_pass"] = extrapolate(mg, sc["H_replace"], ai_g, inkit, rated_g, h_g, args,
                                                           outdir, label)
    return res, fit


def extrapolate(mg, h_replace, ai_g, inkit, rated_g, h_g, args, outdir, label):
    """M3 (optional item 3 of the review), POST HOC, pre-specified 2026-10-07 before any rating existed.

    Starts from H_replace and imputes human fails among the AI-pass sessions that are NOT in the kit (they otherwise
    always keep the AI grade). Per draw, a missed-fail rate is drawn from the Jeffreys posterior Beta(m + 1/2,
    n - m + 1/2), where m of n = rated AI-pass kit sessions got human grade 2 (session level, worse of original and
    repeat); each unrated AI-pass session then fails with that rate. Placement is kept. M1 is refitted per draw with
    the same frozen R helpers. Reported: share of draws with each descriptive label, and ranges of n, beta and p_B."""
    g0, p0 = h_replace
    ap = rated_g & inkit & (ai_g < 2)
    m, n = int((h_g[ap] == 2).sum()), int(ap.sum())
    unr = (~inkit) & (ai_g < 2)
    rng = np.random.default_rng(args.extrap_seed)
    cols, meta = {"sess_key": (mg.participant_id + "_" + mg.session_id).to_numpy()}, []
    for d in range(args.nextrap):
        rate = float(rng.beta(m + 0.5, n - m + 0.5))
        hit = unr & (rng.random(len(g0)) < rate)
        g = np.where(hit, 2.0, g0)
        u, f, _ = reconstruct(mg, g, p0)
        cols[f"usable__x{d}"], cols[f"first_usable__x{d}"] = u, f
        meta.append({"draw": f"x{d}", "rate": rate, "n_imputed_fails": int(hit.sum())})
    scen = os.path.join(outdir, f"{label}extrapolation_scenarios.csv")
    pd.DataFrame(cols).to_csv(scen, index=False)
    rj = os.path.join(outdir, f"{label}m1_refit_extrapolation.json")
    cmd = [args.rscript, os.path.join(HERE, "refit_m1_human_qc.R"), "--merged", args.merged, "--scen", scen,
           "--gates", args.gates, "--rdir", args.rdir, "--out", rj, "--bound", repr(FROZEN_BOUND)]
    pr = subprocess.run(cmd, capture_output=True, text=True)
    out = {"label": "POST HOC extrapolation (review item M3): human fails imputed among unrated AI-pass sessions",
           "n_draws": args.nextrap, "seed": args.extrap_seed, "observed_missed_fails_m": m, "rated_AI_pass_sessions_n": n,
           "unrated_AI_pass_sessions": int(unr.sum()), "rscript_returncode": pr.returncode}
    if pr.returncode != 0:
        out["rscript_stderr"] = pr.stderr[-4000:]
        return out
    fit = json.load(open(rj))["scenarios"]
    rows = []
    for x in meta:
        s = fit.get(x["draw"], {})
        if "error" in s or not s:
            rows.append(dict(x, error=s.get("error", "missing")))
            continue
        rows.append(dict(x, n=s["n"], beta_PS_ppb=s["beta_PS_ppb"], se_hc3=s["se_hc3"],
                         pB_fixed=s["bound_fixed"]["pB_one_sided"], pB_rederived=s["bound_rederived"]["pB"],
                         A=bool(s["A_beta_negative_p05"]), label=s["outcome_label_descriptive"]))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, f"{label}extrapolation_draws.csv"), index=False)
    ok = df[df["error"].isna()] if "error" in df.columns else df
    q = lambda c: [float(ok[c].min()), float(ok[c].quantile(0.025)), float(ok[c].median()),
                   float(ok[c].quantile(0.975)), float(ok[c].max())] if len(ok) else None
    out.update({"n_draws_fitted": int(len(ok)), "n_draws_error": int(len(df) - len(ok)),
                "label_counts": ok.label.astype(str).value_counts().to_dict() if len(ok) else {},
                "share_label_B": float((ok.label.astype(str) == "B").mean()) if len(ok) else float("nan"),
                "share_A": float(ok.A.mean()) if len(ok) else float("nan"),
                "share_pB_fixed_below_05": float((ok.pB_fixed < 0.05).mean()) if len(ok) else float("nan"),
                "n_imputed_fails_min_q025_median_q975_max": q("n_imputed_fails"),
                "n_min_q025_median_q975_max": q("n"), "beta_min_q025_median_q975_max": q("beta_PS_ppb"),
                "pB_fixed_min_q025_median_q975_max": q("pB_fixed")})
    return out


def clean(o):
    """JSON-safe copy: numpy scalars -> Python, NaN/inf -> None, sets -> sorted lists."""
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, set):
        return sorted(clean(v) for v in o)
    if hasattr(o, "item") and not isinstance(o, (str, bytes)):
        o = o.item()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


# ------------------------------------------------------------------ report
def fmt(x, d=3):
    try:
        return "NA" if x is None or x != x else f"{x:.{d}f}"
    except TypeError:
        return str(x)


def report(A, B, fit, meta, path):
    S = meta["simulated"]
    L = []
    if S:
        L += ["# SIMULATED TEST RUN - NOT REAL RATINGS", "",
              "**Every number below comes from a simulated ratings file used to test the pipeline. Do not cite.**", ""]
    L += ["# LEAD-chi: blinded human visual-QC validation of the AI rater (POST HOC)", "",
          "Status: **post hoc / exploratory**. Not part of the registered confirmatory analysis; H1 is not re-decided.", "",
          f"- ratings file: `{meta['ratings']}` (sha256 `{meta['ratings_sha256']}`)",
          f"- kit key sha256: `{meta['key_sha256']}`; run (UTC): {meta['utc']}; bootstrap {meta['nboot']} (seed {meta['seed']})",
          f"- input issues: {json.dumps(meta['issues'])}", "",
          "## A. Agreement, human vs AI (primary items)", "",
          f"Rated primary items: {A['n_rated_primary_items']} ({A['n_by_stratum']}).", "",
          "Confusion matrix (rows human, columns AI):", "", "| | AI 0 | AI 1 | AI 2 |", "|---|---|---|---|"]
    t = A["_tab"]
    for i in range(3):
        L.append(f"| human {i} | " + " | ".join(str(int(t.iloc[i, j])) for j in range(3)) + " |")
    L += ["", "| estimand | linear-weighted kappa [95% CI] | binary fail-vs-pass kappa [95% CI] | exact agreement |",
          "|---|---|---|---|"]
    for lab in ("sample_unweighted", "design_weighted_to_180"):
        r = A[lab]
        L.append(f"| {lab} | {fmt(r['kappa_linear'])} [{fmt(r['kappa_linear_boot']['ci95'][0])}, "
                 f"{fmt(r['kappa_linear_boot']['ci95'][1])}] | {fmt(r['kappa_binary_fail_vs_pass'])} "
                 f"[{fmt(r['kappa_binary_boot']['ci95'][0])}, {fmt(r['kappa_binary_boot']['ci95'][1])}] | "
                 f"{fmt(r['exact_agreement'])} |")
    f0 = A["by_forced_review_flag"].get("0")
    if f0 and "kappa_linear_sample_boot" in f0:      # M2: co-reported with the same prominence as the primary kappa
        for lab, k, b, ex in (("flag0_only_sample (n %d)" % f0["n"], f0["kappa_linear_sample"],
                               f0["kappa_linear_sample_boot"], f0["exact_agreement"]),
                              ("flag0_only_design_weighted (n %d)" % f0["n"], f0["kappa_linear_design_weighted"],
                               f0["kappa_linear_design_weighted_boot"], f0["exact_agreement_design_weighted"])):
            kb = f0["kappa_binary_design_weighted"] if "design" in lab else None
            L.append(f"| {lab} | {fmt(k)} [{fmt(b['ci95'][0])}, {fmt(b['ci95'][1])}] | "
                     f"{fmt(kb) if kb is not None else '-'} | {fmt(ex)} |")
    L += ["", "Pre-specified (README_coordinator section 7, 2026-10-07): the design-weighted linear kappa (primary) and "
          "the design-weighted linear kappa on forced-review flag = 0 items are reported with equal prominence; the "
          "flag is shown to both raters and is a strong cue for the AI grade. Kappa is NA when a rater is constant."]
    dw = A["design_weighted_to_180"]
    f, p = A["AI_fail_stratum"], A["AI_pass_stratum"]
    L += ["", "Design-weighted (to the 180 graded sessions), human as reference for fail = grade 2: "
          f"AI sensitivity {fmt(dw['sensitivity_P(AI2|human2)'])} "
          f"[{fmt(dw['sensitivity_P(AI2|human2)_boot']['ci95'][0])}, {fmt(dw['sensitivity_P(AI2|human2)_boot']['ci95'][1])}], "
          f"specificity {fmt(dw['specificity_P(AI<2|human<2)'])} "
          f"[{fmt(dw['specificity_P(AI<2|human<2)_boot']['ci95'][0])}, {fmt(dw['specificity_P(AI<2|human<2)_boot']['ci95'][1])}], "
          f"PPV {fmt(dw['PPV_P(human2|AI2)'])}, NPV {fmt(dw['NPV_P(human<2|AI<2)'])}; "
          f"estimated human fails among 180: {fmt(dw['est_human_fail_n_of_180'], 1)}.", "",
          f"- AI-fail stratum (n = {f['n']}; kappa undefined, AI constant): human confirms grade 2 in "
          f"{fmt(f['confirmation_rate_P(human=2)'])} [Wilson {fmt(f['confirmation_rate_wilson95'][0])}, "
          f"{fmt(f['confirmation_rate_wilson95'][1])}]; human grades {f['human_grade_counts']}.",
          f"- AI-pass stratum (n = {p['n']}): linear kappa {fmt(p['kappa_linear'])}; exact agreement "
          f"{fmt(p['exact_agreement'])}; human grade 2 (AI missed fail) {fmt(p['missed_fail_rate_P(human=2)'])} "
          f"[Wilson {fmt(p['missed_fail_wilson95'][0])}, {fmt(p['missed_fail_wilson95'][1])}], i.e. about "
          f"{fmt(p['implied_human_fails_among_119_AI_pass'][0], 1)} of the 119 AI-pass sessions.",
          "- By forced-review flag: " + "; ".join(
              (f"flag 0: n {v['n']} ({v['n_AI_grade2']} AI grade 2), linear kappa {fmt(v['kappa_linear_sample'])} "
               f"(design-weighted {fmt(v['kappa_linear_design_weighted'])})") if k == "0" else
              (f"flag 1: n {v['n']} ({v['n_AI_grade2']} AI grade 2), reported as a confirmation rate only: "
               f"P(human = 2 | AI = 2) {fmt(v['confirmation_rate_P(human=2|AI=2)'])} [Wilson "
               f"{fmt(v['confirmation_rate_wilson95'][0])}, {fmt(v['confirmation_rate_wilson95'][1])}], "
               f"P(human = 2) {fmt(v['P(human=2)'])}; human grades of the flagged AI-pass items "
               f"{v['AI_pass_items_human_grades']}")
              for k, v in A["by_forced_review_flag"].items()),
          f"- Placement: human fail {A['placement']['human_fail']} of {A['placement']['n']} items, AI fail "
          f"{A['placement']['AI_fail']}; agreement {fmt(A['placement']['agreement'])}; kappa {fmt(A['placement']['kappa'])}.",
          "- Failure-code presence (both / human only / AI only; kappa): " + "; ".join(
              f"{c} {v['both']}/{v['human_only']}/{v['AI_only']} ({fmt(v['kappa'], 2)})"
              for c, v in A["failure_codes_presence"].items()), ""]
    if "vs_AI_pass1_grade" in A:
        v = A["vs_AI_pass1_grade"]
        L += [f"- Secondary comparator, AI pass-1 grade (the final AI grade is the worse of pass 1 and pass 2 for the "
              f"sessions the AI rated twice; the human looks once): AI pass-1 grades {v['AI_pass1_grade_counts']}; "
              f"linear kappa {fmt(v['kappa_linear_sample'])} (design-weighted {fmt(v['kappa_linear_design_weighted'])} "
              f"[{fmt(v['kappa_linear_design_weighted_boot']['ci95'][0])}, "
              f"{fmt(v['kappa_linear_design_weighted_boot']['ci95'][1])}]); exact agreement {fmt(v['exact_agreement'])}.", ""]
    if A.get("excluded_from_part_A"):
        L += [f"- Excluded from Part A (pre-specified rule): {A['excluded_from_part_A']}.", ""]
    d = A["intra_rater_duplicates"]
    L += ["## Intra-rater reliability (hidden duplicates)", "",
          f"{d['n_pairs_with_both_grades']} of {d['n_pairs']} pairs rated twice: exact grade agreement "
          f"{fmt(d['exact_agreement'])}, linear-weighted kappa {fmt(d['kappa_linear'])}, placement agreement "
          f"{fmt(d['placement_agreement'])}, identical code sets (pairs with both grades) {fmt(d['codes_identical'])}. "
          f"{d['note']}. Repeats were 21-70 items apart (about 15-50 min), so recall is possible and the figure may "
          "be optimistic.", ""]
    ar = A.get("AI_repeatability_pass1_vs_pass2")
    if ar and "error" not in ar:
        L += [f"For comparison, the AI's own repeatability: {ar['identical']} of {ar['n_sessions_rated_twice']} "
              f"re-rated sessions got the same grade in pass 1 and pass 2 (exact {fmt(ar['exact_agreement'])}, "
              f"linear kappa {fmt(ar['kappa_linear'])}).", ""]
    L += ["## B. POST HOC sensitivity: first-usable sample re-derived with the human grade, M1 refitted", "",
          f"Reconstruction of the frozen first-usable list from fail_reason + AI grades: {B['reconstruction_check_AI_frozen']}.",
          f"Sessions in kit: {B['n_sessions_in_kit']}; with a human grade: {B['n_sessions_with_human_grade']}.", ""]
    if fit is None:
        L += [f"R refit FAILED (return code {B['rscript_returncode']}):", "", "```", B.get("rscript_stderr", ""), "```"]
    else:
        L += [f"Refit reproduces the frozen confirmatory estimate under AI_frozen: **{B['refit_reproduces_frozen']}**.", "",
              "| scenario | n | cells TD/ADHD/PS/PS+ADHD | beta_PS ppb (HC3 SE) | 95% CI | p two-sided | "
              "bound test at -6.1318 ppb: t_B, p_B | re-derived bound (ppb): p_B | A | descriptive label |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for name, s in fit["scenarios"].items():
            if "error" in s:
                L.append(f"| {name} | error: {s['error']} |||||||||")
                continue
            c = s["cells"]
            bf, br = s["bound_fixed"], s["bound_rederived"]
            L.append(f"| {name} | {s['n']} | {c['TD']}/{c['ADHD-only']}/{c['PS-only']}/{c['PS+ADHD']} | "
                     f"{fmt(s['beta_PS_ppb'], 2)} ({fmt(s['se_hc3'], 2)}) | [{fmt(s['ci95_ppb'][0], 2)}, "
                     f"{fmt(s['ci95_ppb'][1], 2)}] | {fmt(s['p_two'], 4)} | {fmt(bf['tB'], 2)}, {fmt(bf['pB_one_sided'], 5)} | "
                     f"{fmt(br['bound_ppb'], 2)}: {fmt(br['pB'], 5)} | {s['A_beta_negative_p05']} | {s['outcome_label_descriptive']} |")
        L += ["", "Sample changes relative to the frozen first-usable list:", ""]
        for name, s in B["sample_changes_vs_frozen"].items():
            L.append(f"- {name}: usable sessions {s['n_usable_sessions']}, first-usable {s['n_first_usable']}; "
                     f"subjects added {len(s['subjects_added'])}, removed {len(s['subjects_removed'])}, first session "
                     f"changed {len(s['subjects_with_changed_first_session'])}; post hoc plausibility re-check exclusions "
                     f"{len(s['posthoc_plausibility_recheck_excluded'])}.")
        L += ["", "Labels in the last column are the registered decision arithmetic applied descriptively to a post hoc "
              "sample; they do not replace the registered outcome (B).", ""]
        cv = B.get("coverage", {})
        if cv:
            L += [f"**Coverage (pre-specified wording, README_coordinator section 7.3).** {cv['frozen_first_usable_not_in_kit']}"
                  f" of the {cv['frozen_first_usable_sessions']} frozen first-usable sessions and "
                  f"{cv['AI_pass_sessions_not_in_kit']} of the {cv['AI_pass_sessions']} AI-pass sessions are not in the "
                  "kit and keep the AI grade in every scenario above. A B label under H_replace therefore means "
                  "'robust to replacing the AI grade by the human grade for the re-rated sessions', not 'independent "
                  "of the visual-QC rater'. Before any rating, a sweep of 300 random human gradings of the kit sessions "
                  "gave label B in 300 of 300, so the main value of this exercise is the validation of the AI (Part A).", ""]
        ex = B.get("extrapolation_unrated_AI_pass")
        if ex:
            if ex.get("rscript_returncode") != 0 or "share_label_B" not in ex:
                L += [f"Extrapolation to the unrated AI-pass sessions FAILED (R return code {ex.get('rscript_returncode')}).", ""]
            elif not ex.get("n_draws_fitted"):
                L += [f"**Extrapolation (POST HOC, pre-specified 2026-10-07).** None of the {ex['n_draws']} draws could be "
                      f"fitted (m = {ex['observed_missed_fails_m']} of n = {ex['rated_AI_pass_sessions_n']}; for example "
                      "no usable session left); see the draws CSV for the R error messages.", ""]
            else:
                L += [f"**Extrapolation (POST HOC, pre-specified 2026-10-07).** Starting from H_replace, human fails were "
                      f"imputed among the {ex['unrated_AI_pass_sessions']} AI-pass sessions outside the kit, at a rate "
                      f"drawn per replicate from Beta(m + 1/2, n - m + 1/2) with m = {ex['observed_missed_fails_m']} of "
                      f"n = {ex['rated_AI_pass_sessions_n']} rated AI-pass kit sessions failed by the human; "
                      f"{ex['n_draws_fitted']} of {ex['n_draws']} draws fitted (seed {ex['seed']}). Label counts "
                      f"{ex['label_counts']}; share B {fmt(ex['share_label_B'])}; share A {fmt(ex['share_A'])}; "
                      f"imputed fails (min, 2.5%, median, 97.5%, max) {[fmt(x, 1) for x in ex['n_imputed_fails_min_q025_median_q975_max']]}; "
                      f"n {[fmt(x, 0) for x in ex['n_min_q025_median_q975_max']]}; beta_PS ppb "
                      f"{[fmt(x, 2) for x in ex['beta_min_q025_median_q975_max']]}; p_B at the fixed bound "
                      f"{[fmt(x, 5) for x in ex['pB_fixed_min_q025_median_q975_max']]}.", ""]
    L += ["## Method notes", "",
          "- Kappa: weighted Cohen kappa, linear weights (registered G4 statistic), same formula as `qc_montage.py`; "
          "reported as NA when either rater uses a single category (post hoc tooling fix, review item M1).",
          "- Design weights: AI grade-2 sessions are a census (weight 1); the 40 AI grade-0/1 sessions represent 119 "
          "(weight 119/40). Bootstrap resamples within the two strata (percentile 95% CI).",
          "- Human per-session rating = worse of the original and its hidden duplicate (registered conservative rule).",
          "- Part B keeps every non-visual exclusion component of the frozen `fail_reason` unchanged, including the "
          "automatic grade 2 (`c_auto_grade2`), which no rater can override (rubric section 2).",
          "- M1 refit sources the frozen R helpers; bound test at the frozen registered bound -6.1318 ppb "
          "(t_B = (beta + 6.1318)/SE, one-sided upper p, df = n - p) and at the bound re-derived in each sample."]
    if S:
        L += ["", "**SIMULATED TEST RUN - NOT REAL RATINGS.**"]
    open(path, "w").write("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ratings", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--simulated", action="store_true", help="label every output as a SIMULATED test run")
    for k, v in DEF.items():
        ap.add_argument(f"--{k}", default=v)
    ap.add_argument("--nboot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261008)
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--nextrap", type=int, default=1000,
                    help="draws of the M3 extrapolation to the AI-pass sessions outside the kit (0 = skip)")
    ap.add_argument("--extrap-seed", dest="extrap_seed", type=int, default=20261009)
    ap.add_argument("--exclude-from-a", dest="exclude_from_a", default="",
                    help="comma-separated kit IDs dropped from Part A only, with any repeat of the same session "
                         "(pre-specified: hfaec73 if the rater has ever seen LEADchi/qc_ornekler/)")
    args = ap.parse_args()
    label = "SIMULATED_" if args.simulated else ""
    key = pd.read_csv(args.key, dtype={"code": str, "is_duplicate_of": str}, keep_default_na=False)
    key["ai_grade"] = key.ai_grade.astype(int)
    key["sampling_weight"] = key.sampling_weight.astype(float)
    key["forced_review"] = key.forced_review.astype(int)
    # L5: AI pass-1 grade as a secondary comparator (grade column only; the AI notes are never read)
    try:
        p1 = pd.read_csv(args.pass1, usecols=["qc_id", "grade"]).rename(columns={"grade": "ai_grade_pass1"})
        key = key.merge(p1, on="qc_id", how="left")
    except Exception as e:      # pragma: no cover
        print(f"note: AI pass-1 grades not available ({e}); secondary comparator skipped")
    # read and validate the ratings BEFORE anything is written (M4)
    r, issues = read_ratings(args.ratings, key)
    allit = key.merge(r, on="kit_id", how="left")
    it = allit[allit.is_duplicate_of == ""].reset_index(drop=True)
    dup = allit[allit.is_duplicate_of != ""].reset_index(drop=True)
    if not it.h_grade.notna().any():
        sys.exit(f"no primary item has a valid grade (0, 1 or 2) in {args.ratings} (delimiter read as "
                 f"{issues['file_delimiter']}, encoding {issues['file_encoding']}); nothing was written")
    excl = [x.strip().lower() for x in args.exclude_from_a.split(",") if x.strip()]
    unknown = sorted(set(excl) - set(key.kit_id))
    if unknown:
        sys.exit(f"--exclude-from-a: unknown kit IDs {unknown}; nothing was written")
    os.makedirs(args.out, exist_ok=True)
    if any(f.endswith("human_qc_results.json") for f in os.listdir(args.out)):
        sys.exit(f"{args.out} already holds results: use a new output folder (nothing is overwritten)")
    excl_q = set(key.qc_id[key.kit_id.isin(excl)])
    itA, dupA = it.copy(), dup.copy()
    for d in (itA, dupA):
        d.loc[d.qc_id.isin(excl_q), ["h_grade", "h_place"]] = np.nan
        d.loc[d.qc_id.isin(excl_q), "h_codes"] = ""
    issues["excluded_from_part_A"] = sorted(set(allit.kit_id[allit.qc_id.isin(excl_q)]))
    rng = np.random.default_rng(args.seed)
    A = part_a(itA, dupA, args.nboot, rng, ai_repeatability(args.pass1, args.pass2))
    A["excluded_from_part_A"] = issues["excluded_from_part_A"]
    tab = A["_tab"]
    tab.to_csv(os.path.join(args.out, f"{label}confusion_human_rows_vs_AI_cols.csv"))
    allit.drop(columns=["note"]).to_csv(os.path.join(args.out, f"{label}per_item_human_vs_AI.csv"), index=False)
    B, fit = part_b(it, dup, args, args.out, label)
    meta = {"simulated": args.simulated, "ratings": os.path.abspath(args.ratings), "ratings_sha256": sha256(args.ratings),
            "key_sha256": sha256(args.key), "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "nboot": args.nboot, "seed": args.seed, "issues": {k: (len(v) if isinstance(v, (list, dict)) else v)
                                                              for k, v in issues.items()},
            "nextrap": args.nextrap, "extrap_seed": args.extrap_seed, "exclude_from_part_A": excl,
            "script_sha256": sha256(os.path.abspath(__file__)),
            "refit_script_sha256": sha256(os.path.join(HERE, "refit_m1_human_qc.R"))}
    A_out = {k: v for k, v in A.items() if k != "_tab"}
    res = {"SIMULATED": args.simulated, "label": "POST HOC / EXPLORATORY human visual-QC validation", "meta": meta,
           "input_issues": issues, "agreement": A_out, "sensitivity": B,
           "m1_refit": fit["scenarios"] if fit else None}
    with open(os.path.join(args.out, f"{label}human_qc_results.json"), "w") as fo:
        json.dump(clean(res), fo, indent=1, allow_nan=False)
    report(A, B, fit, meta, os.path.join(args.out, f"{label}human_qc_report.md"))
    print(json.dumps({"simulated": args.simulated, "n_rated": A["n_rated_primary_items"],
                      "kappa_linear_sample": A["sample_unweighted"]["kappa_linear"],
                      "kappa_linear_design": A["design_weighted_to_180"]["kappa_linear"],
                      "intra_exact": A["intra_rater_duplicates"]["exact_agreement"],
                      "refit_ok": B.get("refit_reproduces_frozen"),
                      "H_replace_beta": (fit or {}).get("scenarios", {}).get("H_replace", {}).get("beta_PS_ppb")}, indent=1))


if __name__ == "__main__":
    main()
