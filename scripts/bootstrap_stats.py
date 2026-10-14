"""
bootstrap_stats.py

Per-subject AVY retention with clone-aware confidence intervals.

Resampling unit = distinct CDR3 amino-acid sequence within a subject and isotype
(a crude clone proxy). Each bootstrap draw resamples these groups with replacement,
so sequences in the same group move together instead of counting as independent.
Also reports the IgG minus IgM difference per subject and an exact sign test.
"""

import math
import os

import numpy as np
import pandas as pd

SEED = 42
N_BOOT = 2000

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.abspath(os.path.join(script_dir, ".."))
results_dir = os.path.join(project_dir, "results")

rng = np.random.default_rng(SEED)


def clone_counts(g):
    """Per CDR3 group: number of sequences (n) and number with intact AVY (k)."""
    t = g.groupby("cdr3_aa")["intact"].agg(n="size", k="sum")
    return t["n"].to_numpy(), t["k"].to_numpy()


def boot_percent(n, k):
    """Bootstrap the percent intact by resampling CDR3 groups with replacement."""
    m = len(n)
    out = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = rng.integers(0, m, m)
        out[b] = k[idx].sum() / n[idx].sum() * 100
    return out


def interval(x):
    low, high = np.percentile(x, [2.5, 97.5])
    return round(float(low), 1), round(float(high), 1)


df = pd.read_csv(os.path.join(results_dir, "ighv4_34_sequences.csv.gz"))
df = df[df["avy_status"].isin(["intact", "mutated"])]
df = df.dropna(subset=["cdr3_aa"]).copy()
df["intact"] = (df["avy_status"] == "intact").astype(int)

rows = []
for subject, sg in df.groupby("subject"):
    est, boots = {}, {}
    for iso in ["IGHM", "IGHG"]:
        n, k = clone_counts(sg[sg["isotype"] == iso])
        est[iso] = k.sum() / n.sum() * 100
        boots[iso] = boot_percent(n, k)

    diff = est["IGHG"] - est["IGHM"]
    diff_boot = boots["IGHG"] - boots["IGHM"]
    m_lo, m_hi = interval(boots["IGHM"])
    g_lo, g_hi = interval(boots["IGHG"])
    d_lo, d_hi = interval(diff_boot)
    rows.append({
        "subject": subject,
        "igm_percent_intact": round(est["IGHM"], 1),
        "igm_ci_low": m_lo, "igm_ci_high": m_hi,
        "igg_percent_intact": round(est["IGHG"], 1),
        "igg_ci_low": g_lo, "igg_ci_high": g_hi,
        "igg_minus_igm": round(diff, 1),
        "diff_ci_low": d_lo, "diff_ci_high": d_hi,
    })

out = pd.DataFrame(rows)
print(out.to_string(index=False))
out.to_csv(os.path.join(results_dir, "bootstrap_summary.csv"), index=False)

# Exact two-sided sign test on the per-subject IgG minus IgM differences
n_neg = int((out["igg_minus_igm"] < 0).sum())
n_pos = int((out["igg_minus_igm"] > 0).sum())
n = n_neg + n_pos
k = min(n_neg, n_pos)
p = min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2**n)
print(f"\nSign test: IgG lower than IgM in {n_neg} of {n} subjects, exact two-sided p = {p:.4f}")
