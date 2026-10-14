"""
clone_check.py

First look at clonality: group IGHV4-34 sequences by identical CDR3 amino-acid
sequence within each subject and isotype. Crude proxy for clones (ignores J gene
and CDR3 variants). Reads results/ighv4_34_sequences.csv.gz.
"""

import os

import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.abspath(os.path.join(script_dir, ".."))
results_dir = os.path.join(project_dir, "results")

df = pd.read_csv(os.path.join(results_dir, "ighv4_34_sequences.csv.gz"))
df = df[df["avy_status"].isin(["intact", "mutated"])]
df = df.dropna(subset=["cdr3_aa"]).copy()
df["intact"] = (df["avy_status"] == "intact").astype(int)

rows = []
for (subject, isotype), g in df.groupby(["subject", "isotype"]):
    clones = g.groupby("cdr3_aa").agg(n=("intact", "size"),
                                      frac_intact=("intact", "mean"))
    top3 = clones["n"].nlargest(3)
    rows.append({
        "subject": subject,
        "isotype": isotype,
        "sequences": len(g),
        "distinct_cdr3": len(clones),
        "top_clone_share_pct": round(top3.iloc[0] / len(g) * 100, 1),
        "top3_share_pct": round(top3.sum() / len(g) * 100, 1),
        "pct_intact_sequences": round(g["intact"].mean() * 100, 1),
        "pct_intact_clone_weighted": round(clones["frac_intact"].mean() * 100, 1),
    })

out = pd.DataFrame(rows)
print(out.to_string(index=False))
out.to_csv(os.path.join(results_dir, "clone_check.csv"), index=False)