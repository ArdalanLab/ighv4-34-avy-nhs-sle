"""
analyze_subjects.py

For every subject folder in Data/ (SLE1 ... SLE7), measure AVY motif status in
productive IGHV4-34 sequences of the IGHM and IGHG files, and write:
  results/avy_summary.csv                 one row per subject and isotype
  results/ighv4_34_sequences.csv.gz       one row per IGHV4-34 sequence

IgM vs IgG is used as a maturation proxy (the cells are unsorted PBMC B cells).
The reads cover only the last part of FR1, so only the AVY end of FR1 is assessed.
"""

import math
import os
import sys

import pandas as pd

from motifs import check_avy, is_ighv4_34

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.abspath(os.path.join(script_dir, ".."))
data_dir = os.path.join(project_dir, "Data")
results_dir = os.path.join(project_dir, "results")

# Column names as they appear in the OAS files ("Redundancy" = identical copies
# collapsed into each unique sequence; "v_identity" = V-segment identity to germline)
COLUMNS = ["v_call", "productive", "fwr1_aa", "cdr2_aa", "cdr3_aa",
           "v_identity", "Redundancy"]
ISOTYPES = ["IGHM", "IGHG"]


def find_files():
    """Return sorted (subject, isotype, path) tuples; exit if one is ambiguous."""
    found = {}
    for subject in sorted(os.listdir(data_dir)):
        subject_path = os.path.join(data_dir, subject)
        if not (os.path.isdir(subject_path) and subject.upper().startswith("SLE")):
            continue
        for name in sorted(os.listdir(subject_path)):
            upper = name.upper()
            if not upper.endswith((".CSV", ".CSV.GZ")):
                continue
            matches = [iso for iso in ISOTYPES if iso in upper]
            if len(matches) != 1:
                print(f"WARNING: cannot assign an isotype to {name}; skipped")
                continue
            found.setdefault((subject, matches[0]), []).append(
                os.path.join(subject_path, name)
            )
    for key, paths in found.items():
        if len(paths) != 1:
            sys.exit(f"ERROR: {key} has {len(paths)} files: {paths}")
    return [(subject, iso, paths[0]) for (subject, iso), paths in sorted(found.items())]


def load_oas(path):
    """Read one OAS file (line 1 is metadata, line 2 is the column header)."""
    header = pd.read_csv(path, header=1, nrows=0).columns.tolist()
    missing = [c for c in COLUMNS if c not in header]
    if missing:
        sys.exit(
            f"ERROR: {os.path.basename(path)} lacks columns {missing}.\n"
            f"Available columns: {header}"
        )
    return pd.read_csv(path, header=1, usecols=COLUMNS, low_memory=False)


def wilson_ci(k, n, z=1.96):
    """95% Wilson interval for a proportion k/n (treats sequences as independent)."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return centre - half, centre + half


def percent(k, n):
    return round(k / n * 100, 1) if n else float("nan")


def analyze(path, subject, isotype):
    df = load_oas(path)
    n_total = len(df)

    productive = df["productive"].astype(str).str.upper().isin(["T", "TRUE"])
    df = df[productive]
    n_productive = len(df)

    df = df[df["v_call"].apply(is_ighv4_34)].copy()

    # The reads start partway into FR1, so FR1 length varies. We do NOT filter on
    # it; we only record it (and use the most common length for a sensitivity check).
    lengths = df["fwr1_aa"].dropna().str.len()
    mode_length = int(lengths.mode().iloc[0]) if len(lengths) else None
    print(f"{subject} {isotype}: FR1 length counts (top 5): "
          f"{lengths.value_counts().head(5).to_dict()}  -> most common {mode_length}")

    df["avy_status"] = df["fwr1_aa"].apply(check_avy)
    df["subject"] = subject
    df["isotype"] = isotype

    counts = df["avy_status"].value_counts()
    intact = int(counts.get("intact", 0))
    mutated = int(counts.get("mutated", 0))
    unknown = int(counts.get("unknown", 0))
    known = intact + mutated
    low, high = wilson_ci(intact, known)

    # Sensitivity check: only reads with the most common FR1 length
    same_len = df[(df["avy_status"] != "unknown")
                  & (df["fwr1_aa"].str.len() == mode_length)]
    n_same_len = len(same_len)
    intact_same_len = int((same_len["avy_status"] == "intact").sum())

    row = {
        "subject": subject,
        "isotype": isotype,
        "file": os.path.basename(path),
        "total_sequences": n_total,
        "productive": n_productive,
        "ighv4_34": len(df),
        "intact": intact,
        "mutated": mutated,
        "unknown": unknown,
        "percent_intact": percent(intact, known),
        "ci95_low": round(low * 100, 1),
        "ci95_high": round(high * 100, 1),
        "fr1_most_common_length": mode_length,
        "n_most_common_length": n_same_len,
        "percent_intact_most_common_length": percent(intact_same_len, n_same_len),
    }
    return row, df


def main():
    files = find_files()
    if not files:
        sys.exit(f"No subject files found under {data_dir}")

    rows, frames = [], []
    for subject, isotype, path in files:
        row, df = analyze(path, subject, isotype)
        rows.append(row)
        frames.append(df[["subject", "isotype", "v_call", "fwr1_aa", "cdr2_aa",
                          "cdr3_aa", "v_identity", "Redundancy", "avy_status"]])
        print(f"  {subject} {isotype}: {row['total_sequences']} total | "
              f"{row['ighv4_34']} IGHV4-34 | intact {row['intact']} | "
              f"mutated {row['mutated']} | unknown {row['unknown']} | "
              f"{row['percent_intact']}% intact "
              f"(95% CI {row['ci95_low']}-{row['ci95_high']}) | "
              f"most-common-length only: {row['percent_intact_most_common_length']}% "
              f"(n={row['n_most_common_length']})\n")

    os.makedirs(results_dir, exist_ok=True)
    summary_path = os.path.join(results_dir, "avy_summary.csv")
    pd.DataFrame(rows).to_csv(summary_path, index=False)
    pd.concat(frames).to_csv(
        os.path.join(results_dir, "ighv4_34_sequences.csv.gz"), index=False
    )
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()