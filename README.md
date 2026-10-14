# ighv4-34-avy-nhs-sle

**Do the germline self-reactivity motifs AVY and NHS in IGHV4-34 heavy-chain sequences get retained at different rates in IgM vs IgG repertoires of SLE patients?**

A beginner-level, exploratory antibody-repertoire analysis of public data from the Observed Antibody Space (OAS). It is a multi-subject follow-up to an earlier single-patient analysis.

**Status:** AVY analysis finished for 7 subjects (IgM vs IgG). NHS, mutation load and a paired figure are planned (see [Next steps](#next-steps)). Results are exploratory.

---

## Background

IGHV4-34 is a heavy-chain V gene whose germline sequence already encodes self-reactivity. Antibodies using it are over-represented in SLE, and two germline-encoded features are linked to that self-binding:

- **AVY** (Ala-Val-Tyr) at the end of framework region 1 (FR1)
- **NHS** (Asn-His-Ser) in CDR2

If the motifs are lost (by somatic hypermutation, or by selection against self-reactive cells), the self-reactivity is expected to drop. This project asks how often AVY is still intact in IgM and in class-switched IgG sequences of SLE patients. IgM vs IgG is used as a rough proxy for B-cell maturation.

## Data

- **Source:** Observed Antibody Space (OAS), unpaired heavy-chain data
- **Study:** Tipton et al., 2015
- **Subjects:** 7 SLE patients (Subject-SLE1 to Subject-SLE7)
- **Units used:** IGHM and IGHG for each subject (14 files)
- **Cells:** PBMC, unsorted B cells (OAS labels)
- **Raw data are not stored in this repo.** Download them from OAS (https://opig.stats.ox.ac.uk/webapps/oas/), using the filters Human, Heavy, Disease = SLE.

**Subject-SLE3 has two runs in OAS.** The larger run was used for both isotypes. This was decided before the analysis, and the smaller run was excluded, so SLE3 counts as one subject.

**There are no healthy donors from this study in OAS**, so this is a within-SLE analysis and no rate can be called high or low relative to normal.

### Expected folder layout

```
Data/
  SLE1/
    SLE1_SRR1959703_Heavy_IGHM.csv
    SLE1_SRR1959703_Heavy_IGHG.csv
  SLE2/ ... SLE7/
```

Each subject folder must contain exactly one file per isotype, and each file name must contain `IGHM` or `IGHG`. The script stops with an error otherwise. A manifest of run accessions will be added to the repo.

## Repository contents

```
scripts/
  motifs.py              IGHV4-34 and AVY detection logic
  analyze_subjects.py    per-subject, per-isotype AVY summary
  clone_check.py         clonality check (identical CDR3 groups)
  bootstrap_stats.py     clone-aware bootstrap intervals and sign test
tests/
  test_motifs.py         unit tests for the motif logic
summary_tables/          small result tables (see Results)
requirements.txt
```

## Methods

1. **Filtering.** Productive sequences with an unambiguous IGHV4-34 V call (every call in an ambiguous set must be IGHV4-34).
2. **AVY status.** From the last three residues of the FR1 amino-acid sequence: `intact` if AVY, `mutated` otherwise, `unknown` if FR1 is missing or shorter than three residues.
3. **Read length.** The reads start partway into FR1, so FR1 length varies (the most common length is 9 residues). All lengths are analyzed. A sensitivity check repeats the percentage using only the most common length.
4. **Per-subject summary.** Percent intact among classifiable sequences, for IgM and IgG separately.
5. **Clonality check.** Sequences are grouped by identical CDR3 amino-acid sequence within each subject and isotype. This is a crude clone proxy: it ignores the J gene and CDR3 variants.
6. **Intervals.** Bootstrap, 2,000 resamples, seed 42. The resampling unit is the CDR3 group, so related sequences move together, instead of being treated as independent observations.
7. **Test.** Exact two-sided sign test on the per-subject IgG minus IgM differences.

The motif logic has unit tests (`pytest`).

## Results (exploratory)

Percent of IGHV4-34 sequences with AVY intact, per subject (all read lengths):

| Subject | IgM | IgG | IgG minus IgM | 95% interval of the difference |
|---|---|---|---|---|
| SLE1 | 82.9 | 63.4 | -19.5 | -21.8 to -17.1 |
| SLE2 | 88.7 | 38.3 | -50.4 | -55.0 to -44.5 |
| SLE3 | 88.6 | 70.2 | -18.4 | -21.0 to -15.8 |
| SLE4 | 86.9 | 79.0 | -7.9 | -8.9 to -6.9 |
| SLE5 | 84.4 | 76.4 | -8.1 | -8.9 to -7.3 |
| SLE6 | 78.6 | 62.0 | -16.7 | -22.3 to -11.2 |
| SLE7 | 85.1 | 72.1 | -13.0 | -16.0 to -9.4 |

- **IgG retention was lower than IgM in 7 of 7 subjects** (exact sign test, p = 0.0156, the smallest value possible with 7 pairs). The interval for each subject's difference excludes zero.
- **IgM retention is similar across subjects (78 to 89%), while IgG varies widely (38 to 79%).** The size of the drop differs a lot between patients: about 8 points in SLE4 and SLE5, and 50 points in SLE2.
- **The result does not depend on the read-length rule.** Using only reads of the most common FR1 length gives percentages 0.3 to 2.7 points higher in all 14 rows, with the same IgM > IgG ordering in every subject.
- **It does not depend on a few large clones.** In the CDR3 grouping, the largest group takes at most 13.3% of any file (SLE7 IgG), and SLE2's low IgG retention is spread across its repertoire. The clone-weighted drop agrees with the sequence-level drop within about 4 points per subject.

Full tables are in [`summary_tables/`](summary_tables/): `avy_summary.csv` (counts and percentages per subject and isotype), `bootstrap_summary.csv` (intervals and differences) and `clone_check.csv` (clonality summary).

### What the results do and do not show

They describe how AVY retention differs between IgM and IgG within these SLE patients. They do **not** show why AVY is lost. Selection against self-reactive cells and the chance accumulation of somatic mutations would both produce lower retention in IgG. Separating them needs the mutation-load analysis below.

## Limitations

- **No healthy controls from the same study.** Nothing can be called elevated, reduced or abnormal relative to normal.
- **Only 7 subjects,** so the sign test cannot go below p = 0.0156. Treat the result as exploratory.
- **IgM vs IgG is a maturation proxy.** The cells are unsorted PBMC B cells, so IgM includes naive and memory cells and may include some antibody-secreting cells.
- **Reads cover only the end of FR1,** so only the AVY end is assessed.
- **The clone definition is crude** (identical CDR3). The bootstrap intervals are likely somewhat optimistic, especially for files with mostly single-sequence groups.
- **One sample per subject,** with no clinical metadata (disease activity, treatment, nephritis status) in OAS.
- **Amino-acid level only.** Silent nucleotide changes are not detected.

## Next steps

- [ ] Mutation load vs AVY status (using the `v_identity` column) to separate selection from chance accumulation
- [ ] NHS in CDR2, after checking the exact motif definition and position in the literature
- [ ] Paired IgM-to-IgG figure per subject, with bootstrap intervals
- [ ] Data manifest with run accessions and a download script
- [ ] Optional exploratory comparison with healthy donors from other OAS studies, clearly labeled as cross-study

## Running the analysis

```
pip install -r requirements.txt
pytest
python scripts/analyze_subjects.py
python scripts/clone_check.py
python scripts/bootstrap_stats.py
```

Outputs are written to `results/`. The scripts expect the raw files in `Data/` as described above.

## Data source and citation

- Tipton, C. M., et al. (2015). Antibody-secreting cell data for SLE patients, accessed through OAS.
- Kovaltsuk, A., et al. (2018). Observed Antibody Space. *J. Immunol.*
- Olsen, T. H., Boyles, F. and Deane, C. M. Observed Antibody Space: a diverse database of cleaned, annotated and translated unpaired and paired antibody sequences. *Protein Science.*

Check the OAS website for the current citation wording, and the original paper for the full reference of the Tipton study.

## License

MIT (see `LICENSE`).

