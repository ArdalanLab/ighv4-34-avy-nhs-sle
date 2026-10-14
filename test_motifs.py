import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from motifs import check_avy, is_ighv4_34

# Toy FR1 strings for testing, not a claim about the real germline sequence.
INTACT = "QVQLQQWGAGLLKPSETLSLTCAVY"   # 25 aa, ends in AVY
MUTATED = "QVQLQQWGAGLLKPSETLSLTCAVH"  # 25 aa, last residue changed


def test_exact_gene():
    assert is_ighv4_34("IGHV4-34*01")


def test_two_alleles_of_same_gene():
    assert is_ighv4_34("IGHV4-34*01,IGHV4-34*02")


def test_ambiguous_call_rejected():
    assert not is_ighv4_34("IGHV4-34*01,IGHV4-39*01")


def test_other_gene_rejected():
    assert not is_ighv4_34("IGHV3-23*01")


def test_missing_v_call():
    assert not is_ighv4_34(None)
    assert not is_ighv4_34(float("nan"))
    assert not is_ighv4_34("")


def test_avy_intact():
    assert check_avy(INTACT, expected_length=25) == "intact"


def test_avy_mutated():
    assert check_avy(MUTATED, expected_length=25) == "mutated"


def test_truncated_is_unknown():
    assert check_avy("QVQLQQ", expected_length=25) == "unknown"


def test_missing_fr1_is_unknown():
    assert check_avy(None) == "unknown"
    assert check_avy(float("nan")) == "unknown"