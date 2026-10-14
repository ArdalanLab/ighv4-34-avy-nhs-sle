"""Motif detection for IGHV4-34 heavy-chain sequences (OAS format)."""


def is_ighv4_34(v_call):
    """True only if every (possibly ambiguous) V call is IGHV4-34."""
    if not isinstance(v_call, str) or not v_call.strip():
        return False
    genes = {call.strip().split("*")[0] for call in v_call.split(",")}
    return genes == {"IGHV4-34"}


def check_avy(fwr1, expected_length=None):
    """
    Return 'intact', 'mutated', or 'unknown' for one FR1 amino-acid sequence.

    'unknown' covers missing sequences and sequences whose length differs from
    expected_length (truncated reads or indels are not evidence of mutation).
    """
    if not isinstance(fwr1, str) or len(fwr1) < 3:
        return "unknown"
    if expected_length is not None and len(fwr1) != expected_length:
        return "unknown"
    return "intact" if fwr1[-3:] == "AVY" else "mutated"
