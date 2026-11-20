"""Hybridization energies for miRNA probes (any layout; see designs.py).

Energies use ViennaRNA with DNA parameters (Mathews 2004) as an approximation for
RNA:DNA hybrids, with ViennaRNA's salt correction. Use them for comparison, not as
absolute values.
"""
import math
import RNA

RNA.params_load_DNA_Mathews2004()

R_KCAL = 1.987e-3
MIR153 = "TTGCATAGTCACAAAAGTGA"   # miR-153 (first 20 nt, DNA letters; RNA = UUGCAUAGUCACAAAAGUGA)
SPACER = "T" * 10


def rc(seq):
    """Reverse complement (DNA letters)."""
    return seq.upper().replace("U", "T").translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _md(temp_c=25.0, na_m=0.14):
    md = RNA.md()
    md.temperature = temp_c
    md.salt = na_m
    return md


def duplex_dg(a, b, temp_c=25.0, na_m=0.14):
    """Minimum free energy (kcal/mol) of the a+b dimer."""
    fc = RNA.fold_compound(a + "&" + b, _md(temp_c, na_m))
    return fc.mfe_dimer()[1]


def self_fold(seq, temp_c=25.0, na_m=0.14):
    """(structure, dG) of a single strand."""
    return RNA.fold_compound(seq, _md(temp_c, na_m)).mfe()


def kd(dg, temp_c=25.0):
    """Dissociation constant (M) from dG (kcal/mol)."""
    return math.exp(dg / (R_KCAL * (temp_c + 273.15)))


_SUB = {"A": "C", "C": "A", "G": "T", "T": "G"}


def mismatch_scan(target, probes, temp_c=25.0, na_m=0.14):
    """Single-base substitution (A->C, C->A, G->T, T->G) at every target position.

    The probe whose region covers the position is scored: ddG = dG(mismatch) - dG(match).
    Each dG is rounded to 0.1 kcal/mol before subtracting, as in the reported tables.
    probes: list of dicts with name, seq, start, end (see designs.py).
    """
    target = target.upper().replace("U", "T")
    ref = {p["name"]: round(duplex_dg(p["seq"], target, temp_c, na_m), 1) for p in probes}
    out = []
    for i, b in enumerate(target):
        p = next((p for p in probes if p["start"] <= i < p["end"]), None)
        if p is None:
            out.append(dict(position=i + 1, base=b, probe="-", ddG=0.0))
            continue
        mut = target[:i] + _SUB[b] + target[i + 1:]
        ddg = round(duplex_dg(p["seq"], mut, temp_c, na_m), 1) - ref[p["name"]]
        out.append(dict(position=i + 1, base=b, probe=p["name"], ddG=round(ddg, 1)))
    return out


def cross_check(probes, targets, temp_c=25.0, na_m=0.14):
    """dG of every probe against every target. probes/targets: dict name -> sequence."""
    return {pn: {tn: round(duplex_dg(ps, ts.upper().replace("U", "T"), temp_c, na_m), 1)
                 for tn, ts in targets.items()} for pn, ps in probes.items()}
