"""Probe layouts. Each design returns a list of probes with the target region each one pairs with.

sandwich : capture on the particle (3' half of the target) + labelled reporter (5' half)
full     : one capture probe complementary to the whole target

A probe is a dict: name, seq (5'->3', bases only), mods (text), start/end (0-based target
region it pairs with, end exclusive) and dG (kcal/mol).
"""
from .thermo import SPACER, duplex_dg, rc


def sandwich(target, temp_c=25.0, na_m=0.14, goal=-12.5, min_half=8, spacer=SPACER):
    """Split the target so both halves bind with similar strength near `goal`."""
    target = target.upper().replace("U", "T")
    best = None
    for s in range(min_half, len(target) - min_half + 1):
        cap = spacer + rc(target[s:])
        rep = rc(target[:s])
        gc = duplex_dg(cap, target, temp_c, na_m)
        gr = duplex_dg(rep, target, temp_c, na_m)
        score = abs(gc - gr) + abs(gc - goal) + abs(gr - goal)
        if best is None or score < best[0]:
            best = (score, s, cap, rep, gc, gr)
    _, s, cap, rep, gc, gr = best
    return [dict(name="capture", seq=cap, mods="5' Thiol-C6", start=s, end=len(target), dG=round(gc, 1)),
            dict(name="reporter", seq=rep, mods="5' label", start=0, end=s, dG=round(gr, 1))]


def full(target, temp_c=25.0, na_m=0.14, spacer=SPACER):
    """One capture probe covering the whole target."""
    target = target.upper().replace("U", "T")
    cap = spacer + rc(target)
    return [dict(name="capture", seq=cap, mods="5' Thiol-C6", start=0, end=len(target),
                 dG=round(duplex_dg(cap, target, temp_c, na_m), 1))]


DESIGNS = {"sandwich": sandwich, "full": full}
