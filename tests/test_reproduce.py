"""Check that the code reproduces the reported reference values."""
import numpy as np
from mirna_probe import MIR153, cross_check, full, mismatch_scan, sandwich, simulate_lod


def test_sandwich_design():
    cap, rep = sandwich(MIR153)
    assert cap["seq"] == "TTTTTTTTTTTCACTTTTGT" and cap["dG"] == -11.7
    assert rep["seq"] == "GACTATGCAA" and rep["dG"] == -12.9
    assert (rep["start"], rep["end"]) == (0, 10) and (cap["start"], cap["end"]) == (10, 20)


def test_full_design():
    p, = full(MIR153)
    assert p["seq"] == "TTTTTTTTTT" + "TCACTTTTGTGACTATGCAA"
    assert (p["start"], p["end"]) == (0, 20)


def test_mismatch_scan():
    ddg = [r["ddG"] for r in mismatch_scan(MIR153, sandwich(MIR153))]
    assert ddg == [0.8, 1.8, 4.5, 5.5, 4.1, 2.7, 3.4, 3.4, 1.6, 2.4,
                   1.0, 3.0, 4.5, 4.2, 4.2, 3.9, 5.1, 3.5, 3.3, 1.0]


def test_cross_check():
    probes = {p["name"]: p["seq"] for p in sandwich(MIR153)}
    m = cross_check(probes, {"miR-153": MIR153, "miR-125b": "UCCCUGAGACCCUAACUUGUGA"})
    assert m["capture"]["miR-153"] == -11.7 and m["capture"]["miR-125b"] > -6
    assert m["reporter"]["miR-153"] == -12.9 and m["reporter"]["miR-125b"] > -6


def test_lod_vs_groups():
    med = [round(float(np.median(simulate_lod(n_groups=n)))) for n in (1, 4, 8)]
    assert med == [139, 303, 467]


def test_more_labels_lower_lod():
    one = np.median(simulate_lod(n_draws=120, labels_per_target=1))
    many = np.median(simulate_lod(n_draws=120, labels_per_target=50))
    assert many < one
