"""Reproduce the miR-153 reference values (sandwich design)."""
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from mirna_probe import MIR153, duplex_dg, kd, mismatch_scan, sandwich, self_fold, simulate_lod

out = pathlib.Path("results")
out.mkdir(exist_ok=True)

probes = sandwich(MIR153)
print("miR-153 sandwich design (25 C, 140 mM Na+)")
for p in probes:
    print(f"  {p['name']:9s} 5'-[{p['mods']}]-{p['seq']}-3'  dG {p['dG']}  Kd {kd(p['dG']):.0e} M  "
          f"self-fold {self_fold(p['seq'])[1]:.1f}")
print(f"  capture+reporter {duplex_dg(probes[0]['seq'], probes[1]['seq']):.1f} kcal/mol")

scan = mismatch_scan(MIR153, probes)
with open(out / "mismatch_scan.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=scan[0].keys())
    w.writeheader()
    w.writerows(scan)
print("  mismatch ddG by position:", " ".join(str(r["ddG"]) for r in scan))

rows = []
for n in (1, 3, 4, 6, 8):
    l = simulate_lod(n_groups=n)
    rows.append(dict(groups=n, median_fM=round(float(np.median(l))), p_le_100fM=round(float((l <= 100).mean()), 2)))
    print(f"  LOD, {n} group(s): median {rows[-1]['median_fM']} fM")
with open(out / "lod_vs_groups.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
print("Saved results/mismatch_scan.csv and results/lod_vs_groups.csv")
