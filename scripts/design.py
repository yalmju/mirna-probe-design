"""Design and check probes for any miRNA.

Examples
  python scripts/design.py                                      # miR-153 reference, sandwich
  python scripts/design.py --name miR-7 --target UGGAAGACUAGUGAUUUUGUUGU
  python scripts/design.py --design full --target <SEQ>         # one probe over the whole target
  python scripts/design.py --target <SEQ> --off miR-153=UUGCAUAGUCACAAAAGUGA
  python scripts/design.py --target <SEQ> --temp 4 --salt 0.05  # other conditions
  python scripts/design.py --lod                                # add the Monte Carlo LOD
"""
import argparse
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from mirna_probe import DESIGNS, MIR153, cross_check, duplex_dg, kd, mismatch_scan, self_fold, simulate_lod

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--target", default=MIR153, help="target sequence 5'->3' (RNA or DNA letters)")
ap.add_argument("--name", default="miR-153", help="label used in output file names")
ap.add_argument("--design", default="sandwich", choices=sorted(DESIGNS), help="probe layout")
ap.add_argument("--temp", type=float, default=25.0, help="temperature in C (default 25)")
ap.add_argument("--salt", type=float, default=0.14, help="Na+ in M (default 0.14)")
ap.add_argument("--off", action="append", default=[], metavar="NAME=SEQ", help="off-target to check (repeatable)")
ap.add_argument("--lod", action="store_true", help="run the Monte Carlo LOD vs number of groups")
ap.add_argument("--labels", type=int, default=1, help="labels per captured target for --lod (default 1)")
a = ap.parse_args()

tgt = a.target.upper().replace("U", "T")
out = pathlib.Path("results")
out.mkdir(exist_ok=True)
tag = a.name.replace(" ", "_")

probes = DESIGNS[a.design](tgt, a.temp, a.salt)
print(f"{a.name} ({len(tgt)} nt), {a.design} design, {a.temp:g} C, {a.salt * 1000:g} mM Na+")
for p in probes:
    print(f"  {p['name']:9s} 5'-[{p['mods']}]-{p['seq']}-3'  covers {p['start'] + 1}-{p['end']}  "
          f"dG {p['dG']}  Kd {kd(p['dG'], a.temp):.0e} M  self-fold {self_fold(p['seq'], a.temp, a.salt)[1]:.1f}")
for i in range(len(probes)):
    for j in range(i + 1, len(probes)):
        print(f"  {probes[i]['name']}+{probes[j]['name']}: "
              f"{duplex_dg(probes[i]['seq'], probes[j]['seq'], a.temp, a.salt):.1f} kcal/mol")

scan = mismatch_scan(tgt, probes, a.temp, a.salt)
with open(out / f"{tag}_mismatch_scan.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=scan[0].keys())
    w.writeheader()
    w.writerows(scan)
print("  mismatch ddG by position:", " ".join(str(r["ddG"]) for r in scan))

if a.off:
    targets = {a.name: tgt}
    for o in a.off:
        n, s = o.split("=", 1)
        targets[n] = s
    m = cross_check({p["name"]: p["seq"] for p in probes}, targets, a.temp, a.salt)
    print("  cross-check dG (kcal/mol):")
    for pn, row in m.items():
        print(f"    {pn:9s}" + "  ".join(f"{tn} {v}" for tn, v in row.items()))

if a.lod:
    rows = []
    for n in (1, 3, 4, 6, 8):
        l = simulate_lod(n_groups=n, labels_per_target=a.labels)
        rows.append(dict(groups=n, median_fM=round(float(np.median(l))), p_le_100fM=round(float((l <= 100).mean()), 2)))
        print(f"  LOD, {n} group(s): median {rows[-1]['median_fM']} fM")
    with open(out / "lod_vs_groups.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
print(f"Saved results/{tag}_mismatch_scan.csv" + (" and results/lod_vs_groups.csv" if a.lod else ""))
