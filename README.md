# mirna-probe-design

Design and in-silico checks for **miRNA detection probes** on particles: hybridization energies
with ViennaRNA, a single-base mismatch scan, an off-target cross-check, and a Monte Carlo estimate
of the detection limit for a readout that scores each particle ON or OFF.
The reference target is **miR-153** (first 20 nt, `UUGCAUAGUCACAAAAGUGA`, miRBase).

## Probe layouts (`mirna_probe/designs.py`)

| Design | Probes |
|---|---|
| `sandwich` | a capture strand on the particle binds the target 3' half, a labelled reporter binds the 5' half |
| `full` | one capture strand complementary to the whole target |

```
sandwich, miR-153:
capture  5'-[Thiol-C6]-TTTTTTTTTT-TCACTTTTGT-3'   on the particle, pairs with the 3' half
target   3'-AGUGAAAACA·CUGAUACGUU-5'              miR-153
reporter 5'-[label]-GACTATGCAA-3'                 labelled reporter, pairs with the 5' half
```

Adding a layout is one function that returns probes with the target region each one pairs with;
the mismatch scan, cross-check and Monte Carlo then work on it unchanged.

## What it computes

| Module | Output |
|---|---|
| `mirna_probe/thermo.py` | duplex dG and Kd, self-folding, single-mismatch ddG at every target position, probe x target cross-check |
| `mirna_probe/designs.py` | probe sequences for a layout; `sandwich` picks the split so both halves bind with similar strength |
| `mirna_probe/lod_mc.py` | Monte Carlo distribution of the detection limit, and how it changes with the number of particle groups sharing a fixed number of read positions |

### Thermodynamics
ViennaRNA 2.x, DNA parameters (Mathews 2004) with ViennaRNA's salt correction, 25 °C, 140 mM Na⁺ by default.
Mismatch scan: each target base is substituted (A→C, C→A, G→T, T→G) and ddG = dG(mismatch) − dG(match)
for the probe covering that position.

### Detection limit (Monte Carlo)
Each draw samples uncertain parameters (particles read 10⁴–10⁵, probes per particle 10³–10⁴,
surface on-rate 10⁵–10⁶ M⁻¹s⁻¹, signal heterogeneity, particle-to-particle brightness lognormal σ 0.2–0.6,
nonspecific labels 0.3–20 per particle, label retention after washing 5–100 %), simulates blank and
positive particles, sets the ON threshold at a 1 % or 0.1 % false-positive rate and computes
LOD = (3.29√λ + 3)/slope. Sample 10 µL, 2 h. `labels_per_target` is 1 for a sandwich and larger for an
amplified readout. With N particle groups the read positions are divided by N.

## Run

```bash
pip install -r requirements.txt

# any miRNA (RNA or DNA letters, 5'->3')
python scripts/design.py --name miR-7 --target UGGAAGACUAGUGAUUUUGUUGU
python scripts/design.py --design full --target <SEQ>
python scripts/design.py --target <SEQ> --off other=<SEQ> --off other2=<SEQ>
python scripts/design.py --target <SEQ> --temp 4 --salt 0.05
python scripts/design.py --lod --labels 50          # amplified readout

python scripts/run_all.py   # reproduce the miR-153 reference values, writes results/*.csv
pytest -q                   # checks that the reference values are reproduced
```

Expected (seed 7): capture dG −11.7, reporter dG −12.9 kcal/mol; median LOD 139 / 256 / 303 / 391 / 467 fM
for 1 / 3 / 4 / 6 / 8 groups.

## Limits
- RNA:DNA hybrids are approximated with DNA parameters; surface effects and Mg²⁺ are not modelled.
  Use the energies to compare designs, not as absolute values.
- Monte Carlo ranges are planning assumptions (label brightness, particle collection and readout are
  assumed), so read trends rather than absolute values.
- The tests check reproducibility of the calculations, not agreement with experiment.

## License
MIT
