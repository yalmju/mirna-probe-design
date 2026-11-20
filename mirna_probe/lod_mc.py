"""Monte Carlo estimate of the detection limit for a probe assay read out
particle by particle (each particle scored ON or OFF).

One draw:
  1. capture kinetics of the target on 200 nm particles (diffusion + surface binding)
  2. per-particle signal (chance of a bright site, particle-to-particle brightness, noise)
  3. ON/OFF threshold from blank particles at a false-positive rate alpha
  4. Poisson LOD = (3.29*sqrt(lambda) + 3) / slope, best over the number of particles used

Uncertain parameters are drawn log-uniformly (or uniformly) within the ranges below.
With n_groups > 1 the read positions are shared by that many particle groups.
All values are assumptions for planning, not measurements.
"""
import numpy as np

NA = 6.022e23
V = 1e-5            # sample volume (L) = 10 uL
T_REACT = 7200.0    # reaction time (s) = 2 h
K_DIFF_200 = 1.13e11  # diffusion-limited collision rate for 200 nm particles (M^-1 s^-1)
N_OFF, N_ON = 20000, 4000
ALPHAS = (1e-2, 1e-3)
NP_GRID = np.logspace(4, 10, 25)   # particles in the sample


def _lu(rng, a, b, n=None):
    return 10 ** rng.uniform(np.log10(a), np.log10(b), n)


def _thresholds(off, on):
    out = []
    for a in ALPHAS:
        thr = np.quantile(off, 1 - a)
        out.append((a, (on > thr).mean()))
    return out


def _hot_signal(rng, n, fhot, ehot, g):
    nh = rng.binomial(n, fhot)
    return g * (nh * ehot + (n - nh))


def _lod(slope, lam):
    return (3.29 * np.sqrt(lam) + 3) / np.maximum(slope, 1e-300)


def simulate_lod(n_groups=1, n_draws=600, seed=7, labels_per_target=1):
    """Return an array of LOD values (fM), one per Monte Carlo draw (inf = no LOD).

    labels_per_target: labels one captured target leaves on its particle
    (1 for a sandwich; larger for an amplified readout, drawn as Poisson(m x retention)).
    """
    rng = np.random.default_rng(seed)
    lods = []
    for _ in range(n_draws):
        M = _lu(rng, 1e4, 1e5)                 # read positions
        phi = rng.uniform(0.5, 0.9)            # fraction of positions holding a particle
        fhot = _lu(rng, 0.005, 0.05)           # chance a reporter sits at a bright site
        ehot = _lu(rng, 30, 300)               # brightness gain at that site
        st = rng.uniform(0.2, 0.6)             # lognormal sigma of particle-to-particle brightness
        sn = _lu(rng, 1, 10)                   # readout noise
        G = lambda n: rng.lognormal(0, st, n)
        nh = int(_lu(rng, 1e3, 1e4))           # capture probes per particle
        kt = _lu(rng, 1e5, 1e6)                # surface on-rate (M^-1 s^-1)
        q = _lu(rng, 0.3, 20)                  # nonspecific reporters per particle
        cb = rng.uniform(0.1, 0.5)             # reporter brightness scale
        b0 = _lu(rng, 0.5, 5)                  # particle background
        theta = _lu(rng, 0.05, 1)              # reporter retained after washing
        noff = rng.poisson(q, N_OFF)
        base = rng.poisson(q, N_ON)   # draw order kept so seed 7 reproduces the reported values
        extra = rng.binomial(1, theta, N_ON) if labels_per_target == 1 else rng.poisson(labels_per_target * theta, N_ON)
        non = base + extra
        off = cb * _hot_signal(rng, noff, fhot, ehot, 1) * G(N_OFF)
        gg = G(N_ON)
        on = cb * _hot_signal(rng, non, fhot, ehot, 1) * gg
        off = off + b0 * G(N_OFF) + rng.normal(0, sn, N_OFF)
        on = on + b0 * gg + rng.normal(0, sn, N_ON)
        tb = _thresholds(off, on)
        best = np.inf
        for npart in NP_GRID:
            P = npart / (V * NA)
            k = 1 / (1 / K_DIFF_200 + 1 / (nh * kt))
            fcap = 1 - np.exp(-k * P * T_REACT)
            m_group = M * phi / n_groups
            r = min(1, m_group / npart)
            nread = min(npart, m_group)
            for a, beta in tb:
                if beta > a:
                    best = min(best, _lod(NA * V * fcap * r * (beta - a), nread * a))
        lods.append(best * 1e15)
    return np.array(lods)
