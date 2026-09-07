"""
Kelp population model - MATH 142 project.

Two stages (canopy and below-surface), discrete monthly time steps:
    x_{t+1} = M @ x_t + w_t
"""

import numpy as np

# calm-condition dynamics, trace fixed at 1.175 by the project constraint
M = np.array([
    [0.543, 0.33],
    [-0.05, 0.632],
])

# weather coefficients - most biomass a fully severe event can move in a month
aC = 0.5    # storm hitting canopy
hC = 0.1    # heatwave hitting canopy
aB = 0.15   # storm helping below-surface (more light gets through)
hB = 0.4    # heatwave hurting below-surface

X0 = np.array([500.0, 100.0])


def weather(x, s_storm, s_heat):
    """Extra change from this month's weather, on top of M @ x."""
    C, B = x
    return np.array([
        -aC * s_storm * C - hC * s_heat * C,
         aB * s_storm * C - hB * s_heat * B,
    ])


def simulate(n_months=12, x0=X0, event_prob=1.0, rng=None):
    """
    One run of the model.

    event_prob is the chance a storm/heatwave happens in a given month.
    If it happens, severity is uniform on [0,1], otherwise 0.
    event_prob=1.0 gives the original version where something happens every month.
    """
    if rng is None:
        rng = np.random.default_rng()

    P = np.zeros((n_months + 1, 2))
    P[0] = x0

    for k in range(n_months):
        s_storm = rng.uniform(0, 1) if rng.random() < event_prob else 0.0
        s_heat = rng.uniform(0, 1) if rng.random() < event_prob else 0.0
        P[k + 1] = M @ P[k] + weather(P[k], s_storm, s_heat)

    return P


def eigenvalues():
    """Eigenvalues and |lambda|. |lambda| < 1 means the population shrinks."""
    vals = np.linalg.eigvals(M)
    return vals, np.max(np.abs(vals))


def collapse_month(P, frac=0.01):
    """First month total biomass drops below frac of where it started."""
    total = P.sum(axis=1)
    hit = np.flatnonzero(total < frac * total[0])
    return int(hit[0]) if len(hit) else None


def ensemble(n_runs=1000, n_months=60, event_prob=1.0, seed=None):
    """Run the model a bunch of times with different weather."""
    rng = np.random.default_rng(seed)
    paths = np.zeros((n_runs, n_months + 1, 2))
    times = []

    for i in range(n_runs):
        paths[i] = simulate(n_months, event_prob=event_prob, rng=rng)
        times.append(collapse_month(paths[i]))

    return paths, times


def frequency_sweep(probs=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), n_runs=500):
    """Does it matter how often weather events happen? (Not much, it turns out.)"""
    out = []
    for p in probs:
        _, times = ensemble(n_runs=n_runs, event_prob=p, seed=0)
        got = [t for t in times if t is not None]
        out.append((p, np.median(got) if got else None))
    return out


if __name__ == "__main__":
    vals, radius = eigenvalues()
    print(f"eigenvalues: {vals[0]:.4f}, {vals[1]:.4f}")
    print(f"|lambda| = {radius:.4f}")
    print("complex, so decay oscillates. |lambda| < 1, so it decays.\n")

    print("median collapse month vs. how often weather happens:")
    for p, median in frequency_sweep():
        print(f"  p={p:.1f}  ->  {median:.0f} months")
