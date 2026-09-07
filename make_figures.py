"""Makes the plots in the README."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from kelp_model import simulate, ensemble


def single_run():
    rng = np.random.default_rng(4)
    P = simulate(12, rng=rng)

    plt.figure(figsize=(6, 3.6))
    plt.plot(P[:, 0], label="Canopy")
    plt.plot(P[:, 1], label="Below-surface")
    plt.xlabel("Time (months)")
    plt.ylabel("Biomass")
    plt.title("One run with random weather")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/single_realization.png", dpi=140)
    plt.close()


def ensemble_plot():
    paths, _ = ensemble(n_runs=300, n_months=12, seed=1)
    total = paths.sum(axis=2)
    months = np.arange(total.shape[1])
    lo, mid, hi = np.percentile(total, [5, 50, 95], axis=0)

    plt.figure(figsize=(6, 3.6))
    plt.fill_between(months, lo, hi, alpha=0.25, label="5th-95th percentile")
    plt.plot(months, mid, "r", lw=2, label="Median")
    plt.xlabel("Time (months)")
    plt.ylabel("Total biomass")
    plt.title("300 different weather sequences")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/ensemble.png", dpi=140)
    plt.close()


def decay_check():
    # calm case is deterministic so one run is enough
    P = simulate(10, event_prob=0.0)
    total = P.sum(axis=1)
    months = np.arange(len(total))

    plt.figure(figsize=(6, 3.6))
    plt.semilogy(months, total, "o-", label="Simulated")
    plt.semilogy(months, total[0] * 0.5997 ** months, "--", label="|lambda|^t")
    plt.xlabel("Time (months)")
    plt.ylabel("Total biomass (log)")
    plt.title("Calm decay vs. eigenvalue prediction")
    plt.legend()
    plt.grid(alpha=0.3, which="both")
    plt.tight_layout()
    plt.savefig("figures/decay_vs_theory.png", dpi=140)
    plt.close()


if __name__ == "__main__":
    single_run()
    ensemble_plot()
    decay_check()
    print("figures saved")
