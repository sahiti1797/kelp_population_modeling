# Kelp Population Model

A two-stage discrete-time linear model of giant kelp (Macrocystis pyrifera) biomass under random storm and heatwave conditions. Done for MATH 142 at UCLA, then extended a bit to add the Monte Carlo part.

## Basics of the Model

The overall kelp population is split into canopy (at the surface) and below-surface growth. Storms rip up the
canopy but let more light through to the layer below, hence aiding the growth of the below surface growth, while heatwaves hurt both.

State is `x_t = [canopy, below-surface]`:

```
x_{t+1} = M x_t + w_t
```

`M` is what happens in calm months, `w_t` is the extra effect of that month's weather. Storm and heat severity are drawn from [0,1] each month.

```
M = [ 0.543   0.33  ]
    [-0.05    0.632 ]
```

Parameters come from published field studies rather than fitting. Frond half-lives were used for survival rates, canopy-removal experiments for the shading term, storm and heatwave studies for the weather coefficients. Full details and sources are in `report.pdf`.

## Results

Eigenvalues of M are 0.5875 ± 0.1205i, so |λ| = 0.60.

Two things follow. |λ| < 1 means the population shrinks even with no weather at all. The eigenvalues being complex means it doesn't shrink smoothly — biomass swings back and forth between canopy-heavy and below-surface-heavy because the two stages mature at different
speeds.

![decay](figures/decay_vs_theory.png)

Simulating calm conditions matches |λ|^t, which checks the eigenvalue math.

![ensemble](figures/ensemble.png)

Running 300 random weather sequences, every single one declines. Weather makes it faster or slower but never causes growth.

I also swept how often weather events happen to see if that changed time to collapse:

| chance of event per month | median collapse |
|---|---|
| 0.0 | 7 months |
| 0.4 | 6 months |
| 1.0 | 6 months |

Barely anything. At |λ| = 0.6 the population is already losing 40% per month on its own, so the weather term doesn't have much room to matter.

## Limitations

If this model were right, kelp would be extinct. That comes from assuming a closed population with no spores coming in from outside, which is what pushes |λ| below 1. Adding recruitment (that I did not consider) would fix it.

Also: the model is linear so there's no carrying capacity, weather is equally likely in every month which ignores seasons, and sea urchins aren't in it at all even though they are a big cause of real kelp die-offs.

## Running it

```bash
pip install -r requirements.txt
python kelp_model.py      # eigenvalues + the frequency sweep
python make_figures.py    # plots
```

## Files

- `kelp_model.py` — the model and analysis
- `make_figures.py` — plots
- `kelp_simulation.ipynb` — original notebook from the class
- `report.pdf` — full writeup with derivation and sources
