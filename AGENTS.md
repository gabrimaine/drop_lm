# drop_lm — fitting lipid-droplet lens profiles with the Deslandes–Thiam–Forêt model

A Lightcone/ASTRA project. `astra.yaml` is the source of truth for inputs, outputs,
decisions; `universes/baseline.yaml` selects the defaults; outputs land in
`results/<universe>/<output>.<format>` and are made only by `lc materialize`.

## Commands

```bash
lc compute launch --wait                 # local cluster, named `local`
lc status                                 # what is current / stale
lc materialize local                      # (re)make everything stale
lc run local -- python src/fit.py ...     # probe (write to results/probe/, delete before committing)
uvx astra-tools@0.2.18 validate           # after any spec edit
```

## Code map

- `src/lens.py` — shape equations (SI Eq. 18) in units lambda = sqrt(kappa_m/gamma_m);
  shooting from the axis (series start psi = a s + (a - P/2)/8 s^3), contact line = psi
  returning to 0, contact condition kappa_m u(R)^2 = 2 gamma_m - gamma_b.
- `src/fitting.py` — least squares in theta = (ln lambda, g = gamma_b/gamma_m, ln a [, z0]).
- `src/fit.py`, `fit_profile.py`, `plot_fit.py`, `identifiability.py`, `validate_solver.py` — one script per output.

## Project Notes

- Request (2026-11-27): fit the curve `td_reset2.xvg` with the SI model of Deslandes,
  Thiam & Forêt (BJ 2017), extract V, gamma_b/gamma_m, gamma_m, gamma_b, with
  kappa_m = 6 kBT; tensions expected ~1 mN/m. User accepted all default decisions.
- Assumed (not yet confirmed by the user): z(r) is the oil/monolayer interface of the
  upper half of an up-down symmetric lens, in nm, measured from the flat bilayer plane;
  r is distance from the axis. Simulation temperature unknown -> default 310 K.
- Physics reductions used: symmetric model (SI Sect. 5) => bilayer flat, psi_contact = 0,
  kappa_b and Gaussian moduli drop out. SI text states t_i(R) = gamma_i R - kappa_i R u^2;
  Eq. (19) and the small-volume result R^2 = 8 sqrt(kappa Delta gamma)/p require the
  factor 1/2 (t = gamma R - kappa R u^2 / 2), which is what the code uses.
- Solver structure: for fixed apex curvature a, contact solutions exist only for
  P > 2a (P = 2a is the closed sphere); search on x = P/(2a) - 1 on a log grid.
  Two roots usually: the lens (psi_max < 90 deg, largest x) and an overhanging
  near-sphere; the code keeps the lens. Very large lenses (R/lambda >~ 20) need
  x < 1e-9 and become ill-conditioned for shooting.
- Small volumes (V << lambda^3) only constrain 2 gamma_m - gamma_b; this droplet has
  V/lambda^3 ~ 550, so both tensions are separable (see `identifiability`).
- Literature pass offered and skipped (defaults). No prior_insights recorded.
