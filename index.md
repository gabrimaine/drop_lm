---
title: Tensions of a lipid-droplet lens from its shape
subtitle: Fitting the Deslandes–Thiam–Forêt membrane model to a droplet profile
short_title: Lipid-droplet lens fit
keywords: [lipid droplets, membrane elasticity, Helfrich, shape equations, surface tension]
---

+++ {"part": "abstract"}

A lipid droplet that nucleates between the two leaflets of the endoplasmic reticulum bilayer forms a lens. Its outline is set by the bending rigidity and tensions of the two monolayers that wrap the oil and by the tension of the bilayer it sits in. The model of Deslandes, Thiam and Forêt gives the free energy of this lens but no closed formula for its shape. We solve the model's Euler–Lagrange shape equations by shooting from the symmetry axis. We then fit the solutions to a measured droplet profile by least squares and read off the droplet volume, the bilayer-to-monolayer tension ratio, and the two tensions at an assumed monolayer bending rigidity. The fit describes the profile across its whole radial range. The droplet is large enough that the two tensions are determined separately, not only their combination $2\gamma_m-\gamma_b$. The tension ratio and the volume change by about 0.1 % or less under every analysis choice we varied. The absolute tensions scale with the assumed rigidity. Every number in this report is read from the analysis outputs. The tabs in [](#sec-universes) switch between alternative analysis choices.

+++

## Introduction

Lipid droplets (LDs) store neutral lipids. They are thought to form when oil accumulates between the two leaflets of the endoplasmic-reticulum (ER) bilayer, which then bulge into a lens [@thiam2016]. @deslandes2017 treat the two monolayers around the oil and the surrounding bilayer as elastic fluid membranes joined at a contact line. They find that, beyond a critical volume, the symmetric lens abruptly turns into a spherical protrusion. No curvature-inducing machinery is needed: the budding is driven by the droplet's large interfacial tension compared with the bilayer tension.

Here we use their model the other way round. Given a measured profile $z(r)$ of a droplet, we ask which volume and tensions produce it. This is not a curve fit in the usual sense. The theory specifies a free-energy functional, and $z(r)$ follows only by minimising it under the boundary conditions at the symmetry axis and at the contact line. Each evaluation of the model is therefore a boundary-value problem.

The paper's Supporting Material puts the monolayer bending rigidity at a few $k_BT$, following @marsh2006 {astra:cite}`prior_insights.monolayer_rigidity` and the interfacial tensions of a forming droplet around $10^{-3}$ N/m {astra:cite}`prior_insights.tension_scale`. We fix $\kappa_m$ = {astra:value col=value where="quantity=kappa_m" sig=1}`outputs.fit_summary` $k_BT$ and fit everything else.

## Model

### Free energy

Following @deslandes2017, the droplet is bounded by two monolayers, $c$ (cytosolic) and $l$ (luminal), which meet the bilayer $b$ at a circular contact line of radius $R$. The free energy is

$$
F = \sum_{i\in\{c,l,b\}} \left[ \int_{A_i} \left( \frac{\kappa_i}{2}\,(2H - C_{0i})^2 + \kappa^i_G K \right) dA + \gamma_i A_i \right] + 2\pi R\,\tau - P_l V_l - P_c V_c ,
$$ (eq-free-energy)

with bending rigidities $\kappa_i$, Gaussian moduli $\kappa^i_G$, spontaneous curvatures $C_{0i}$, tensions $\gamma_i$ and line tension $\tau$ [@helfrich1973]. $H$ and $K$ are the mean and Gaussian curvatures.

We use the paper's minimal symmetric model {astra:cite}`prior_insights.symmetric_reduction`: identical monolayers ($\kappa_c=\kappa_l\equiv\kappa_m$, $\gamma_c=\gamma_l\equiv\gamma_m$), equal pressures on both sides of the bilayer, and $\tau = C_0 = 0$. The bilayer then stays flat, the lens is up–down symmetric, and the monolayers meet the bilayer tangentially.

:::{astra} prior_insights.symmetric_reduction
:::

### Shape equations

In arc-length parametrisation $(r(s), z(s), \psi(s))$ the Euler–Lagrange equations of [](#eq-free-energy) are those of axisymmetric vesicles [@seifert1991; @julicher1994; @seifert1997]. For one monolayer they read (SI Eq. 18)

$$
\dot r = \cos\psi,\quad
\dot z = -\sin\psi,\quad
\dot\psi = u,\quad
\dot u = -\frac{u\cos\psi}{r} + \frac{\sin\psi\cos\psi}{r^2} - \frac{p\,r\cos\psi}{2\kappa_m} + \frac{t\,\sin\psi}{\kappa_m r},\quad
\dot t = \frac{\kappa_m}{2}\left(u^2 - \frac{\sin^2\psi}{r^2}\right) + \gamma_m - p\,r\sin\psi ,
$$ (eq-shape)

where $u$ is the meridional curvature, $p$ the Laplace pressure and $t$ a Lagrange-multiplier function. The equations conserve a Hamiltonian (SI Eq. 19), which vanishes on physical solutions.

The boundary conditions sit at both ends of the contour. At the axis $r=\psi=0$. At the contact line $\psi$ returns to $0$ (tangency with the flat bilayer), and the lateral forces balance, $2t_m(R) = \gamma_b R$. Evaluated with the Hamiltonian, the force balance becomes

$$
\kappa_m\,u(R)^2 = 2\gamma_m - \gamma_b \equiv \Delta\gamma .
$$ (eq-contact)

:::{dropdown} A factor ½ in the contact-line tension
The Supporting Material writes the contact-line tension as $t_i(R) = \gamma_i R - \kappa_i R\, u_i(R)^2$. Evaluating the Hamiltonian at $\psi = 0$ gives half the curvature term, $t_i(R) = \gamma_i R - \tfrac{1}{2}\kappa_i R\, u_i(R)^2$. Only the halved form agrees with the paper's own small-volume solution (SI Eqs. 30 and 33), which gives $\kappa_c u(R)^2 = \Delta\gamma$ exactly. The printed form would require $2\kappa_c u(R)^2=\Delta\gamma$. Our solver uses the halved form, i.e. [](#eq-contact).

```{astra} prior_insights.contact_tension_factor
```
:::

### Scales and limits

Lengths scale with the monolayer bending length $\lambda_m = \sqrt{\kappa_m/\gamma_m}$. In these units the symmetric shape depends on $V/\lambda_m^3$, $\gamma_b/\gamma_m$ and $\kappa_b/\kappa_m$ {astra:cite}`prior_insights.three_params`. With a flat bilayer $\kappa_b$ does not enter. A profile therefore fixes two dimensionless numbers and the length $\lambda_m$. Only an assumed $\kappa_m$ turns $\lambda_m$ into tensions, via $\gamma_m = \kappa_m/\lambda_m^2$ and $\gamma_b = (\gamma_b/\gamma_m)\,\gamma_m$.

Two limits bracket the problem:

- **Small droplets** ($V \ll \lambda_m^3$). The radius and thickness depend on the tensions only through $\Delta\gamma$ (SI Eqs. 33–34), so $\gamma_m$ and $\gamma_b$ cannot be told apart {astra:cite}`prior_insights.small_volume_degeneracy`.
- **Large droplets** ($V \gg \lambda_m^3$). The lens becomes two spherical caps with $\cos\psi_c = \gamma_b/2\gamma_m$, so the shape gives the ratio but no length scale {astra:cite}`prior_insights.large_volume_caps`.

A droplet between the two limits carries information on both.

## Methods

### Solving the shape equations

We integrate [](#eq-shape) in dimensionless form from a series expansion at the axis outward, with an 8th-order Runge–Kutta integrator [@virtanen2020]. The integration stops when $\psi$ first returns to zero. This is the shooting method of the Supporting Material {astra:cite}`prior_insights.shooting_method`. In the symmetric model only one unknown must be tuned for a given apex curvature: the pressure. We find it by root-bracketing [](#eq-contact), taking the single-valued lens branch.

{astra:ref}`outputs.solver_validation` checks the solver against the paper's analytic limits ({astra}`findings.solver_validated`).

:::{astra} outputs.solver_validation
:caption: Solver against the analytic limits of the Supporting Material, for $\gamma_b/\gamma_m = 0.1$ and $1$. Top row: contact radius and thickness against SI Eqs. 33–34, and the largest computed shapes against a Laplace cap. Bottom row: maximum slope and apparent cap angle against the large-volume angle of SI Eq. 40, the small-volume profile against SI Eq. 31, and the curvature reaching the contact condition [](#eq-contact).
:::

### Fitting

The profile ({astra}`inputs.profile`) gives the height of the oil/monolayer interface above the bilayer plane. We treat it as the upper half of a symmetric lens. The model has three parameters: $\ln\lambda_m$, $\gamma_b/\gamma_m$ and the apex curvature. An optional fourth, a vertical offset of the bilayer plane, is controlled by {astra}`decisions.baseline`. Beyond the contact line the model is the flat bilayer. Points just outside the contact line therefore constrain $R$.

We start from a coarse grid over shape parameters and then minimise the squared residuals in $z$ with a trust-region least-squares solver [@virtanen2020]. Volume, tensions, contact radius and thickness are derived from the best-fit shape. Their uncertainties come from the method set by {astra}`decisions.uncertainty`.

### Analysis choices

Four choices are declared in the analysis. Each is shown below with its rationale and options; the option used for the numbers in this report is marked. [](#sec-universes) shows the result of changing each one.

:::{astra} decisions.temperature
:::

:::{astra} decisions.baseline
:::

:::{astra} decisions.fit_range
:::

:::{astra} decisions.uncertainty
:::

## Results

### The fit

{astra:ref}`outputs.fit_plot` compares the best-fit lens with the data. The model follows the profile from the apex through the contact line to the flat bilayer. The rms residual is {astra:value col=value where="quantity=rms" sig=2}`outputs.fit_summary` pm over {astra:value col=value where="quantity=n_points" sig=3}`outputs.fit_summary` points, against a lens thickness of {astra:value col=value where="quantity=h" err=value_std sig=4}`outputs.fit_summary` nm. The largest residual is {astra:value col=value where="quantity=max_resid" sig=2}`outputs.fit_summary` pm. The residuals are not white: they oscillate slowly along $r$ and dip just outside the contact line, where the model's curvature jumps.

:::{astra} outputs.fit_plot
:caption: Measured profile (black) and its mirror image (grey) with the best-fit monolayer (red) and flat bilayer (blue). The dotted line marks the contact radius. Lower panel: residuals in picometres.
:::

:::{astra} findings.model_reproduces_profile
:::

### Volume and tensions

The fitted droplet volume is $V$ = {astra:value col=value where="quantity=V" err=value_std sig=5}`outputs.fit_summary` nm³, and the tension ratio is $\gamma_b/\gamma_m$ = {astra:value col=value where="quantity=gamma_ratio" err=value_std sig=4}`outputs.fit_summary`. With $\kappa_m$ = {astra:value col=value where="quantity=kappa_m" sig=1}`outputs.fit_summary` $k_BT$ at {astra:value}`decisions.temperature`, the tensions are

- monolayer: $\gamma_m$ = {astra:value col=value where="quantity=gamma_m" err=value_std sig=3}`outputs.fit_summary` mN/m,
- bilayer: $\gamma_b$ = {astra:value col=value where="quantity=gamma_b" err=value_std sig=3}`outputs.fit_summary` mN/m,
- spreading tension: $\Delta\gamma = 2\gamma_m-\gamma_b$ = {astra:value col=value where="quantity=delta_gamma" err=value_std sig=3}`outputs.fit_summary` mN/m.

These tensions are of the order the paper anticipates {astra:cite}`prior_insights.tension_scale`. Their bending length $\lambda_m$ = {astra:value col=value where="quantity=lambda_m" err=value_std sig=3}`outputs.fit_summary` nm. The contact radius is $R$ = {astra:value col=value where="quantity=R" err=value_std sig=4}`outputs.fit_summary` nm.

The reduced volume is $V/\lambda_m^3$ = {astra:value col=value where="quantity=V_over_lambda3" sig=3}`outputs.fit_summary`. This sits between the two limits of the model. The large-volume contact angle implied by the ratio, $\arccos(\gamma_b/2\gamma_m)$ = {astra:value col=value where="quantity=psi_large_volume" sig=3}`outputs.fit_summary`°, is larger than the steepest slope of the fitted monolayer, {astra:value col=value where="quantity=psi_max" sig=3}`outputs.fit_summary`°. Bending still flattens the droplet noticeably.

:::{astra} outputs.fit_summary
:caption: All fitted and derived quantities. Uncertainties are 1σ and use the method selected by the uncertainty decision.
:::

:::{dropdown} How the uncertainties should be read
Both uncertainty methods treat the residuals as independent. The residuals in {astra:ref}`outputs.fit_plot` are visibly correlated along $r$, so the quoted errors are lower bounds. The *z0-free* tab in [](#sec-universes) gives a better sense of the systematic error on the tensions: it shifts them by roughly one quoted σ. None of the choices moves the ratio or the volume appreciably.
:::

### Are the two tensions separately determined?

At a small reduced volume only $\Delta\gamma$ would be measurable {astra:cite}`prior_insights.small_volume_degeneracy`. {astra:ref}`outputs.identifiability` maps the profile $\Delta\chi^2$ over $(\gamma_m, \gamma_b/\gamma_m)$, re-optimising the apex curvature at each node. The minimum is closed and does not run along the dash-dotted lines of constant $\Delta\gamma$. Instead, the valley stretches along $\gamma_m$ at almost fixed ratio. The profile pins the ratio more tightly than either tension, but it does constrain both.

:::{astra} outputs.identifiability
:caption: Profile $\Delta\chi^2$ around the best fit (star), with residual scale taken from the best fit. Dashed: $\Delta\chi^2$ contours. Dash-dotted: lines of constant $2\gamma_m-\gamma_b$, the degeneracy direction of the small-volume limit.
:::

:::{astra} findings.tensions_separable
:::

(sec-universes)=
## Changing the analysis choices

Each tab below is a *universe*: the analysis rerun with exactly one choice changed from the baseline. Every figure is that universe's own output, regenerated from the same specification, and its legend shows the fitted values. The numbers in the text above are those of the baseline universe.

::::{tab-set}

:::{tab-item} Baseline
All defaults: {astra}`decisions.temperature.t310`, {astra}`decisions.baseline.fixed`, {astra}`decisions.fit_range.all`, {astra}`decisions.uncertainty.jacobian`.

![Baseline fit](results/baseline/fit_plot.png)

Values: [results/baseline/fit_summary.csv](results/baseline/fit_summary.csv)
:::

:::{tab-item} 298 K
{astra}`decisions.temperature.t298` instead of 310 K for converting $\kappa_m$ to joules. Lengths, volume and ratio are unchanged. Both tensions scale by the ratio of temperatures, since $\gamma_m = \kappa_m/\lambda_m^2$.

![Fit with kappa_m at 298 K](results/temp-298/fit_plot.png)

Values: [results/temp-298/fit_summary.csv](results/temp-298/fit_summary.csv)
:::

:::{tab-item} Free baseline
{astra}`decisions.baseline.free`. The bilayer plane gets a fitted offset. This is the largest shift of the tensions among the choices, about one quoted σ.

![Fit with free z0](results/z0-free/fit_plot.png)

![Identifiability with free z0](results/z0-free/identifiability.png)

Values: [results/z0-free/fit_summary.csv](results/z0-free/fit_summary.csv)
:::

:::{tab-item} r ≤ 40 nm
{astra}`decisions.fit_range.r_le_40`. The drifting far tail is dropped. With the bilayer plane fixed, the model is $z=0$ at those points whatever the parameters, so they cannot move the fit. The best fit is identical, and only the rms and the point count change.

![Fit restricted to r <= 40 nm](results/r-le-40/fit_plot.png)

Values: [results/r-le-40/fit_summary.csv](results/r-le-40/fit_summary.csv)
:::

:::{tab-item} Bootstrap
{astra}`decisions.uncertainty.residual_bootstrap`. Same best fit. The errors come from refitting resampled residuals and agree with the Jacobian estimate.

![Fit with bootstrap errors](results/bootstrap/fit_plot.png)

Values: [results/bootstrap/fit_summary.csv](results/bootstrap/fit_summary.csv)
:::

::::

:::{astra} findings.robust_to_choices
:::

:::{tip} Trying other combinations
The report renders the first universe in `universes/` (here `baseline`). To explore any combination of choices, add a file to `universes/` and run `lc materialize local`. `lc view` then shows every universe side by side, with each output's provenance. See `README.md` for step-by-step instructions.
:::

## Discussion

The symmetric lens model of @deslandes2017 accounts for the measured profile with three parameters and gives a consistent set of droplet parameters. The tension ratio and volume come straight from the shape. They do not depend on the bending rigidity, the temperature, or any analysis choice we varied. The absolute tensions inherit the assumed $\kappa_m$ linearly. Since monolayer rigidities are known only within a factor of a few {astra:cite}`prior_insights.monolayer_rigidity`, that assumption, rather than the fit, dominates their uncertainty. A measured $\kappa_m$ for the lipids in question would remove it.

The main modelling assumptions are up–down symmetry, a flat bilayer, and no line tension or spontaneous curvature. The paper argues these are sub-dominant for budding. A profile of the lower monolayer, or of the bilayer away from the droplet, would test the symmetry directly. The structured residuals near the contact line point to physics the continuum model leaves out at molecular scales, or to smoothing in how the profile was extracted.

## Data and code availability

The analysis is specified in `astra.yaml` and runs with [Lightcone](https://github.com/LightconeResearch). Inputs, code and every output of every universe live in the repository, with large files in git-annex; `README.md` explains how to get and rerun them. All inputs:

:::{astra} inputs
:::
