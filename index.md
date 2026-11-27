---
title: Tensions of a lipid-droplet lens from its shape
subtitle: Fitting the Deslandes–Thiam–Forêt membrane model to a droplet profile
short_title: Lipid-droplet lens fit
keywords: [lipid droplets, membrane elasticity, Helfrich, shape equations, surface tension]
---

+++ {"part": "abstract"}

A lipid droplet that nucleates between the two leaflets of the endoplasmic reticulum bilayer forms a lens. Its outline is set by the bending rigidity and tensions of the two monolayers that wrap the oil and by the tension of the bilayer it sits in. The model of Deslandes, Thiam and Forêt gives the free energy of this lens but no closed formula for its shape. We solve the model's Euler–Lagrange shape equations by shooting from the symmetry axis. We then fit the solutions to a measured droplet profile by least squares and read off the droplet volume, the bilayer-to-monolayer tension ratio, and the two tensions at an assumed monolayer bending rigidity. The fit describes the profile across its whole radial range. The droplet is large enough that the two tensions are determined separately, not only their combination $2\gamma_m-\gamma_b$. The tension ratio and the volume change by about 0.1 % or less under every analysis choice we varied. The absolute tensions scale with the assumed rigidity. A separately written code reproduces the fit. Comparing the two exposed a small sampling bias in our first version, which is now fixed. Every number in this report is read from the analysis outputs. The tabs in [](#sec-universes) switch between alternative analysis choices.

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

(sec-solver)=
### Solving the shape equations

We integrate [](#eq-shape) in dimensionless form from a series expansion at the axis outward, with an 8th-order Runge–Kutta integrator [@virtanen2020]. The integration stops when $\psi$ first returns to zero. This is the shooting method of the Supporting Material {astra:cite}`prior_insights.shooting_method`. In the symmetric model only one unknown must be tuned for a given apex curvature: the pressure. We find it by root-bracketing [](#eq-contact), taking the single-valued lens branch. The model height at the data radii is read from the integrator's continuous solution, sampled at 4001 points along the contour, which keeps the interpolation error below $10^{-6}$ nm.

[](#output-solver_validation) checks the solver against the paper's analytic limits ({astra}`findings.solver_validated`).

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

[](#output-fit_plot) compares the best-fit lens with the data. The model follows the profile from the apex through the contact line to the flat bilayer. The rms residual is {astra:value col=value where="quantity=rms" sig=2}`outputs.fit_summary` pm over {astra:value col=value where="quantity=n_points" sig=3}`outputs.fit_summary` points, against a lens thickness of {astra:value col=value where="quantity=h" err=value_std sig=4}`outputs.fit_summary` nm. The largest residual is {astra:value col=value where="quantity=max_resid" sig=2}`outputs.fit_summary` pm. The residuals are not white: they oscillate slowly along $r$ and dip just outside the contact line, where the model's curvature jumps.

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
Both uncertainty methods treat the residuals as independent. The residuals in [](#output-fit_plot) are visibly correlated along $r$, so the quoted errors are lower bounds. The *z0-free* tab in [](#sec-universes) gives a sense of the systematic error on the tensions: it shifts them by about half a quoted σ. None of the choices moves the ratio or the volume appreciably.
:::

### Are the two tensions separately determined?

At a small reduced volume only $\Delta\gamma$ would be measurable {astra:cite}`prior_insights.small_volume_degeneracy`. [](#output-identifiability) maps the profile $\Delta\chi^2$ over $(\gamma_m, \gamma_b/\gamma_m)$, re-optimising the apex curvature at each node. The minimum is closed and does not run along the dash-dotted lines of constant $\Delta\gamma$. Instead, the valley stretches along $\gamma_m$ at almost fixed ratio. The profile pins the ratio more tightly than either tension, but it does constrain both.

:::{astra} outputs.identifiability
:caption: Profile $\Delta\chi^2$ around the best fit (star), with residual scale taken from the best fit. Dashed: $\Delta\chi^2$ contours. Dash-dotted: lines of constant $2\gamma_m-\gamma_b$, the degeneracy direction of the small-volume limit.
:::

:::{astra} findings.tensions_separable
:::

(sec-cross-check)=
## Comparison with an independent implementation

A fit like this one can go wrong in many places without anything looking wrong: a sign in the shape equations, a missed branch of solutions, a loose tolerance, an interpolation. Matching the paper's analytic limits ([](#output-solver_validation)) tests the solver but not the fit. The profile was therefore fitted a second time by a separately written code, {astra}`inputs.other_code` ([`data/raw/compare_fit.py`](#code-compare-fit)). Nothing from our fit goes into it. [`src/cross_check.py`](#code-cross-check) calls the script's own fitting routine from its own five starting guesses, in parallel, and tabulates its result next to ours. It skips the script's additional coarse grid search, because on this profile all of the script's starts, grid ones included, end at the same optimum.

### Two codes, one model

The codes share the model and the data, and little else. Both solve the shape equations [](#eq-shape) with the contact condition [](#eq-contact) and take the model to be the flat bilayer beyond the contact line. Both interpolate the model linearly onto the data radii and minimise the unweighted squared residuals at the same $\kappa_m$. [](#tbl-two-codes) lists where they differ.

:::{list-table} How the two codes solve the same problem.
:label: tbl-two-codes
:header-rows: 1
* -
  - This analysis ([`src/lens.py`](#code-lens), [`src/fitting.py`](#code-fitting))
  - Independent code ([`data/raw/compare_fit.py`](#code-compare-fit))
* - Unknown tuned when shooting
  - Pressure, at fixed apex curvature: 90-point log scan of $p/2a-1$, then Brent's method
  - Apex curvature, at fixed pressure: 48-point log scan, then bisection, including brackets where the contour stops being a lens
* - Fit parameters
  - $\ln\lambda_m$, $\gamma_b/\gamma_m$, $\ln(a\lambda_m)$
  - $\ln\lambda_m$, $\operatorname{logit}(\delta/2)$ with $\delta = 2-\gamma_b/\gamma_m$, $\ln(p\lambda_m^3/\kappa_m)$
* - Integrator
  - 8th-order Runge–Kutta (DOP853), relative tolerance $10^{-10}$
  - 5th-order Runge–Kutta (RK45), relative tolerance $10^{-10}$
* - Lens branch
  - Largest-pressure root whose slope stays below 90°
  - Integration stops if the slope reaches 90°; root nearest the previous one
* - Contour for the model height
  - 4001 points of the integrator's continuous solution
  - 2500 points of the integrator's continuous solution
* - Volume
  - Trapezoid rule over the sampled contour
  - Integrated along the contour as an extra equation
* - Starting point
  - Grid of 12 × 16 × 80 over ($\gamma_b/\gamma_m$, $a$, $\lambda_m$)
  - Five guesses from the small-volume formulas (SI Eqs. 33–34) at the measured $R$ and $h$
* - Errors
  - Jacobian covariance, propagated to every quantity with all correlations
  - Jacobian covariance; $V$, $R$, $h$ with correlations, $\gamma_b$ in quadrature
* - Bilayer-plane offset
  - Optional fitted $z_0$ ({astra}`decisions.baseline`)
  - None: with a free plane, our fitted $z_0$ is subtracted from the data first
:::

Two of these differences make the comparison a real test. The codes tune different unknowns when shooting, so a root-finding error in one would not repeat in the other. They also start from unrelated guesses, so landing on the same optimum shows that neither stopped in a local minimum.

### The same best fit

```{list-table} This analysis and the independent code on the same profile. The last column is the difference in units of our 1σ. The independent code does not report an error on $\Delta\gamma$; the one shown propagates its covariance with all correlations.
:label: tbl-cross-check
:header-rows: 1
* - Quantity
  - This analysis
  - Independent code
  - Difference (σ)
* - $V$ (nm³)
  - {astra:value col=ours err=ours_std where="quantity=V" sig=6}`outputs.cross_check`
  - {astra:value col=other err=other_std where="quantity=V" sig=6}`outputs.cross_check`
  - {astra:value col=diff_sigma where="quantity=V" sig=1}`outputs.cross_check`
* - $\gamma_b/\gamma_m$
  - {astra:value col=ours err=ours_std where="quantity=gamma_ratio" sig=5}`outputs.cross_check`
  - {astra:value col=other err=other_std where="quantity=gamma_ratio" sig=5}`outputs.cross_check`
  - {astra:value col=diff_sigma where="quantity=gamma_ratio" sig=1}`outputs.cross_check`
* - $\gamma_m$ (mN/m)
  - {astra:value col=ours err=ours_std where="quantity=gamma_m" sig=5}`outputs.cross_check`
  - {astra:value col=other err=other_std where="quantity=gamma_m" sig=5}`outputs.cross_check`
  - {astra:value col=diff_sigma where="quantity=gamma_m" sig=1}`outputs.cross_check`
* - $\gamma_b$ (mN/m)
  - {astra:value col=ours err=ours_std where="quantity=gamma_b" sig=5}`outputs.cross_check`
  - {astra:value col=other err=other_std where="quantity=gamma_b" sig=5}`outputs.cross_check`
  - {astra:value col=diff_sigma where="quantity=gamma_b" sig=1}`outputs.cross_check`
* - $\Delta\gamma = 2\gamma_m-\gamma_b$ (mN/m)
  - {astra:value col=ours err=ours_std where="quantity=delta_gamma" sig=5}`outputs.cross_check`
  - {astra:value col=other err=other_std_full where="quantity=delta_gamma" sig=5}`outputs.cross_check`
  - {astra:value col=diff_sigma where="quantity=delta_gamma" sig=1}`outputs.cross_check`
* - $\lambda_m$ (nm)
  - {astra:value col=ours err=ours_std where="quantity=lambda_m" sig=5}`outputs.cross_check`
  - {astra:value col=other err=other_std where="quantity=lambda_m" sig=5}`outputs.cross_check`
  - {astra:value col=diff_sigma where="quantity=lambda_m" sig=1}`outputs.cross_check`
* - $R$ (nm)
  - {astra:value col=ours err=ours_std where="quantity=R" sig=6}`outputs.cross_check`
  - {astra:value col=other err=other_std where="quantity=R" sig=6}`outputs.cross_check`
  - {astra:value col=diff_sigma where="quantity=R" sig=1}`outputs.cross_check`
* - $h$ (nm)
  - {astra:value col=ours err=ours_std where="quantity=h" sig=6}`outputs.cross_check`
  - {astra:value col=other err=other_std where="quantity=h" sig=6}`outputs.cross_check`
  - {astra:value col=diff_sigma where="quantity=h" sig=1}`outputs.cross_check`
* - reduced pressure $p\lambda_m^3/\kappa_m$
  - {astra:value col=ours where="quantity=pressure" sig=6}`outputs.cross_check`
  - {astra:value col=other where="quantity=pressure" sig=6}`outputs.cross_check`
  -
* - apex curvature $\times\,\lambda_m$
  - {astra:value col=ours where="quantity=apex_curvature" sig=6}`outputs.cross_check`
  - {astra:value col=other where="quantity=apex_curvature" sig=6}`outputs.cross_check`
  -
* - rms residual (pm)
  - {astra:value col=ours where="quantity=rms" sig=5}`outputs.cross_check`
  - {astra:value col=other where="quantity=rms" sig=5}`outputs.cross_check`
  -
```

[](#tbl-cross-check) and [](#output-cross_check_plot) show the result. The two models cannot be told apart on the data (panel a), and their residuals coincide point by point (panel b). Panel c measures the distance between the two best-fit models directly. It is at most {astra:value col=other where="quantity=model_diff_max" sig=2}`outputs.cross_check` pm, more than a thousand times below the residual rms. No quantity differs by more than {astra:value col=other where="quantity=max_abs_diff_sigma" sig=2}`outputs.cross_check` of our σ (panel d, blue). The shooting variables agree too, although each code fixes the one the other solves for: the reduced pressure and the apex curvature of the two codes match to within one part in $10^5$.

:::{astra} outputs.cross_check_plot
:caption: This analysis against the independent code. (a) Data with both best-fit models. (b) Residuals of both. (c) Distance between the two best-fit models (blue) and, for scale, the interpolation error of our old contour sampling on the same shape (orange) and the residual rms (dashed). (d) Differences from our best fit in units of our σ: the independent code (blue) and our own code before the sampling fix (orange). (e) The independent code's errors relative to ours, as it reports them (grey) and with all correlations (blue).
:::

### The error bars

Both codes estimate their errors in the same way, from the Jacobian of the residuals at the best fit with the residual variance as the noise level. They agree for every quantity except $\gamma_b$. There the independent code reports {astra:value col=other_std where="quantity=gamma_b" sig=2}`outputs.cross_check` mN/m against our {astra:value col=ours_std where="quantity=gamma_b" sig=2}`outputs.cross_check` mN/m. The fits are not at fault; the propagation is. The script forms $\gamma_b = (\gamma_b/\gamma_m)\,\gamma_m$ and adds the errors of the two factors in quadrature. That assumes the ratio and $\gamma_m$, i.e. $\lambda_m$, are uncorrelated, but the fit correlates them strongly (the correlation matrix is in {astra}`outputs.fit_params`). Propagating the independent code's own covariance with the correlation gives {astra:value col=other_std_full where="quantity=gamma_b" sig=3}`outputs.cross_check` mN/m, the same as ours (panel e). The same propagation gives the error on $\Delta\gamma$ that the script does not print, {astra:value col=other_std_full where="quantity=delta_gamma" sig=3}`outputs.cross_check` mN/m against our {astra:value col=ours_std where="quantity=delta_gamma" sig=3}`outputs.cross_check` mN/m.

:::{astra} findings.independent_reproduction
:::

### What the comparison caught

The first run of the comparison did not agree this well. The tensions, $\lambda_m$ and $R$ differed by 0.3–0.4σ. Nobody would have noticed that in a single fit, but two correct codes should agree far better. The cause was in our code. Its solver integrated the shape equations as accurately as it does now, but returned the contour only at the integrator's own steps, capped at {astra:value col=coarse where="quantity=max_step" sig=2}`outputs.sampling_bias` $\lambda_m$. The model height at the data radii was then interpolated linearly between those points. On this lens that interpolation is off by up to {astra:value col=coarse where="quantity=interp_max" sig=2}`outputs.sampling_bias` pm ({astra:value col=coarse where="quantity=interp_rms" sig=2}`outputs.sampling_bias` pm rms; panel c, orange). That is small against the 8 pm residuals, but it is not noise. It is fixed by the shape and grows towards the contact line, where the contour bends most and where the data fix $R$, so the optimiser can absorb it into the parameters.

The solver now reads the height off the integrator's continuous solution at 4001 points ([](#sec-solver)). To measure what the old sampling cost, {astra}`outputs.sampling_bias` reruns the fit exactly as it ran before the fix. [`src/old_sampling.py`](#code-old-sampling) keeps the old sampling function line for line. In the baseline universe the rerun reproduces the archived pre-fix results to the last digit. [](#tbl-sampling-bias) lists the shifts.

```{list-table} Our fit with the current contour sampling and with the sampling used before the fix. The last column is the shift in units of our current 1σ.
:label: tbl-sampling-bias
:header-rows: 1
* - Quantity
  - Current sampling
  - Old sampling (every {astra:value col=coarse where="quantity=max_step" sig=2}`outputs.sampling_bias` $\lambda_m$)
  - Shift (σ)
* - $V$ (nm³)
  - {astra:value col=dense err=dense_std where="quantity=V" sig=6}`outputs.sampling_bias`
  - {astra:value col=coarse err=coarse_std where="quantity=V" sig=6}`outputs.sampling_bias`
  - {astra:value col=shift_sigma where="quantity=V" sig=2}`outputs.sampling_bias`
* - $\gamma_b/\gamma_m$
  - {astra:value col=dense err=dense_std where="quantity=gamma_ratio" sig=5}`outputs.sampling_bias`
  - {astra:value col=coarse err=coarse_std where="quantity=gamma_ratio" sig=5}`outputs.sampling_bias`
  - {astra:value col=shift_sigma where="quantity=gamma_ratio" sig=2}`outputs.sampling_bias`
* - $\gamma_m$ (mN/m)
  - {astra:value col=dense err=dense_std where="quantity=gamma_m" sig=5}`outputs.sampling_bias`
  - {astra:value col=coarse err=coarse_std where="quantity=gamma_m" sig=5}`outputs.sampling_bias`
  - {astra:value col=shift_sigma where="quantity=gamma_m" sig=2}`outputs.sampling_bias`
* - $\gamma_b$ (mN/m)
  - {astra:value col=dense err=dense_std where="quantity=gamma_b" sig=5}`outputs.sampling_bias`
  - {astra:value col=coarse err=coarse_std where="quantity=gamma_b" sig=5}`outputs.sampling_bias`
  - {astra:value col=shift_sigma where="quantity=gamma_b" sig=2}`outputs.sampling_bias`
* - $\Delta\gamma = 2\gamma_m-\gamma_b$ (mN/m)
  - {astra:value col=dense err=dense_std where="quantity=delta_gamma" sig=5}`outputs.sampling_bias`
  - {astra:value col=coarse err=coarse_std where="quantity=delta_gamma" sig=5}`outputs.sampling_bias`
  - {astra:value col=shift_sigma where="quantity=delta_gamma" sig=2}`outputs.sampling_bias`
* - $\lambda_m$ (nm)
  - {astra:value col=dense err=dense_std where="quantity=lambda_m" sig=5}`outputs.sampling_bias`
  - {astra:value col=coarse err=coarse_std where="quantity=lambda_m" sig=5}`outputs.sampling_bias`
  - {astra:value col=shift_sigma where="quantity=lambda_m" sig=2}`outputs.sampling_bias`
* - $R$ (nm)
  - {astra:value col=dense err=dense_std where="quantity=R" sig=6}`outputs.sampling_bias`
  - {astra:value col=coarse err=coarse_std where="quantity=R" sig=6}`outputs.sampling_bias`
  - {astra:value col=shift_sigma where="quantity=R" sig=2}`outputs.sampling_bias`
* - $h$ (nm)
  - {astra:value col=dense err=dense_std where="quantity=h" sig=6}`outputs.sampling_bias`
  - {astra:value col=coarse err=coarse_std where="quantity=h" sig=6}`outputs.sampling_bias`
  - {astra:value col=shift_sigma where="quantity=h" sig=2}`outputs.sampling_bias`
* - rms residual (pm)
  - {astra:value col=dense where="quantity=rms" sig=5}`outputs.sampling_bias`
  - {astra:value col=coarse where="quantity=rms" sig=5}`outputs.sampling_bias`
  -
```

The old sampling raised both tensions by {astra:value col=shift_sigma where="quantity=gamma_m" sig=2}`outputs.sampling_bias`σ, lowered $\lambda_m$ and $R$ by about as much, and moved the ratio by {astra:value col=shift_sigma where="quantity=gamma_ratio" sig=2}`outputs.sampling_bias`σ. It also made the Jacobian errors of the tensions about 10 % smaller. The independent code agrees with the current fit and not with the old one (panel d). None of this changes a conclusion: the shift is well inside the statistical error and far below the uncertainty from the assumed $\kappa_m$. It was a bias of the code, though, not of the data, and only the second code revealed it.

:::{astra} findings.sampling_bias_caught
:::

### Across the analysis choices

The comparison runs in every universe of [](#sec-universes). The independent code always gets the same profile, fit range and $\kappa_m$ as ours. With a free bilayer plane it is given the data minus our fitted offset, since it has no offset of its own. The two codes agree equally well in every universe. The old sampling's bias differs between them: with a free plane it shifts the tensions the other way. With bootstrap errors, panels d and e compare against the bootstrap σ.

::::{tab-set}

:::{tab-item} Baseline
![Cross-check, baseline](results/baseline/cross_check_plot.png)

Values: [cross_check.csv](results/baseline/cross_check.csv), [sampling_bias.csv](results/baseline/sampling_bias.csv)
:::

:::{tab-item} 298 K
![Cross-check, 298 K](results/temp-298/cross_check_plot.png)

Values: [cross_check.csv](results/temp-298/cross_check.csv), [sampling_bias.csv](results/temp-298/sampling_bias.csv)
:::

:::{tab-item} Free baseline
![Cross-check, free bilayer plane](results/z0-free/cross_check_plot.png)

Values: [cross_check.csv](results/z0-free/cross_check.csv), [sampling_bias.csv](results/z0-free/sampling_bias.csv)
:::

:::{tab-item} r ≤ 40 nm
![Cross-check, r <= 40 nm](results/r-le-40/cross_check_plot.png)

Values: [cross_check.csv](results/r-le-40/cross_check.csv), [sampling_bias.csv](results/r-le-40/sampling_bias.csv)
:::

:::{tab-item} Bootstrap
![Cross-check, bootstrap errors](results/bootstrap/cross_check_plot.png)

Values: [cross_check.csv](results/bootstrap/cross_check.csv), [sampling_bias.csv](results/bootstrap/sampling_bias.csv)
:::

::::

### What the comparison does not test

Both codes implement the same reading of the paper: the symmetric model, with the factor ½ in the contact-line tension {astra:cite}`prior_insights.contact_tension_factor`. Both read the profile as the upper half of a symmetric lens measured from the bilayer plane, and both assume the same $\kappa_m$. The comparison therefore tests the numerics: solver, root finding, optimiser and error propagation. It does not test the model or the assumptions behind it; those are discussed in [](#sec-discussion).

:::{astra} outputs.cross_check
:caption: The full comparison. `other_std` is the independent code's error as it reports it, `other_std_full` the same covariance propagated with all correlations, and `diff_sigma` = (independent − ours) / our σ.
:::

:::{astra} outputs.sampling_bias
:caption: Our fit with the current and the old contour sampling. `shift_sigma` = (old − current) / our current σ. The last rows give the old sampling's interpolation error at the current best fit.
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
{astra}`decisions.baseline.free`. The bilayer plane gets a fitted offset. Apart from the temperature, which only rescales them, this is the one choice that moves the tensions: by about half a quoted σ.

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
{astra}`decisions.uncertainty.residual_bootstrap`. Same best fit. The errors come from refitting resampled residuals and are close to the Jacobian estimate, at most about 12 % smaller.

![Fit with bootstrap errors](results/bootstrap/fit_plot.png)

Values: [results/bootstrap/fit_summary.csv](results/bootstrap/fit_summary.csv)
:::

::::

:::{astra} findings.robust_to_choices
:::

:::{tip} Trying other combinations
The report renders the first universe in `universes/` (here `baseline`). The decision boxes above show its selections and cannot be switched in the page. To try another combination of choices, add a file to `universes/` and run `lc materialize local`. Its outputs then appear in `results/<universe>/`. See `README.md` for step-by-step instructions.
:::

(sec-discussion)=
## Discussion

The symmetric lens model of @deslandes2017 accounts for the measured profile with three parameters and gives a consistent set of droplet parameters. The tension ratio and volume come straight from the shape. They do not depend on the bending rigidity, the temperature, or any analysis choice we varied. The absolute tensions inherit the assumed $\kappa_m$ linearly. Since monolayer rigidities are known only within a factor of a few {astra:cite}`prior_insights.monolayer_rigidity`, that assumption, rather than the fit, dominates their uncertainty. A measured $\kappa_m$ for the lipids in question would remove it.

The main modelling assumptions are up–down symmetry, a flat bilayer, and no line tension or spontaneous curvature. The paper argues these are sub-dominant for budding. A profile of the lower monolayer, or of the bilayer away from the droplet, would test the symmetry directly. The structured residuals near the contact line point to physics the continuum model leaves out at molecular scales, or to smoothing in how the profile was extracted.

## Data and code availability

The analysis is specified in `astra.yaml` and runs with [Lightcone](https://github.com/LightconeResearch). Inputs, code and every output of every universe live in the repository, with large files in git-annex; `README.md` explains how to get and rerun them. All inputs:

:::{astra} inputs
:::

(sec-code)=
### Code

Every script is below. Click a file name to open it as plain text in a new tab, and save it from there (your browser may add a hash to the file name). Or expand a box to read the code in the page. The analysis specification itself is {download}`astra.yaml <astra.yaml>`, and the Python environment is pinned in {download}`pyproject.toml <pyproject.toml>` and {download}`uv.lock <uv.lock>`.

:::{list-table} Analysis code. Each output is made by exactly one script, run by `lc materialize` with the command recorded in `astra.yaml`.
:label: tbl-code
:header-rows: 1
* - File
  - What it does
  - Makes
* - {download}`src/lens.py <src/lens.py>`
  - Shape solver: shape equations, axis series start, contact line, lens-branch selection
  - used by every output
* - {download}`src/fitting.py <src/fitting.py>`
  - Least-squares model, initial grid search, derived physical quantities
  - used by `fit_params`, `identifiability`
* - {download}`src/fit.py <src/fit.py>`
  - Runs the fit and its uncertainty estimate
  - `fit_params`
* - {download}`src/fit_profile.py <src/fit_profile.py>`
  - Writes the best-fit contour
  - `fit_profile`
* - {download}`src/plot_fit.py <src/plot_fit.py>`
  - Data vs model figure
  - `fit_plot`
* - {download}`src/identifiability.py <src/identifiability.py>`
  - Profile $\Delta\chi^2$ map over the two tensions
  - `identifiability`
* - {download}`src/validate_solver.py <src/validate_solver.py>`
  - Solver against the SI's analytic limits
  - `solver_validation`
* - {download}`src/summarize.py <src/summarize.py>`
  - Flattens the fit into one row per quantity
  - `fit_summary`
* - {download}`src/cross_check.py <src/cross_check.py>`
  - Refits with the independent code and tabulates both
  - `cross_check`
* - {download}`src/old_sampling.py <src/old_sampling.py>`
  - The contour sampling used before the fix, kept to measure its bias
  - used by `sampling_bias`, `cross_check_plot`
* - {download}`src/sampling_bias.py <src/sampling_bias.py>`
  - Reruns the fit with the old sampling and tabulates the shifts
  - `sampling_bias`
* - {download}`src/plot_cross_check.py <src/plot_cross_check.py>`
  - Comparison figure for the two codes
  - `cross_check_plot`
* - {download}`data/raw/compare_fit.py <data/raw/compare_fit.py>`
  - The independent implementation (standalone script)
  - used by `cross_check`, `cross_check_plot`
:::

(code-lens)=
:::{dropdown} `src/lens.py`: shape solver: shape equations, axis series start, contact line, lens-branch selection
```{literalinclude} src/lens.py
:language: python
```
:::

(code-fitting)=
:::{dropdown} `src/fitting.py`: least-squares model, initial grid search, derived physical quantities
```{literalinclude} src/fitting.py
:language: python
```
:::

(code-fit)=
:::{dropdown} `src/fit.py`: runs the fit and its uncertainty estimate
```{literalinclude} src/fit.py
:language: python
```
:::

(code-fit-profile)=
:::{dropdown} `src/fit_profile.py`: writes the best-fit contour
```{literalinclude} src/fit_profile.py
:language: python
```
:::

(code-plot-fit)=
:::{dropdown} `src/plot_fit.py`: data vs model figure
```{literalinclude} src/plot_fit.py
:language: python
```
:::

(code-identifiability)=
:::{dropdown} `src/identifiability.py`: profile $\Delta\chi^2$ map over the two tensions
```{literalinclude} src/identifiability.py
:language: python
```
:::

(code-validate-solver)=
:::{dropdown} `src/validate_solver.py`: solver against the SI's analytic limits
```{literalinclude} src/validate_solver.py
:language: python
```
:::

(code-summarize)=
:::{dropdown} `src/summarize.py`: flattens the fit into one row per quantity
```{literalinclude} src/summarize.py
:language: python
```
:::

(code-cross-check)=
:::{dropdown} `src/cross_check.py`: refits with the independent code and tabulates both
```{literalinclude} src/cross_check.py
:language: python
```
:::

(code-old-sampling)=
:::{dropdown} `src/old_sampling.py`: the contour sampling used before the fix
```{literalinclude} src/old_sampling.py
:language: python
```
:::

(code-sampling-bias)=
:::{dropdown} `src/sampling_bias.py`: reruns the fit with the old sampling
```{literalinclude} src/sampling_bias.py
:language: python
```
:::

(code-plot-cross-check)=
:::{dropdown} `src/plot_cross_check.py`: comparison figure for the two codes
```{literalinclude} src/plot_cross_check.py
:language: python
```
:::

(code-compare-fit)=
:::{dropdown} `data/raw/compare_fit.py`: the independent implementation (standalone script)
```{literalinclude} data/raw/compare_fit.py
:language: python
```
:::
