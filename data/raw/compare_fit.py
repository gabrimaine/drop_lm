#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fit_ld_thiam17.py
=================
Fit a measured lipid-droplet (LD) budding profile with the continuum theory of

  F. Deslandes, A. R. Thiam, L. Foret,
  "Lipid Droplets Can Spontaneously Bud Off from a Symmetric Bilayer",
  Biophysical Journal 113, 15-18 (2017) - Supporting Information.

THE MODEL (symmetric configuration analyzed in the paper)
--------------------------------------------------------
An oil lens is embedded between the two identical monolayers of a flat,
symmetric bilayer:

  * monolayers c (top) and l (bottom): bending rigidity k_m (default
    6 kBT), surface tension gamma_m, zero spontaneous curvature;
  * bilayer b: surface tension gamma_b - exactly flat in the symmetric
    state (its bending rigidity then plays no role);
  * no line tension (tau = 0); cytosol and ER-lumen pressures equal.

Free energy (SI Eqs. 1-2):

  F = sum_{i=c,l,b} [ INT ( k_i/2 (2H - C0i)^2 + kG_i K ) dA + gamma_i A_i ]
      + 2 pi R tau - Pc Vc - Pl Vl .

Minimization at fixed oil volume V yields the shape equations (SI Eq. 18)
for each membrane plus boundary conditions at the symmetry axis and at the
contact line (SI Eqs. 5, 6, 23, 24, 26).

For the symmetric state the problem reduces to the upper monolayer alone.
Its contour is parametrized by arclength s and the tangent angle psi
(SI Eq. 3); z is the height above the bilayer plane. With C0 = 0:

  dr/ds   = cos(psi)
  dz/ds   = -sin(psi)
  dpsi/ds = u
  dt/ds   = (k_m/2)(u^2 - sin^2(psi)/r^2) + gamma_m - p r sin(psi)
  du/ds   = -u cos(psi)/r + cos(psi) sin(psi)/r^2
            - p r cos(psi)/(2 k_m) + t sin(psi)/(k_m r)

where p = P_oil - P_cytosol (Lagrange multiplier of the volume constraint)
and t(s) is the Lagrange multiplier of the constraint dr/ds = cos(psi).
The first integral (SI Eq. 19) is conserved along the contour.

Boundary conditions of the symmetric lens:
  * apex (s=0):    r=0, psi=0, t=0, u=u0  (regularity -> series start);
  * contact (s=L): psi=0  (monolayer leaves the flat bilayer tangentially),
                   z=0    (contact line lies in the bilayer plane),
                   u(L) = -sqrt( (2 gamma_m - gamma_b) / k_m )
                   (contact-line force balance, SI Eq. 23, combined with the
                    first integral, SI Eq. 19; the torque balance SI Eq. 24
                    is satisfied automatically by symmetry).

Outputs of a solved shape: contact radius R, dome height h, and the total
lens volume

  V = 2 pi INT_0^L r^2 sin(psi) ds        (upper + lower half).

DIMENSIONLESS FORM
------------------
Lengths are measured in the elastic length lambda = sqrt(k_m/gamma_m),
t in gamma_m*lambda, p in gamma_m/lambda. The dimensionless shape depends
only on

  p_bar = p lambda / gamma_m   and   delta = 2 - gamma_b/gamma_m ,

while lambda sets the absolute scale. With k_m fixed, (lambda, delta, p_bar)
are the three fit parameters; V, gamma_m, gamma_b and gamma_b/gamma_m follow.

FITTING
-------
The measured profile z_data(r) (r = distance from the droplet axis, assumed
in nm) is fitted by nonlinear least squares over (lambda, delta, p_bar);
the model is the dome for r <= R and the flat bilayer (z = 0) for r > R.
The inner boundary-value problem is solved by shooting on the apex
curvature u0 (1-D root finding); the outer problem is solved by a coarse
grid search followed by trust-region-reflective least squares from several
starting points. Parameter uncertainties come from the least-squares
Jacobian.

Usage:
  python3 fit_ld_thiam17.py --test
  python3 fit_ld_thiam17.py --data td_reset2.xvg --outdir results

Assumptions to remember:
  * data length unit = nm  (rescale with --unit-in-nm if different);
  * kBT at T = 300 K (--temperature); k_m = 6 kBT (--km-kbt);
  * symmetric configuration: gamma_c = gamma_l = gamma_m, flat bilayer,
    tau = 0, C0 = 0 - the "simplest case" of the paper.  If the measured
    droplet is strongly asymmetric (budded), this model does not apply.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

KB = 1.380649e-23  # J/K

# ---------------------------------------------------------------------------
# 1. Dimensionless shape equations and shooting solver
# ---------------------------------------------------------------------------

S0 = 1.0e-3     # dimensionless arclength where the apex series is matched
SMAX = 500.0    # integration cutoff


def _rhs(s, y, pbar):
    """Dimensionless RHS of SI Eq. (18) for membrane c (+ volume state W)."""
    r, z, psi, t, u, W = y
    sp = math.sin(psi)
    cp = math.cos(psi)
    r2 = r * r
    return (cp,
            -sp,
            u,
            0.5 * (u * u - sp * sp / r2) + 1.0 - pbar * r * sp,
            -u * cp / r + cp * sp / r2 - 0.5 * pbar * r * cp + t * sp / r,
            2.0 * math.pi * r2 * sp)


def _apex_series(u0, pbar, s0=S0):
    """Regular Taylor start at the apex (r = psi = t = 0). Dimensionless.

    psi = u0 s + c3 s^3,  t = s + t3 s^3, with
    c3 = (u0 - pbar/2)/8 and t3 = u0^2/12 - 3 pbar u0/8  (from regularity
    of the shape equations at the symmetry axis).
    """
    c3 = (u0 - 0.5 * pbar) / 8.0
    t3 = u0 * u0 / 12.0 - 3.0 * pbar * u0 / 8.0
    return np.array([
        s0 - u0 * u0 * s0 ** 3 / 6.0,                # r
        -(u0 * s0 ** 2) / 2.0 - c3 * s0 ** 4 / 4.0,  # z (apex at z=0)
        u0 * s0 + c3 * s0 ** 3,                      # psi
        s0 + t3 * s0 ** 3,                           # t
        u0 + 3.0 * c3 * s0 ** 2,                     # u
        0.0])                                        # W (volume accumulator)


def _events():
    """Terminal events: psi crossing 0 downward (contact line) or psi
    passing pi/2 (overhang -> not a simple-dome lens branch)."""
    def ev_psi0(s, y, pbar):
        return y[2]
    ev_psi0.terminal = True
    ev_psi0.direction = -1.0

    def ev_over(s, y, pbar):
        return y[2] - (0.5 * math.pi + 1.0e-9)
    ev_over.terminal = True
    ev_over.direction = 1.0
    return [ev_psi0, ev_over]


def integrate_dome(u0, pbar, rtol=1e-10, atol=1e-13, dense=True, n=2500):
    """Integrate the upper monolayer from the apex.

    Returns a dict with the end state (and the dense contour if dense=True),
    or None when the contour does not return to psi = 0 (e.g. it closes into
    a sphere or overhangs past psi = pi/2): no simple-dome solution.
    """
    if not (u0 > 0.0):
        return None
    sol = solve_ivp(_rhs, (S0, SMAX), _apex_series(u0, pbar), args=(pbar,),
                    method="RK45", events=_events(), dense_output=dense,
                    rtol=rtol, atol=atol)
    if sol.status != 1 or len(sol.t_events[0]) == 0:
        return None
    s_end = float(sol.t_events[0][0])
    y_end = sol.y[:, -1]
    out = dict(u0=u0, pbar=pbar, s_end=s_end,
               r_end=float(y_end[0]), z_end=float(y_end[1]),
               psi_end=float(y_end[2]), t_end=float(y_end[3]),
               u_end=float(y_end[4]), W_end=float(y_end[5]))
    if dense:
        s_grid = np.linspace(S0, s_end, n)
        Y = sol.sol(s_grid)
        out.update(s=s_grid, r=Y[0], z=Y[1], psi=Y[2], t=Y[3], u=Y[4], W=Y[5])
        out["psi_max"] = float(np.max(Y[2]))
    return out


def _probe(u0, pbar, delta, rtol):
    """Classify a shooting attempt.

    Returns (0, u_end) for a valid dome (psi returns to 0 without
    overhang) and (1, nan) for a failure (overhang / closing / cutoff).
    """
    d = integrate_dome(u0, pbar, rtol=rtol, dense=False)
    if d is None:
        return 1, math.nan
    return 0, d["u_end"]


def _root_in_interval(a, b, pbar, delta, rtol_fine):
    """Bisection on the classification 'u_end > -sqrt(delta) and valid'.

    Requires: point a is 'above' (valid, u_end > target), b is 'not above'
    (failed, or valid with u_end <= target).  Handles the exponentially
    narrow windows of large cap-like domes, where u_end plunges to -inf
    as u0 approaches the overhang threshold and brentq cannot be used.
    Returns the crossing u0 or None."""
    target = -math.sqrt(delta)
    for _ in range(100):
        m = 0.5 * (a + b)
        if not (a < m < b):
            break
        st, ue = _probe(m, pbar, delta, rtol_fine)
        if st == 0 and ue > target:
            a = m
        else:
            b = m
        if b - a <= max(1e-15, 4.0 * np.finfo(float).eps * abs(b)):
            break
    root = 0.5 * (a + b)
    st, ue = _probe(root, pbar, delta, rtol_fine)
    if st != 0 or abs(ue - target) > 1e-3 * math.sqrt(delta):
        return None
    return root


def find_roots(pbar, delta, rtol_scan=2e-8, rtol_fine=1e-11, n_scan=48):
    """All roots u0 of the contact condition u(L) = -sqrt(delta).

    Two kinds of brackets are searched in a coarse scan of u0:
      (i)  both endpoints valid with u_end on either side of the target
           (regular sign change);
      (ii) a valid->failed transition, across which u_end plunges to
           -infinity (cap-like domes; the root sits just below the
           overhang threshold in an exponentially narrow window).
    """
    target = -math.sqrt(delta)
    grid = np.geomspace(1.5e-3, 4.0, n_scan)
    st = np.empty(len(grid), dtype=int)
    ue = np.empty(len(grid))
    for i, u in enumerate(grid):
        st[i], ue[i] = _probe(u, pbar, delta, rtol_scan)
    roots = []
    for i in range(len(grid) - 1):
        a_ok = st[i] == 0 and ue[i] > target
        b_ok = st[i + 1] == 0 and ue[i + 1] > target
        if a_ok and not b_ok:
            r = _root_in_interval(grid[i], grid[i + 1], pbar, delta, rtol_fine)
        elif b_ok and not a_ok:
            r = _root_in_interval(grid[i + 1], grid[i], pbar, delta, rtol_fine)
        else:
            r = None
        if r is not None:
            roots.append(r)
    # a root may also sit exactly at a grid point
    for i, u in enumerate(grid):
        if st[i] == 0 and abs(ue[i] - target) < 1e-9:
            roots.append(float(u))
    uniq = []
    for r in sorted(roots):
        if not uniq or r - uniq[-1] > 1e-9:
            uniq.append(r)
    return uniq


def find_root_seeded(pbar, delta, u_seed, rtol_scan=2e-8, rtol_fine=1e-11):
    """Root closest to u_seed; tries a cheap local bracket before the full scan."""
    target = -math.sqrt(delta)
    if u_seed and u_seed > 0:
        for widen in (1.3, 2.0, 6.0, 20.0):
            lo = max(u_seed / widen, 1e-4)
            hi = min(u_seed * widen, 4.0)
            if hi <= lo:
                continue
            g = np.geomspace(lo, hi, 15)
            st = np.empty(len(g), dtype=int)
            ue = np.empty(len(g))
            for i, u in enumerate(g):
                st[i], ue[i] = _probe(u, pbar, delta, rtol_scan)
            for i in range(len(g) - 1):
                a_ok = st[i] == 0 and ue[i] > target
                b_ok = st[i + 1] == 0 and ue[i + 1] > target
                if a_ok and not b_ok:
                    r = _root_in_interval(g[i], g[i + 1], pbar, delta, rtol_fine)
                    if r is not None:
                        return r
                elif b_ok and not a_ok:
                    r = _root_in_interval(g[i + 1], g[i], pbar, delta, rtol_fine)
                    if r is not None:
                        return r
    roots = find_roots(pbar, delta, rtol_scan, rtol_fine)
    if not roots:
        return None
    if u_seed:
        return min(roots, key=lambda r: abs(r - u_seed))
    return roots[0]


# ---------------------------------------------------------------------------
# 2. Full shape solution (dimensional) and model evaluation on the data
# ---------------------------------------------------------------------------

def build_profile(sol, lam, delta, pbar):
    """Dimensional quantities from a dimensionless contour."""
    h_bar = -sol["z_end"]                # dome height (dimensionless)
    z_dome = sol["z"] + h_bar            # 0 at contact, h_bar at apex
    return dict(
        lam=lam, delta=delta, pbar=pbar, u0=sol["u0"],
        R=lam * sol["r_end"], h=lam * h_bar,
        V=lam ** 3 * sol["W_end"],
        psi_max=sol["psi_max"], s_end=lam * sol["s_end"],
        r=lam * sol["r"], z=lam * z_dome,
        psi=sol["psi"], u=sol["u"], t=sol["t"],
        r_bar=sol["r"], z_dome_bar=z_dome)


def solve_shape(lam, delta, pbar, u_seed=None, rtol=1e-10):
    """Solve the symmetric-lens BVP; return the dimensional profile or None."""
    u0 = find_root_seeded(pbar, delta, u_seed) if u_seed else None
    if u0 is None:
        roots = find_roots(pbar, delta)
        if not roots:
            return None
        u0 = min(roots, key=lambda r: abs(r - u_seed)) if u_seed else roots[0]
    sol = integrate_dome(u0, pbar, rtol=rtol, dense=True)
    if sol is None:
        return None
    return build_profile(sol, lam, delta, pbar)


def solve_shape_all(lam, delta, pbar, rtol=1e-9):
    """All branch solutions (used by the coarse grid search)."""
    profs = []
    for u0 in find_roots(pbar, delta):
        sol = integrate_dome(u0, pbar, rtol=rtol, dense=True)
        if sol is not None:
            profs.append(build_profile(sol, lam, delta, pbar))
    return profs


def model_z(prof, r_data):
    """Model height at the data radii: dome for r <= R, flat bilayer beyond."""
    rb = np.asarray(r_data, float) / prof["lam"]
    zb = np.interp(rb, prof["r_bar"], prof["z_dome_bar"])
    zb[rb > prof["r_bar"][-1]] = 0.0
    return prof["lam"] * zb


def theta_to_pars(th):
    lam = math.exp(th[0])
    d = math.exp(th[1])
    delta = 2.0 * d / (1.0 + d)
    pbar = math.exp(th[2])
    return lam, delta, pbar


def pars_to_theta(lam, delta, pbar):
    return np.array([math.log(lam),
                     math.log(delta / (2.0 - delta)),
                     math.log(pbar)])


# ---------------------------------------------------------------------------
# 3. Fitting machinery
# ---------------------------------------------------------------------------

def load_profile(path):
    rs, zs = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in "#@;%":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                r, z = float(parts[0]), float(parts[1])
            except ValueError:
                continue
            rs.append(r)
            zs.append(z)
    if len(rs) < 10:
        raise ValueError(f"could not parse enough data points from {path}")
    order = np.argsort(rs)
    return np.asarray(rs, float)[order], np.asarray(zs, float)[order]


def smart_starts(r_data, z_data, km):
    """Initial guesses from the small-volume asymptotics of the SI
    (R^2 = 8 sqrt(k Dg)/p, h = p R^4 / (64 k)), evaluated at the measured
    contact radius and apex height."""
    zmax = float(np.max(z_data))
    rmax = float(np.max(r_data))
    flat = z_data[r_data > 0.6 * rmax]
    noise = float(np.std(flat)) if len(flat) > 5 else 1e-3
    thr = max(5.0 * noise, 0.02)
    above = np.where(z_data > thr)[0]
    R0 = float(r_data[above[-1]]) if len(above) else 0.8 * rmax
    dr = float(np.median(np.diff(r_data))) if len(r_data) > 2 else 0.4
    R0 += 0.5 * dr
    h0 = zmax
    p0 = 64.0 * km * h0 / R0 ** 4            # J/nm^3
    dgamma0 = (p0 * R0 ** 2 / 8.0) ** 2 / km  # J/nm^2
    starts = []
    for lam in (2.5, 3.2, 4.0, 5.0, 6.5):
        delta = min(max(dgamma0 * lam ** 2 / km, 0.05), 1.9)
        pbar = min(max(p0 * lam ** 3 / km, 0.01), 30.0)
        starts.append(pars_to_theta(lam, delta, pbar))
    info = dict(R0=R0, h0=h0, p0=p0, dgamma0=dgamma0, noise=noise)
    return starts, info


def grid_starts(r_data, z_data, km, top=6, verbose=True):
    """Coarse grid over (lambda, delta, p_bar); keep the best starts."""
    lam_list = (1.8, 2.5, 3.5, 5.0, 7.0)
    del_list = (0.7, 1.0, 1.3, 1.6, 1.9)
    pbar_list = (0.1, 0.25, 0.6, 1.5, 3.5)
    results = []
    t0 = time.time()
    for lam in lam_list:
        for delta in del_list:
            for pbar in pbar_list:
                best = None
                for prof in solve_shape_all(lam, delta, pbar, rtol=1e-8):
                    resid = model_z(prof, r_data) - z_data
                    ssr = float(np.sum(resid ** 2))
                    if best is None or ssr < best[0]:
                        best = (ssr, prof)
                if best is not None:
                    results.append((best[0], pars_to_theta(lam, delta, pbar),
                                    best[1]["u0"]))
    results.sort(key=lambda t: t[0])
    if verbose:
        if results:
            print(f"[grid] {len(results)} feasible combos evaluated in "
                  f"{time.time() - t0:.1f} s; best SSR = {results[0][0]:.4g}")
        else:
            print("[grid] no feasible combo")
    return [(th, u0) for _, th, u0 in results[:top]]


def fit_profile(r_data, z_data, starts, verbose=True):
    """Outer nonlinear least squares over theta = (ln lam, logit delta, ln pbar)."""
    N = len(r_data)
    bounds = (np.array([math.log(1.0), math.log(0.02 / 1.98), math.log(2e-3)]),
              np.array([math.log(25.0), math.log(1.98 / 0.02), math.log(80.0)]))
    state = {"u0": None, "th_ref": None}

    def residuals(th):
        lam, delta, pbar = theta_to_pars(th)
        prof = solve_shape(lam, delta, pbar, u_seed=state["u0"])
        if prof is None:
            ref = state["th_ref"] if state["th_ref"] is not None else np.asarray(th)
            base = 1e2 + 0.1 * (np.asarray(th) - ref)
            return np.tile(base, N // 3 + 1)[:N]
        state["u0"] = prof["u0"]
        state["th_ref"] = np.asarray(th, float)
        return model_z(prof, r_data) - z_data

    best = None
    for i, (th0, u0_seed) in enumerate(starts):
        th0 = np.clip(np.asarray(th0, float), bounds[0] + 1e-9, bounds[1] - 1e-9)
        state.update(u0=u0_seed, th_ref=None)
        try:
            res = least_squares(residuals, th0, bounds=bounds, method="trf",
                                x_scale="jac", diff_step=1e-5,
                                ftol=1e-14, xtol=1e-14, gtol=1e-14,
                                max_nfev=400)
        except Exception as exc:  # pragma: no cover
            if verbose:
                print(f"[fit] start {i + 1}/{len(starts)} failed: {exc}")
            continue
        if verbose:
            print(f"[fit] start {i + 1}/{len(starts)}: "
                  f"SSR = {2.0 * res.cost:.6g}, theta = {np.round(res.x, 4)}")
        if best is None or res.cost < best.cost:
            best = res
    if best is None:
        raise RuntimeError("all least-squares starts failed")
    return best


def derived_outputs(prof, km):
    """Physical parameters from a solved profile (units: nm, mN/m, ...)."""
    lam, delta, pbar = prof["lam"], prof["delta"], prof["pbar"]
    gm_J = km / lam ** 2               # J/nm^2
    gm = gm_J * 1e21                   # mN/m
    gb = (2.0 - delta) * gm
    ratio = 2.0 - delta
    p_J = pbar * km / lam ** 3         # J/nm^3
    p_Pa = p_J * 1e27
    V = prof["V"]
    RLD = (3.0 * V / (4.0 * math.pi)) ** (1.0 / 3.0)
    return dict(lam=lam, delta=delta, pbar=pbar,
                gamma_m=gm, gamma_b=gb, ratio=ratio,
                dgamma=(2.0 * gm - gb),
                p_J=p_J, p_Pa=p_Pa, V=V, R_LD=RLD,
                R=prof["R"], h=prof["h"], psi_max=prof["psi_max"],
                R_LD_over_lam=RLD / lam, u0=prof["u0"])


def parameter_errors(res, r_data, z_data, km, prof):
    """1-sigma uncertainties from the least-squares Jacobian + error
    propagation to the physical parameters (V, R, h via numeric gradient)."""
    J = np.atleast_2d(res.jac)
    N, P = J.shape
    rss = 2.0 * res.cost
    sig2 = rss / max(N - P, 1)
    try:
        cov = np.linalg.inv(J.T @ J) * sig2
    except np.linalg.LinAlgError:
        cov = np.zeros((P, P))
    sig_th = np.sqrt(np.maximum(np.diag(cov), 0.0))
    th = res.x
    lam, delta, pbar = theta_to_pars(th)
    sig_lam = lam * sig_th[0]
    sig_delta = delta * (1.0 - delta / 2.0) * sig_th[1]
    sig_pbar = pbar * sig_th[2]
    gm = km / lam ** 2 * 1e21
    sig_gm = 2.0 * gm * sig_lam / lam
    sig_gb = math.hypot(gm * sig_delta, (2.0 - delta) * sig_gm)

    # numeric gradient of derived quantities (V, R, h) wrt theta
    keys = ("V", "R", "h")
    grad = {k: np.zeros(3) for k in keys}
    h = 1e-4
    for j in range(3):
        tp = np.array(th, float); tp[j] += h
        tm = np.array(th, float); tm[j] -= h
        pp = solve_shape(*theta_to_pars(tp), u_seed=prof["u0"])
        pm = solve_shape(*theta_to_pars(tm), u_seed=prof["u0"])
        if pp is None or pm is None:
            continue
        for k in keys:
            grad[k][j] = (pp[k] - pm[k]) / (2.0 * h)
    sig = {}
    for k in keys:
        sig[k] = float(math.sqrt(max(float(grad[k] @ cov @ grad[k]), 0.0)))
    sig.update(sig_lam=sig_lam, sig_delta=sig_delta, sig_pbar=sig_pbar,
               sig_gm=sig_gm, sig_gb=sig_gb, sig_ratio=sig_delta,
               sig_V=sig["V"], sig_R=sig["R"], sig_h=sig["h"],
               cov=cov, sig_th=sig_th, sigma_data=math.sqrt(sig2))
    return sig


# ---------------------------------------------------------------------------
# 4. Validation tests
# ---------------------------------------------------------------------------

def run_tests():
    print("=" * 70)
    print("TEST 1: exact spherical cap (pbar = 2*u0 -> Laplace)")
    u0, pbar = 0.3, 0.6
    sol = solve_ivp(_rhs, (S0, 3.0), _apex_series(u0, pbar), args=(pbar,),
                    dense_output=True, rtol=1e-12, atol=1e-14)
    s = np.linspace(S0, 3.0, 200)
    Y = sol.sol(s)
    psi_ex = u0 * s
    r_ex = np.sin(psi_ex) / u0
    t_ex = np.sin(psi_ex) * np.cos(psi_ex) / u0
    z_ex = (np.cos(psi_ex) - 1.0) / u0
    # W(s) = 2 pi INT r^2 sin(psi) ds = (2 pi/u0^3) INT sin^3 = ...
    W_ex = (2.0 * math.pi / u0 ** 3) * (-np.cos(psi_ex)
                                        + np.cos(psi_ex) ** 3 / 3.0 + 2.0 / 3.0)
    err = max(np.max(np.abs(Y[0] - r_ex)), np.max(np.abs(Y[2] - psi_ex)),
              np.max(np.abs(Y[3] - t_ex)), np.max(np.abs(Y[4] - u0)),
              np.max(np.abs(Y[1] - z_ex)), np.max(np.abs(Y[5] - W_ex)))
    print(f"  max |numerical - exact| = {err:.3e}")
    assert err < 1e-7, "sphere test FAILED"

    print("TEST 2: first integral (SI Eq. 19) conserved along a lens solution")
    roots = find_roots(0.4, 0.8)
    assert roots, "no lens root found for (pbar=0.4, delta=0.8)"
    d = integrate_dome(roots[0], 0.4, dense=True)
    r, psi, t, u = d["r"], d["psi"], d["t"], d["u"]
    H = (0.5 * r * u ** 2 - np.sin(psi) ** 2 / (2.0 * r) - r
         + 0.4 * r ** 2 * np.sin(psi) / 2.0 + t * np.cos(psi))
    print(f"  max |Hamiltonian| = {np.max(np.abs(H)):.3e} (root u0 = {roots[0]:.6f})")
    assert np.max(np.abs(H)) < 1e-6, "Hamiltonian test FAILED"

    print("TEST 3: small-volume asymptotics (SI Eqs. 31-34)")
    cases = [(0.2, 200.0, 0.05), (1.0, 200.0, 0.10), (0.2, 1000.0, 0.02)]
    for delta, pbar, tol in cases:
        roots = find_roots(pbar, delta)
        u_pred = math.sqrt(delta) / 2.0
        u0 = min(roots, key=lambda rr: abs(rr - u_pred))
        d = integrate_dome(u0, pbar, dense=True)
        R_pred = math.sqrt(8.0 * math.sqrt(delta) / pbar)
        h_pred = pbar * R_pred ** 4 / 64.0
        V_pred = (16.0 * math.pi / 3.0) * delta ** 1.5 / pbar ** 2
        eR = abs(d["r_end"] - R_pred) / R_pred
        eh = abs(-d["z_end"] - h_pred) / h_pred
        eV = abs(d["W_end"] - V_pred) / V_pred
        eu = abs(u0 - u_pred) / u_pred
        print(f"  delta={delta:4.2f} pbar={pbar:6.1f}: "
              f"dR/R={eR:7.4f}  dh/h={eh:7.4f}  dV/V={eV:7.4f}  du0/u0={eu:7.4f}")
        assert max(eR, eh, eV, eu) < tol, "small-volume test FAILED"

    print("TEST 4: full pipeline solve_shape(lam=4 nm, delta=0.8, pbar=0.4)")
    prof = solve_shape(4.0, 0.8, 0.4)
    assert prof is not None, "solve_shape failed"
    ok = (prof["R"] > 0 and prof["h"] > 0 and prof["V"] > 0
          and abs(prof["z"][-1]) < 1e-6
          and np.all(np.diff(prof["r_bar"]) > 0)
          and prof["psi_max"] < 0.5 * math.pi)
    print(f"  R = {prof['R']:.3f} nm, h = {prof['h']:.3f} nm, "
          f"V = {prof['V']:.1f} nm^3, psi_max = {math.degrees(prof['psi_max']):.1f} deg")
    assert ok, "pipeline test FAILED"
    print("ALL TESTS PASSED")
    print("=" * 70)


# ---------------------------------------------------------------------------
# 5. Reporting, figure, main
# ---------------------------------------------------------------------------

def make_figure(path, r_data, z_data, prof, out, sig, meta, data_name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(16.0, 12.0))
    lam = prof["lam"]

    # (a) full profile
    a = ax[0, 0]
    a.plot(prof["r"], -prof["z"], "--", color="0.65", lw=1.2,
           label="model, lower half (mirror)")
    a.plot(prof["r"], prof["z"], "-", color="crimson", lw=2.2,
           label="model, upper monolayer")
    a.plot(r_data, z_data, "o", ms=3.5, color="steelblue", alpha=0.8,
           label="data")
    a.axhline(0.0, color="k", lw=0.8)
    a.axvline(prof["R"], color="gray", ls=":", lw=1.0)
    a.annotate(f"R = {out['R']:.2f} nm", xy=(prof["R"], 0.0),
               xytext=(prof["R"] + 1.0, 1.2), fontsize=10, color="gray")
    a.annotate(f"h = {out['h']:.2f} nm", xy=(0.0, out["h"]),
               xytext=(1.5, out["h"] + 0.3), fontsize=10, color="gray")
    a.set_xlabel("r (nm)")
    a.set_ylabel("z (nm)")
    a.set_title("(a) LD profile: data vs. symmetric-lens theory")
    a.legend(loc="upper right", fontsize=10)
    txt = (f"V = {out['V']:.0f} nm$^3$  ($R_{{LD}}$ = {out['R_LD']:.2f} nm)\n"
           f"$\\gamma_m$ = {out['gamma_m']:.3f} mN/m\n"
           f"$\\gamma_b$ = {out['gamma_b']:.3f} mN/m\n"
           f"$\\gamma_b/\\gamma_m$ = {out['ratio']:.3f}\n"
           f"$\\lambda$ = {lam:.2f} nm,  $p$ = {out['p_Pa'] / 1e3:.1f} kPa")
    a.text(0.985, 0.60, txt, transform=a.transAxes, ha="right", va="top",
           fontsize=11, bbox=dict(boxstyle="round", fc="ivory", ec="0.7"))

    # (b) zoom on the contact region
    b = ax[0, 1]
    m = (r_data > prof["R"] - 8.0) & (r_data < prof["R"] + 8.0)
    b.plot(r_data[m], z_data[m], "o", ms=5, color="steelblue",
           label="data")
    b.plot(prof["r"], prof["z"], "-", color="crimson", lw=2.0, label="model")
    b.axhline(0.0, color="k", lw=0.8)
    b.axvline(prof["R"], color="gray", ls=":", lw=1.0)
    b.set_xlabel("r (nm)")
    b.set_ylabel("z (nm)")
    b.set_title("(b) contact-line region (bending boundary layer)")
    b.legend(fontsize=10)

    # (c) residuals
    c = ax[1, 0]
    zm = model_z(prof, r_data)
    resid = z_data - zm
    c.plot(r_data, resid, "o-", ms=3.5, lw=0.8, color="seagreen")
    c.axhline(0.0, color="k", lw=0.8)
    sd = sig["sigma_data"]
    c.axhspan(-2 * sd, 2 * sd, color="gold", alpha=0.35,
              label=f"$\\pm 2\\sigma$ ($\\sigma$ = {sd:.4f} nm)")
    c.axvline(prof["R"], color="gray", ls=":", lw=1.0)
    c.set_xlabel("r (nm)")
    c.set_ylabel("z(data) - z(model) (nm)")
    c.set_title("(c) fit residuals")
    c.legend(fontsize=10)

    # (d) tangent angle of the fitted monolayer
    d_ = ax[1, 1]
    d_.plot(prof["r"], np.degrees(prof["psi"]), color="darkorchid", lw=2.0)
    d_.axvline(prof["R"], color="gray", ls=":", lw=1.0)
    d_.set_xlabel("r (nm)")
    d_.set_ylabel(r"$\psi$ (deg)")
    d_.set_ylim(-3, math.degrees(prof["psi_max"]) * 1.25)
    d_.set_title("(d) tangent angle of the fitted upper monolayer\n"
                 r"$\psi = 0$ at apex and at the contact line")
    d_.annotate(f"$\\psi_{{max}}$ = {math.degrees(prof['psi_max']):.1f}$^\\circ$",
                xy=(float(prof["r"][np.argmax(prof["psi"])]),
                    math.degrees(prof["psi_max"])),
                xytext=(prof["R"] * 0.45, math.degrees(prof["psi_max"]) * 0.45),
                arrowprops=dict(arrowstyle="->", color="0.4"), fontsize=11)

    fig.suptitle("Lipid-droplet shape fit - theory of Deslandes, Thiam & Foret, "
                 f"Biophys. J. 113:15 (2017), symmetric model | data: {data_name}",
                 fontsize=14)
    fig.text(0.5, 0.012,
             f"fixed $k_m$ = {meta['km_kbt']:.1f} $k_BT$ (T = {meta['T']:.0f} K), "
             f"$k_BT$ = {meta['kBT']:.3e} J | RMS residual = {meta['rms']:.4f} nm | "
             f"$R^2$ = {meta['r2']:.5f}",
             ha="center", fontsize=11, color="0.35")
    fig.tight_layout(rect=(0, 0.025, 1, 0.97))
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Fit an LD budding profile with the Deslandes-Thiam-Foret "
                    "(2017) symmetric-lens theory.")
    ap.add_argument("--data", required=False, default=None)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--km-kbt", type=float, default=6.0,
                    help="monolayer bending rigidity in units of kBT (default 6)")
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--unit-in-nm", type=float, default=1.0,
                    help="length unit of the data expressed in nm (default 1)")
    args = ap.parse_args(argv)

    if args.test:
        run_tests()
        return

    if not args.data:
        ap.error("--data is required (or use --test)")

    kBT = KB * args.temperature
    km = args.km_kbt * kBT

    r_data, z_data = load_profile(args.data)
    r_data = r_data * args.unit_in_nm
    z_data = z_data * args.unit_in_nm

    outdir = args.outdir or os.path.dirname(os.path.abspath(__file__))
    os.makedirs(outdir, exist_ok=True)

    print(f"data   : {args.data}  ({len(r_data)} points, unit = "
          f"{args.unit_in_nm} nm)")
    print("theory : Deslandes, Thiam & Foret, Biophys. J. 113:15 (2017), SI - "
          "symmetric model")
    print(f"fixed  : k_m = {args.km_kbt:.1f} kBT = {km:.4e} J "
          f"(T = {args.temperature:.1f} K)")

    smart, sinfo = smart_starts(r_data, z_data, km)
    print(f"[data] contact radius ~ {sinfo['R0']:.2f} nm, apex height ~ "
          f"{sinfo['h0']:.2f} nm, flat-region noise ~ {sinfo['noise']:.4f} nm")
    print(f"[data] small-volume guess: p ~ {sinfo['p0']:.3e} J/nm^3, "
          f"Dgamma = 2*gamma_m - gamma_b ~ {sinfo['dgamma0'] * 1e21:.3f} mN/m")

    starts = [(th, None) for th in smart]
    starts += grid_starts(r_data, z_data, km, top=6)
    print(f"[fit] {len(starts)} starting points "
          f"({len(smart)} smart + 6 grid)")

    t0 = time.time()
    res = fit_profile(r_data, z_data, starts)
    print(f"[fit] done in {time.time() - t0:.1f} s")

    # final evaluation: among all branch solutions at the optimum, take the
    # one with the lowest SSR (robust to branch multiplicity)
    lam, delta, pbar = theta_to_pars(res.x)
    roots = find_roots(pbar, delta)
    print(f"[fit] optimum: lambda = {lam:.4f} nm, delta = {delta:.4f}, "
          f"p_bar = {pbar:.4f}; {len(roots)} branch root(s) "
          f"{[round(rr, 6) for rr in roots]}")
    prof = None
    best_ssr = None
    for cand in solve_shape_all(lam, delta, pbar):
        ssr = float(np.sum((model_z(cand, r_data) - z_data) ** 2))
        if best_ssr is None or ssr < best_ssr:
            best_ssr, prof = ssr, cand
    if prof is None:
        raise RuntimeError("failed to rebuild the profile at the optimum")
    zm = model_z(prof, r_data)
    resid = z_data - zm
    rss = float(np.sum(resid ** 2))
    rms = math.sqrt(rss / len(r_data))
    tss = float(np.sum((z_data - z_data.mean()) ** 2))
    r2 = 1.0 - rss / tss

    out = derived_outputs(prof, km)
    sig = parameter_errors(res, r_data, z_data, km, prof)

    laplace_apex = 2.0 * (km / lam ** 2) * (prof["u0"] / lam)  # J/nm^3
    print()
    print("=" * 74)
    print("FIT RESULTS  (symmetric model; k_m = "
          f"{args.km_kbt:.1f} kBT fixed)")
    print("-" * 74)
    print(f"  Droplet volume        V        = {out['V']:10.1f} +/- "
          f"{sig['sig_V']:.1f}  nm^3")
    print(f"  Effective LD radius   R_LD     = {out['R_LD']:10.3f} +/- "
          f"{sig['sig_R']:.3f}  nm   (R_LD/lambda = {out['R_LD_over_lam']:.2f})")
    print(f"  gamma_b / gamma_m              = {out['ratio']:10.4f} +/- "
          f"{sig['sig_ratio']:.4f}")
    print(f"  Monolayer tension     gamma_m  = {out['gamma_m']:10.4f} +/- "
          f"{sig['sig_gm']:.4f}  mN/m")
    print(f"  Bilayer tension       gamma_b  = {out['gamma_b']:10.4f} +/- "
          f"{sig['sig_gb']:.4f}  mN/m")
    print(f"  (2 gamma_m - gamma_b) Dgamma   = {out['dgamma']:10.4f} mN/m")
    print(f"  Elastic length        lambda   = {out['lam']:10.4f} +/- "
          f"{sig['sig_lam']:.4f}  nm   (sqrt(k_m/gamma_m))")
    print(f"  Oil pressure offset   p        = {out['p_J']:.4e} J/nm^3 = "
          f"{out['p_Pa'] / 1e3:.2f} kPa   (Laplace at apex: "
          f"{laplace_apex * 1e27 / 1e3:.2f} kPa)")
    print(f"  Contact radius        R        = {out['R']:10.3f} +/- "
          f"{sig['sig_R']:.3f}  nm   (data ~ {sinfo['R0']:.2f} nm)")
    print(f"  Apex height           h        = {out['h']:10.3f} +/- "
          f"{sig['sig_h']:.3f}  nm   (data ~ {sinfo['h0']:.2f} nm)")
    print(f"  Max tangent angle     psi_max  = "
          f"{math.degrees(out['psi_max']):10.2f} deg  (< 90 deg: simple dome)")
    print(f"  Fit quality: RMS = {rms:.4f} nm (data noise ~ "
          f"{sinfo['noise']:.4f} nm), R^2 = {r2:.5f}, N = {len(r_data)}")
    if out["ratio"] < 0.17:
        note = ("gamma_b/gamma_m < 0.17: symmetric lens would be UNSTABLE -> "
                "budded branch expected (paper Fig. 3B)")
    else:
        note = ("gamma_b/gamma_m > 0.17: symmetric lens is the STABLE shape "
                "at this size (paper Fig. 3B)")
    print(f"  Phase-diagram position ({out['ratio']:.3f}, "
          f"{out['R_LD_over_lam']:.2f}): {note}")
    print("=" * 74)

    # ---------------- outputs ----------------
    meta = dict(km_kbt=args.km_kbt, T=args.temperature, kBT=kBT,
                rms=rms, r2=r2)
    fig_path = os.path.join(outdir, "ld_fit_figure.png")
    data_name = os.path.basename(args.data)
    make_figure(fig_path, r_data, z_data, prof, out, sig, meta, data_name)

    results = dict(
        data_file=args.data, n_points=int(len(r_data)),
        unit_in_nm=args.unit_in_nm, temperature_K=args.temperature,
        km_J=float(km), km_kbt=float(args.km_kbt), kBT_J=float(kBT),
        theory="Deslandes, Thiam & Foret, Biophys. J. 113:15 (2017), SI - "
               "symmetric model (gamma_c=gamma_l=gamma_m, flat bilayer, "
               "tau=0, C0=0)",
        fit_parameters=dict(lam_nm=float(lam), delta=float(delta),
                            p_bar=float(pbar), u0_dimless=float(prof["u0"])),
        extracted=dict(
            V_nm3=float(out["V"]), V_nm3_err=float(sig["sig_V"]),
            R_LD_nm=float(out["R_LD"]), R_LD_nm_err=float(sig["sig_R"]),
            gamma_b_over_gamma_m=float(out["ratio"]),
            gamma_b_over_gamma_m_err=float(sig["sig_ratio"]),
            gamma_m_mN_m=float(out["gamma_m"]),
            gamma_m_mN_m_err=float(sig["sig_gm"]),
            gamma_b_mN_m=float(out["gamma_b"]),
            gamma_b_mN_m_err=float(sig["sig_gb"]),
            Delta_gamma_mN_m=float(out["dgamma"]),
            lambda_nm=float(lam), lambda_nm_err=float(sig["sig_lam"]),
            p_J_per_nm3=float(out["p_J"]), p_Pa=float(out["p_Pa"]),
            R_contact_nm=float(out["R"]), h_nm=float(out["h"]),
            psi_max_deg=float(math.degrees(out["psi_max"])),
        ),
        quality=dict(rms_nm=float(rms), R2=float(r2),
                     sigma_data_nm=float(sig["sigma_data"]),
                     SSR=float(rss)),
        diagnostics=dict(
            n_branch_roots=len(roots),
            data_R_est_nm=sinfo["R0"], data_h_est_nm=sinfo["h0"],
            data_noise_nm=sinfo["noise"],
            theta=res.x.tolist(),
            sigma_theta=sig["sig_th"].tolist(),
        ),
        files=dict(figure=fig_path),
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
    )
    json_path = os.path.join(outdir, "ld_fit_results.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)

    prof_path = os.path.join(outdir, "fitted_profile_model_curve.txt")
    with open(prof_path, "w") as f:
        f.write("# Fitted symmetric-lens monolayer profile "
                "(Deslandes-Thiam-Foret 2017 model)\n")
        f.write(f"# lambda={lam:.6f} nm  delta={delta:.6f}  p_bar={pbar:.6f}  "
                f"R={out['R']:.6f} nm  h={out['h']:.6f} nm  V={out['V']:.3f} "
                f"nm^3\n")
        f.write("# columns: r_nm  z_nm  psi_rad  mean_curvature_2H_nm^-1\n")
        for rr, zz, ps, uu in zip(prof["r"], prof["z"], prof["psi"], prof["u"]):
            twoH = uu + math.sin(ps) / rr  # 2H in 1/nm
            f.write(f"{rr:.6f} {zz:.6f} {ps:.8f} {twoH:.8f}\n")

    comp_path = os.path.join(outdir, "data_with_model.txt")
    with open(comp_path, "w") as f:
        f.write("# data vs model at data radii\n"
                "# columns: r_nm  z_data_nm  z_model_nm  residual_nm\n")
        for rr, zz, mm in zip(r_data, z_data, zm):
            f.write(f"{rr:.6f} {zz:.6f} {mm:.6f} {zz - mm:.6f}\n")

    print(f"\noutputs:\n  {fig_path}\n  {json_path}\n  {prof_path}\n  {comp_path}")
    return results


if __name__ == "__main__":
    main()

