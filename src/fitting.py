"""Least-squares fit of the symmetric-lens model to a measured z(r) profile.

Free parameters theta = (ln lambda[nm], g = gamma_b/gamma_m, ln a[1/lambda]) and,
optionally, a baseline offset z0[nm].  With kappa_m fixed, lambda gives
gamma_m = kappa_m / lambda^2; the apex curvature a labels the position along
the equilibrium branch (equivalently the volume; P is solved for).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from lens import lens, profile_on_grid, to_physical, KB

TEMPERATURES = {"t310": 310.0, "t298": 298.15}
FIT_RANGES = {"all": np.inf, "r_le_40": 40.0}


def load_profile(path, fit_range):
    r, z = np.loadtxt(path, comments=("#", "@"), usecols=(0, 1), unpack=True)
    m = r <= FIT_RANGES[fit_range]
    return r[m], z[m]


class LensModel:
    def __init__(self, r, z, free_baseline: bool):
        self.r, self.z = r, z
        self.free_baseline = free_baseline
        self._x = {}  # cache of x = P/2a - 1 keyed on rounded (g, a) for continuation
        self._last_x = None
        self.n_eval = 0

    def shape(self, g, a):
        sh, x = lens(a, g, x_guess=self._last_x)
        if sh is not None:
            self._last_x = x
        return sh

    def predict(self, theta, r=None):
        r = self.r if r is None else r
        lam, g, a = np.exp(theta[0]), theta[1], np.exp(theta[2])
        z0 = theta[3] if self.free_baseline else 0.0
        sh = self.shape(g, a)
        if sh is None:
            return None, None
        return profile_on_grid(sh, lam, r) + z0, sh

    def residuals(self, theta):
        self.n_eval += 1
        zm, _ = self.predict(theta)
        if zm is None:
            return np.full_like(self.z, 10.0 * np.abs(self.z).max() + 10.0)
        return zm - self.z


def initial_guess(model: LensModel):
    """Coarse grid over (g, a); for each shape scan lambda; keep the best."""
    r, z = model.r, model.z
    best = (np.inf, None)
    for g in np.linspace(0.1, 1.8, 12):
        x = None
        for a in np.geomspace(0.95 * np.sqrt(2 - g) / 2, 0.06, 16):
            sh, x_new = lens(a, g, x_guess=x)
            if sh is None:
                continue
            x = x_new
            for lam in np.geomspace(1.0, 30.0, 80):
                res = profile_on_grid(sh, lam, r) - z
                c = float(res @ res)
                if c < best[0]:
                    best = (c, (np.log(lam), g, np.log(a)))
    theta = list(best[1])
    if model.free_baseline:
        theta.append(0.0)
    return np.array(theta)


def fit(model: LensModel, theta0=None):
    theta0 = initial_guess(model) if theta0 is None else np.asarray(theta0, float)
    lo = [np.log(0.3), 1e-3, np.log(1e-3)] + ([-1.0] if model.free_baseline else [])
    hi = [np.log(100.0), 2.0 - 1e-3, np.log(5.0)] + ([1.0] if model.free_baseline else [])
    model._last_x = None
    model.shape(theta0[1], np.exp(theta0[2]))  # seed continuation
    res = least_squares(model.residuals, theta0, bounds=(lo, hi), method="trf",
                        x_scale="jac", diff_step=1e-6, xtol=1e-12, ftol=1e-12, gtol=1e-12,
                        max_nfev=400)
    return res


def derived(model: LensModel, theta, kappa):
    """Physical quantities for parameter vector theta (kappa in J)."""
    lam_nm, g, a = np.exp(theta[0]), theta[1], np.exp(theta[2])
    sh = model.shape(g, a)
    if sh is None:
        return None
    ph = to_physical(sh, lam_nm, kappa)
    gm = ph["gamma_m_mN_per_m"]
    out = dict(
        V_nm3=ph["V_nm3"],
        gamma_ratio=g,
        gamma_m_mN_per_m=gm,
        gamma_b_mN_per_m=g * gm,
        delta_gamma_mN_per_m=(2 - g) * gm,
        lambda_nm=lam_nm,
        R_nm=ph["R_nm"],
        h_nm=ph["h_nm"],
        R_LD_nm=(3 * ph["V_nm3"] / (4 * np.pi)) ** (1 / 3),
        V_over_lambda3=ph["V_nm3"] / lam_nm**3,
        p_bar=ph["p_Pa"] / 1e5,
        psi_max_deg=ph["psi_max_deg"],
        large_volume_contact_angle_deg=float(np.degrees(np.arccos(g / 2))),
        apex_curvature_inv_nm=a / lam_nm,
    )
    if model.free_baseline:
        out["z0_nm"] = float(theta[3])
    return out


def delta_method(model, theta, cov, kappa, keys, h=1e-5):
    base = derived(model, theta, kappa)
    J = np.zeros((len(keys), len(theta)))
    for j in range(len(theta)):
        tp, tm = theta.copy(), theta.copy()
        tp[j] += h
        tm[j] -= h
        dp, dm = derived(model, tp, kappa), derived(model, tm, kappa)
        for i, k in enumerate(keys):
            J[i, j] = (dp[k] - dm[k]) / (2 * h)
    C = J @ cov @ J.T
    return base, C
