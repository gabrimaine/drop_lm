"""Symmetric lipid-droplet lens shapes from the Deslandes-Thiam-Foret (BJ 2017) model.

Symmetric model of SI Sect. 5 (identical leaflets, P_c = P_l, tau = 0, C0 = 0):
the bilayer b stays flat (psi_b = 0, z_b = 0), psi_contact = 0, and the shape is
set by ONE monolayer obeying the shape equations SI Eq. (18):

    r' = cos psi,  z' = -sin psi,  psi' = u,
    u' = -u cos psi / r + sin psi cos psi / r^2 - p r cos psi / (2 kappa) + t sin psi / (kappa r),
    t' = kappa/2 (u^2 - sin^2 psi / r^2) + gamma_m - p r sin psi,

(' = d/ds, s = arclength), with axis conditions r = psi = 0 (SI Eq. 5) and the
contact line where psi returns to 0 (SI Eq. 6 with psi_contact = 0). The lateral
force balance at the contact line, SI Eq. (23), 2 t_m(R) = t_b(R) = gamma_b R,
combined with the Hamiltonian constraint SI Eq. (19) at psi = 0,
t_m(R) = gamma_m R - kappa R u(R)^2 / 2, gives

    kappa u(R)^2 = 2 gamma_m - gamma_b .                       (contact condition)

The torque balance SI Eq. (24) is satisfied identically by symmetry, and the
bilayer rigidity kappa_b and the Gaussian moduli drop out of the symmetric shape.

Everything below is dimensionless: lengths in lambda = sqrt(kappa/gamma_m),
tensions in gamma_m, pressure P = p lambda^3 / kappa.  g = gamma_b / gamma_m.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

KB = 1.380649e-23  # J/K


def rhs(s, y, P):
    r, z, psi, u, t = y
    c, sn = np.cos(psi), np.sin(psi)
    du = -u * c / r + sn * c / r**2 - 0.5 * P * r * c + t * sn / r
    dt = 0.5 * (u**2 - (sn / r) ** 2) + 1.0 - P * r * sn
    return [c, -sn, u, du, dt]


def axis_start(a, P, s0=1e-3):
    """Series solution near the axis: psi = a s + b s^3 with b = (a - P/2)/8."""
    b = (a - 0.5 * P) / 8.0
    return [s0 - a**2 * s0**3 / 6.0, 0.0, a * s0 + b * s0**3, a + 3 * b * s0**2, s0]


@dataclass
class Shape:
    a: float          # apex curvature (1/lambda)
    P: float          # dimensionless pressure
    ok: bool          # reached a contact line psi -> 0 with psi' < 0
    s: np.ndarray
    r: np.ndarray
    z: np.ndarray     # height above the bilayer plane (z(R) = 0)
    psi: np.ndarray
    u: np.ndarray
    t: np.ndarray

    @property
    def R(self):
        return self.r[-1]

    @property
    def uR(self):
        return self.u[-1]

    @property
    def psi_max(self):
        return float(np.max(self.psi))

    def volume(self):
        """Total oil volume V = 2 * 2 pi * int r z dr (both leaflets), in lambda^3."""
        return 2.0 * np.trapezoid(2 * np.pi * self.r * self.z * np.cos(self.psi), self.s)

    def hamiltonian(self, P):
        """SI Eq. (19); must stay ~0 along an exact solution."""
        r, psi, u, t = self.r, self.psi, self.u, self.t
        return 0.5 * r * (u**2 - (np.sin(psi) / r) ** 2) - r + 0.5 * P * r**2 * np.sin(psi) + t * np.cos(psi)


def shoot(a, P, s_max=200.0, dense=False, rtol=1e-10, atol=1e-12):
    """Integrate from the axis until psi returns to 0 (contact line)."""

    def contact(s, y, P):
        return y[2]

    contact.terminal = True
    contact.direction = -1

    def overturn(s, y, P):  # psi reaching pi: not a lens
        return y[2] - np.pi

    overturn.terminal = True

    def axis_hit(s, y, P):  # r back to ~0
        return y[0] - 1e-6

    axis_hit.terminal = True
    axis_hit.direction = -1

    y0 = axis_start(a, P)
    sol = solve_ivp(rhs, (y0[0], s_max), y0, args=(P,), method="DOP853",
                    rtol=rtol, atol=atol, events=(contact, overturn, axis_hit),
                    dense_output=False, max_step=0.05 if dense else np.inf)
    ok = sol.status == 1 and len(sol.t_events[0]) > 0
    s, Y = sol.t, sol.y
    if ok:  # replace last point by the exact event state
        s = np.append(s[:-1], sol.t_events[0][0])
        Y = np.column_stack([Y[:, :-1], sol.y_events[0][0]])
    z = Y[1] - Y[1, -1]
    return Shape(a, P, ok, s, Y[0], z, Y[2], Y[3], Y[4])


def contact_residual(a, P, g):
    sh = shoot(a, P)
    if not sh.ok:
        return np.nan
    return sh.uR + np.sqrt(2.0 - g)  # u(R) < 0 on the lens branch


def _res_x(x, a, g):
    return contact_residual(a, 2.0 * a * (1.0 + x), g)


def solve_pressure(a, g, x_lo=1e-9, x_hi=1e3, n=90):
    """For a given apex curvature a, find all P with a valid contact line.

    P = 2a is the closed sphere of radius 1/a; for P < 2a the contour overturns
    (psi -> pi) and for P > 2a it comes back to psi = 0 at a contact line. So
    the search runs on x = P/(2a) - 1 > 0 on a log grid, then bisects.
    """
    xs = np.geomspace(x_lo, x_hi, n)
    f = np.array([_res_x(x, a, g) for x in xs])
    roots = []
    for i in range(len(xs) - 1):
        if np.isfinite(f[i]) and np.isfinite(f[i + 1]) and np.sign(f[i]) != np.sign(f[i + 1]):
            x = brentq(_res_x, xs[i], xs[i + 1], args=(a, g), xtol=1e-15, rtol=1e-13)
            roots.append(2.0 * a * (1.0 + x))
    return roots


def lens(a, g, x_guess=None):
    """Equilibrium symmetric LENS for apex curvature a and tension ratio g.

    Of the roots found by solve_pressure, the lens is the one with the largest
    x = P/(2a) - 1 and a single-valued profile (psi_max < pi/2); the other root
    is an overhanging, nearly closed sphere.  Returns (Shape, x) or (None, None).
    x_guess enables a cheap local bracket around a nearby solution (fitting).
    """
    if not (0.0 < g < 2.0) or a <= 0:
        return None, None
    if x_guess is not None and x_guess > 0:
        for w in (1.02, 1.1, 1.5, 3.0, 10.0):
            lo, hi = x_guess / w, x_guess * w
            flo, fhi = _res_x(lo, a, g), _res_x(hi, a, g)
            if np.isfinite(flo) and np.isfinite(fhi) and np.sign(flo) != np.sign(fhi):
                x = brentq(_res_x, lo, hi, args=(a, g), xtol=1e-15, rtol=1e-13)
                sh = shoot(a, 2.0 * a * (1.0 + x), dense=True)
                if sh.ok and sh.psi_max < np.pi / 2:
                    return sh, x
                break
    roots = solve_pressure(a, g)
    for P in sorted(roots, reverse=True):
        sh = shoot(a, P, dense=True)
        if sh.ok and sh.psi_max < np.pi / 2:
            return sh, P / (2.0 * a) - 1.0
    return None, None


# ---------------------------------------------------------------- dimensional

def kappa_J(kappa_kT, T):
    return kappa_kT * KB * T


def to_physical(shape: Shape, lam_nm, kappa):
    """Return dimensional quantities. lam in nm, kappa in J."""
    lam = lam_nm * 1e-9
    gm = kappa / lam**2
    return dict(
        r_nm=shape.r * lam_nm,
        z_nm=shape.z * lam_nm,
        R_nm=shape.R * lam_nm,
        h_nm=2 * shape.z[0] * lam_nm,
        V_nm3=shape.volume() * lam_nm**3,
        gamma_m_mN_per_m=gm * 1e3,
        p_Pa=shape.P * kappa / lam**3,
        psi_max_deg=np.degrees(shape.psi_max),
    )


def profile_on_grid(shape: Shape, lam_nm, r_nm):
    """Model z(r) on the data grid: z = 0 (flat bilayer) beyond the contact line.

    Requires a single-valued profile (psi < pi/2 everywhere), which holds for the
    lenses considered here; overhanging shapes are rejected by the caller.
    """
    rm, zm = shape.r * lam_nm, shape.z * lam_nm
    out = np.zeros_like(r_nm, dtype=float)
    inside = r_nm <= rm[-1]
    out[inside] = np.interp(r_nm[inside], rm, zm)
    return out
