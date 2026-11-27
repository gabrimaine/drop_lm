"""Validate the shape solver against the analytic limits of the SI.

Small volume (SI Sect. 6, Eqs. 31-34), in units lambda = sqrt(kappa/gamma_m):
    R = (12 V / pi)^(1/4) Delta^(-1/8),  h = sqrt(3/(4 pi)) Delta^(1/4) V^(1/2),
    z(r) = P/64 (R^2 - r^2)^2,           Delta = 2 - gamma_b/gamma_m.
Large volume (SI Sect. 7, Eq. 40): spherical caps meeting at cos psi_c = g/2,
approached by psi_max as V/lambda^3 -> infinity.
Also reports the Hamiltonian first integral (SI Eq. 19) along every solution.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from lens import lens  # noqa: E402


def branch(g, a_min=0.06, n=40):
    a_hi = 0.98 * np.sqrt(2 - g) / 2  # small-volume apex limit is sqrt(Delta)/2
    out, x = [], None
    for a in np.geomspace(a_hi, a_min, n):
        sh, x_new = lens(a, g, x_guess=x)
        if sh is None:
            continue
        x = x_new
        out.append(sh)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--g-values", default="0.1,1.0")
    args = ap.parse_args()
    gs = [float(v) for v in args.g_values.split(",")]

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    worstH = 0.0
    for g, col in zip(gs, ["C0", "C3", "C2", "C1"]):
        D = 2 - g
        shapes = branch(g)
        V = np.array([s.volume() for s in shapes])
        R = np.array([s.R for s in shapes])
        h = np.array([2 * s.z[0] for s in shapes])
        pm = np.degrees([s.psi_max for s in shapes])
        worstH = max(worstH, max(np.abs(s.hamiltonian(s.P)).max() for s in shapes))
        Vg = np.geomspace(V.min(), V.max(), 100)
        ax[0, 0].loglog(V, R, "o", ms=3, color=col, label=f"numeric, g={g}")
        ax[0, 0].loglog(Vg, (12 * Vg / np.pi) ** 0.25 * D ** (-1 / 8), "--", color=col, label="SI Eq. 33")
        ax[0, 1].loglog(V, h, "o", ms=3, color=col, label=f"numeric, g={g}")
        ax[0, 1].loglog(Vg, np.sqrt(3 / (4 * np.pi)) * D**0.25 * Vg**0.5, "--", color=col, label="SI Eq. 34")
        # Laplace sphere of radius 2/P through the apex meets z=0 at cos(theta) = 1 - z0 P/2
        th = np.degrees(np.arccos(np.clip(1 - np.array([s.z[0] * s.P / 2 for s in shapes]), -1, 1)))
        ax[1, 0].semilogx(V, pm, "o-", ms=3, color=col, label=f"psi_max, g={g}")
        ax[1, 0].semilogx(V, th, "s-", ms=3, mfc="none", color=col, label=f"apparent cap angle, g={g}")
        ax[1, 0].axhline(np.degrees(np.arccos(g / 2)), ls="--", color=col, label=f"SI Eq. 40: {np.degrees(np.arccos(g/2)):.1f} deg")
        # small-volume profile shape vs Eq. 31 (normalised)
        for s_, ls in ((shapes[0], "-"), (shapes[min(8, len(shapes) - 1)], "-.")):
            rr = np.linspace(0, 1, 100)
            ax[1, 1].plot(s_.r / s_.R, s_.z / s_.z[0], ls, color=col, label=f"g={g}, V={s_.volume():.2g}")
        ax[1, 1].plot(rr, (1 - rr**2) ** 2, ":", color="k", lw=2)
        # large volume profile vs spherical cap of Laplace radius 2/P
        s1 = shapes[-1]
        Rs = 2 / s1.P
        zc = s1.z[0] - Rs
        tt = np.linspace(0, np.pi / 2, 200)
        xc, yc = Rs * np.sin(tt), zc + Rs * np.cos(tt)
        m = yc >= 0
        ax[0, 2].plot(s1.r, s1.z, color=col, label=f"numeric g={g}, V={s1.volume():.3g}")
        ax[0, 2].plot(xc[m], yc[m], ":", color="k")
        ax[1, 2].plot(s1.r, s1.u, color=col, label=f"g={g}")
        ax[1, 2].axhline(-np.sqrt(D), ls="--", color=col, lw=0.8)
    ax[0, 0].set(xlabel=r"$V/\lambda^3$", ylabel=r"$R/\lambda$", title="contact radius")
    ax[0, 1].set(xlabel=r"$V/\lambda^3$", ylabel=r"$h/\lambda$", title="total thickness")
    ax[1, 0].set(xlabel=r"$V/\lambda^3$", ylabel="angle (deg)", title="large-V limit: contact angle")
    ax[1, 1].set(xlabel=r"$r/R$", ylabel=r"$z/z(0)$", title="small-V shape (dotted: SI Eq. 31)")
    ax[0, 2].set(xlabel=r"$r/\lambda$", ylabel=r"$z/\lambda$", title="largest V (dotted: Laplace cap 2/P)", aspect="equal")
    ax[1, 2].set(xlabel=r"$r/\lambda$", ylabel=r"$u=d\psi/ds$ ($1/\lambda$)", title=r"curvature; dashed: contact cond. $-\sqrt{2-g}$")
    for a_ in ax.flat:
        a_.legend(fontsize=7)
    fig.suptitle(f"Shape solver vs SI analytic limits; max |Hamiltonian (SI Eq. 19)| = {worstH:.1e}")
    fig.tight_layout()
    fig.savefig(args.output, dpi=130)


if __name__ == "__main__":
    main()
