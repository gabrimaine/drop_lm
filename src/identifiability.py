"""Profile chi^2 over (gamma_m, gamma_b/gamma_m): apex curvature (and z0) re-optimised at each node.

Shows whether one profile separates the two tensions (in the small-volume
limit only 2 gamma_m - gamma_b would be constrained, SI Eqs. 33-34).
"""
import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).parent))
from fitting import TEMPERATURES, LensModel, load_profile  # noqa: E402
from lens import KB  # noqa: E402

CTX = {}


def row(args):
    g, lams, theta_best = args
    model = LensModel(CTX["r"], CTX["z"], CTX["free"])
    out = []
    p = list(theta_best[2:])  # ln a (, z0)
    for lam in lams:
        def f(q):
            return model.residuals(np.r_[np.log(lam), g, q])
        try:
            model._last_x = None
            res = least_squares(f, p, diff_step=1e-6, xtol=1e-10, ftol=1e-10, max_nfev=60)
            if np.all(np.isfinite(res.fun)):
                p = list(res.x)
            out.append(float(res.fun @ res.fun))
        except Exception:
            out.append(np.nan)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--fit-params", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--temperature", required=True, choices=sorted(TEMPERATURES))
    ap.add_argument("--baseline", required=True, choices=["fixed", "free"])
    ap.add_argument("--fit-range", required=True, choices=["all", "r_le_40"])
    ap.add_argument("--kappa-kT", type=float, default=6.0)
    ap.add_argument("--n", type=int, default=21)
    ap.add_argument("--span", type=float, default=0.15, help="+- fractional span around best fit")
    ap.add_argument("--procs", type=int, default=16)
    args = ap.parse_args()

    d = json.loads(Path(args.fit_params).read_text())
    th = d["theta"]
    theta_best = np.array([th["ln_lambda_nm"], th["g"], th["ln_a"]] + ([th["z0_nm"]] if "z0_nm" in th else []))
    kappa = args.kappa_kT * KB * TEMPERATURES[args.temperature]
    r, z = load_profile(args.profile, args.fit_range)
    CTX.update(r=r, z=z, free=args.baseline == "free")

    gm0, g0 = d["best_fit"]["gamma_m_mN_per_m"], th["g"]
    gms = gm0 * np.linspace(1 - args.span, 1 + args.span, args.n)
    gs = np.clip(g0 * np.linspace(1 - args.span, 1 + args.span, args.n), 1e-3, 2 - 1e-3)
    lams = np.sqrt(kappa / (gms * 1e-3)) * 1e9
    with Pool(args.procs) as pool:
        S = np.array(pool.map(row, [(g, lams, theta_best) for g in gs]))
    N, k = len(z), len(theta_best)
    smin = np.nanmin(S)
    s2 = smin / (N - k)
    dchi2 = (S - smin) / s2

    fig, ax = plt.subplots(figsize=(7, 6))
    pc = ax.pcolormesh(gms, gs, np.log10(1 + dchi2), shading="auto", cmap="viridis_r")
    lev = [2.30, 11.83, 1e2, 1e3, 1e4]
    cs = ax.contour(gms, gs, dchi2, levels=lev, colors="w", linestyles=["-", ":", "--", "--", "--"])
    ax.clabel(cs, fmt={2.30: "68%", 11.83: "99.7%", 1e2: "Δχ²=1e2", 1e3: "1e3", 1e4: "1e4"}, fontsize=7)
    # lines of constant Delta gamma = 2 gamma_m - gamma_b (the only small-volume combination)
    G, Gm = np.meshgrid(gs, gms, indexing="ij")
    cd = ax.contour(gms, gs, (2 - G) * Gm, levels=8, colors="0.85", linewidths=0.6, linestyles="-.")
    ax.clabel(cd, fmt="Δγ=%.2f", fontsize=6)
    ax.plot(gm0, g0, "r*", ms=12, label="best fit")
    fig.colorbar(pc, label=r"$\log_{10}(1+\Delta\chi^2)$")
    ax.set(xlabel=r"$\gamma_m$ (mN/m)", ylabel=r"$\gamma_b/\gamma_m$",
           title=f"Profile Δχ² (apex curvature{' and z0' if CTX['free'] else ''} re-optimised)\n"
                 f"σ_z = {np.sqrt(s2)*1e3:.1f} pm from best fit")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(args.output, dpi=140)


if __name__ == "__main__":
    main()
