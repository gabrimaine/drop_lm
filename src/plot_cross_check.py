"""Figure comparing our fit with the independent code (compare_fit.py).

Both best-fit shapes are rebuilt from the cross_check table: ours with
lens.lens at our (apex curvature, gamma_b/gamma_m), theirs with compare_fit's
own solve_shape at its (lambda, delta = 2 - gamma_b/gamma_m, reduced pressure),
seeded with its apex curvature. Nothing is refitted here.

Panels: (a) data and both models; (b) residuals of both; (c) |difference
between the two models| against the interpolation error of our pre-fae3bc2
sampling and the residual rms; (d) differences in units of our sigma, now and
with the old sampling (sampling_bias); (e) ratio of the error bars.
"""
import argparse
import csv
import importlib.util
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import lens as L  # noqa: E402
from fitting import load_profile  # noqa: E402
from old_sampling import shoot_fixed_step  # noqa: E402

QUANT = [("V", "$V$"), ("gamma_ratio", r"$\gamma_b/\gamma_m$"), ("gamma_m", r"$\gamma_m$"),
         ("gamma_b", r"$\gamma_b$"), ("delta_gamma", r"$\Delta\gamma$"),
         ("lambda_m", r"$\lambda_m$"), ("R", "$R$"), ("h", "$h$")]


def read_table(path):
    rows = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            rows[row["quantity"]] = {k: (float(v) if v not in ("", None) and k not in ("quantity", "label", "unit")
                                         else v) for k, v in row.items()}
    return rows


def num(x):
    return x if isinstance(x, float) else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--fit-params", required=True)
    ap.add_argument("--cross-check", required=True)
    ap.add_argument("--sampling-bias", required=True)
    ap.add_argument("--other-code", required=True)
    ap.add_argument("--fit-range", required=True, choices=["all", "r_le_40"])
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    spec = importlib.util.spec_from_file_location("compare_fit", args.other_code)
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)

    r, z = load_profile(args.profile, args.fit_range)
    ours = json.loads(Path(args.fit_params).read_text())
    cc = read_table(args.cross_check)
    sb = read_table(args.sampling_bias)
    th = ours["theta"]
    z0 = th.get("z0_nm", 0.0)
    lam, g, a = math.exp(th["ln_lambda_nm"]), th["g"], math.exp(th["ln_a"])

    sh, _ = L.lens(a, g)
    z_ours = L.profile_on_grid(sh, lam, r) + z0
    lam_t, ratio_t = cc["lambda_m"]["other"], cc["gamma_ratio"]["other"]
    prof = C.solve_shape(lam_t, 2.0 - ratio_t, cc["pressure"]["other"], u_seed=cc["apex_curvature"]["other"])
    z_other = C.model_z(prof, r) + z0
    old = shoot_fixed_step(a, sh.P, dense=True)
    z_old = L.profile_on_grid(old, lam, r) + z0

    pm = 1e3
    res_o, res_t = pm * (z - z_ours), pm * (z - z_other)
    rms = float(np.sqrt(np.mean(res_o ** 2)))
    R = sh.R * lam

    fig = plt.figure(figsize=(14, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(r, z, "o", ms=3, color="k", label="data")
    rr = np.linspace(0, r.max(), 800)
    ax.plot(rr, L.profile_on_grid(sh, lam, rr) + z0, "-", color="tab:red", lw=2.2, label="this analysis")
    ax.plot(rr, C.model_z(prof, rr) + z0, "--", color="tab:blue", lw=1.4, label="independent code")
    ax.axvline(R, ls=":", color="0.4")
    ax.set(xlabel="r (nm)", ylabel="z (nm)", title="(a) Both best fits on the data")
    ax.legend(frameon=False)

    ax = fig.add_subplot(gs[0, 1])
    ax.axhline(0, color="0.6", lw=0.8)
    ax.plot(r, res_o, "o", ms=4, mfc="none", color="tab:red", label=f"this analysis, rms {rms:.4f} pm")
    ax.plot(r, res_t, ".-", ms=3, lw=0.8, color="tab:blue",
            label=f"independent code, rms {float(np.sqrt(np.mean(res_t ** 2))):.4f} pm")
    ax.axvline(R, ls=":", color="0.4")
    ax.set(xlabel="r (nm)", ylabel="data - model (pm)", title="(b) Residuals")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[0, 2])
    inside = r < R
    d_now = np.abs(pm * (z_other - z_ours))
    d_old = np.abs(pm * (z_old - z_ours))
    ax.semilogy(r[inside], np.maximum(d_now[inside], 1e-7), ".-", color="tab:blue", lw=0.8,
                label=f"independent - ours, each at its best fit (max {d_now.max():.1e} pm)")
    ax.semilogy(r[inside], np.maximum(d_old[inside], 1e-7), ".-", color="tab:orange", lw=0.8,
                label=f"old 0.05 $\\lambda_m$ sampling - exact, same shape (max {d_old.max():.2f} pm)")
    ax.axhline(rms, ls="--", color="k", lw=1, label=f"residual rms ({rms:.2f} pm)")
    ax.axvline(R, ls=":", color="0.4")
    ax.set(xlabel="r (nm)", ylabel="|model difference| (pm)", ylim=(1e-5, 5e3),
           title="(c) How far apart the models are")
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax = fig.add_subplot(gs[1, 0:2])
    y = np.arange(len(QUANT))
    now = np.array([abs(num(cc[q]["diff_sigma"])) for q, _ in QUANT])
    then = np.array([abs(num(sb[q]["shift_sigma"])) for q, _ in QUANT])
    ax.barh(y + 0.2, then, 0.38, color="tab:orange", label="old sampling (before fae3bc2) - ours")
    ax.barh(y - 0.2, now, 0.38, color="tab:blue", label="independent code - ours")
    for yi, v in zip(y, then):
        ax.text(v * 1.15, yi + 0.2, f"{v:.2f}", va="center", fontsize=8)
    for yi, v in zip(y, now):
        ax.text(v * 1.15, yi - 0.2, f"{v:.4f}", va="center", fontsize=8)
    ax.axvline(1.0, color="k", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(1e-4, 60)
    ax.set_yticks(y, [lab for _, lab in QUANT])
    ax.invert_yaxis()
    ax.set(xlabel="|difference| / our 1σ", title="(d) Differences from our best fit, in units of our error")
    ax.legend(frameon=False, loc="lower right")

    ax = fig.add_subplot(gs[1, 2])
    rep = np.array([num(cc[q]["other_std"]) / num(cc[q]["ours_std"]) for q, _ in QUANT])
    full = np.array([num(cc[q]["other_std_full"]) / num(cc[q]["ours_std"]) for q, _ in QUANT])
    ax.barh(y - 0.2, rep, 0.38, color="0.65", label="as compare_fit.py reports")
    ax.barh(y + 0.2, full, 0.38, color="tab:blue", label="its covariance, all correlations")
    ax.axvline(1.0, color="k", lw=0.8, ls="--")
    ax.set_yticks(y, [lab for _, lab in QUANT])
    ax.invert_yaxis()
    ax.set_xlim(0.8, 1.3)
    ax.set(xlabel="independent σ / our σ", title="(e) Error bars")
    ax.legend(frameon=False, fontsize=8, loc="center right")

    fig.suptitle("This analysis against the independent implementation (data/raw/compare_fit.py)")
    fig.savefig(args.output, dpi=130)


if __name__ == "__main__":
    main()
