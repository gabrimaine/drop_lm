"""Refit the profile with the pre-fae3bc2 contour sampling and compare with fit_params.

The independent cross-check (cross_check) first disagreed with our fit by
0.3-0.4 sigma on the tensions. The cause was our sampling of the model contour
every 0.05 lambda (see old_sampling.py). This output reproduces the old fit
exactly as it was run then (same grid start, same optimiser, same Jacobian
errors), so the size of that bias is a number in the report, not a remark.

Writes one row per quantity: dense (fit_params, the current code), dense_std
(our quoted sigma), coarse (refit with the old sampling), coarse_std (its
Jacobian sigma), shift_sigma = (coarse - dense) / dense_std, unit. Two extra
rows give the interpolation error of the old sampling at the current best fit:
its maximum and rms over the data radii, in pm.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import lens as L  # noqa: E402
from fitting import TEMPERATURES, LensModel, delta_method, derived, fit, load_profile  # noqa: E402
from old_sampling import OLD_MAX_STEP, shoot_fixed_step, use_old_sampling  # noqa: E402

KEYS = [("V_nm3", "V", "Droplet volume", "nm^3"),
        ("gamma_ratio", "gamma_ratio", "Tension ratio gamma_b/gamma_m", ""),
        ("gamma_m_mN_per_m", "gamma_m", "Monolayer tension gamma_m", "mN/m"),
        ("gamma_b_mN_per_m", "gamma_b", "Bilayer tension gamma_b", "mN/m"),
        ("delta_gamma_mN_per_m", "delta_gamma", "Spreading tension 2 gamma_m - gamma_b", "mN/m"),
        ("lambda_nm", "lambda_m", "Monolayer bending length lambda_m", "nm"),
        ("R_nm", "R", "Contact-line radius R", "nm"),
        ("h_nm", "h", "Lens thickness h = 2 z(0)", "nm")]


def interpolation_error(theta, r):
    """Old-sampling minus exact model height at the data radii, at fixed theta (nm)."""
    lam, g, a = math.exp(theta["ln_lambda_nm"]), theta["g"], math.exp(theta["ln_a"])
    exact, _ = L.lens(a, g)
    coarse = shoot_fixed_step(a, exact.P, dense=True)
    return L.profile_on_grid(coarse, lam, r) - L.profile_on_grid(exact, lam, r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--fit-params", required=True)
    ap.add_argument("--kappa-kT", type=float, default=6.0)
    ap.add_argument("--temperature", required=True, choices=sorted(TEMPERATURES))
    ap.add_argument("--baseline", required=True, choices=["fixed", "free"])
    ap.add_argument("--fit-range", required=True, choices=["all", "r_le_40"])
    ap.add_argument("--max-step", type=float, default=OLD_MAX_STEP)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    kappa = args.kappa_kT * L.KB * TEMPERATURES[args.temperature]
    r, z = load_profile(args.profile, args.fit_range)
    ours = json.loads(Path(args.fit_params).read_text())
    err = interpolation_error(ours["theta"], r)  # before the swap: needs the exact lens

    use_old_sampling(args.max_step)
    model = LensModel(r, z, free_baseline=args.baseline == "free")
    res = fit(model)  # grid start and optimiser as in fit.py
    N, k = len(z), len(res.x)
    s2 = float(res.fun @ res.fun) / (N - k)
    cov = s2 * np.linalg.inv(res.jac.T @ res.jac)
    best, C = delta_method(model, res.x, cov, kappa, [key for key, *_ in KEYS])
    sd = np.sqrt(np.diag(C))

    b, s = ours["best_fit"], ours["uncertainty"]["sigma"]
    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "label", "dense", "dense_std", "coarse", "coarse_std",
                    "shift_sigma", "unit"])
        for i, (key, qid, label, unit) in enumerate(KEYS):
            w.writerow([qid, label, repr(float(b[key])), repr(float(s[key])), repr(float(best[key])),
                        repr(float(sd[i])), repr(float((best[key] - b[key]) / s[key])), unit])
        w.writerow(["rms", "RMS residual", repr(1e3 * ours["fit_quality"]["rms_nm"]), "",
                    repr(1e3 * float(np.sqrt(np.mean(res.fun ** 2)))), "", "", "pm"])
        w.writerow(["interp_max", "Largest interpolation error of the old sampling at our best fit",
                    "", "", repr(1e3 * float(np.abs(err).max())), "", "", "pm"])
        w.writerow(["interp_rms", "RMS interpolation error of the old sampling at our best fit",
                    "", "", repr(1e3 * float(np.sqrt(np.mean(err ** 2)))), "", "", "pm"])
        w.writerow(["max_step", "Old sampling step", "", "", repr(args.max_step), "", "", "lambda"])


if __name__ == "__main__":
    main()
