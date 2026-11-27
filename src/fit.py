"""Fit the Deslandes-Thiam-Foret symmetric lens model to a droplet profile.

Writes one JSON with best-fit physical parameters and their uncertainties.
"""
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from fitting import (TEMPERATURES, LensModel, delta_method, derived, fit,  # noqa: E402
                     load_profile)
from lens import KB  # noqa: E402

KEYS = ["V_nm3", "gamma_ratio", "gamma_m_mN_per_m", "gamma_b_mN_per_m",
        "delta_gamma_mN_per_m", "R_nm", "h_nm", "lambda_nm"]


def _refit(zb, r, theta, kappa, free_baseline):
    """One bootstrap refit, warm-started at the best fit; None if no lens is found."""
    mb = LensModel(r, zb, free_baseline=free_baseline)
    rb = fit(mb, theta0=theta)
    d = derived(mb, rb.x, kappa)
    return None if d is None else [d[kk] for kk in KEYS]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--kappa-kT", type=float, default=6.0)
    ap.add_argument("--temperature", required=True, choices=sorted(TEMPERATURES))
    ap.add_argument("--baseline", required=True, choices=["fixed", "free"])
    ap.add_argument("--fit-range", required=True, choices=["all", "r_le_40"])
    ap.add_argument("--uncertainty", required=True, choices=["jacobian", "residual_bootstrap"])
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--procs", type=int, default=1, help="worker processes for the bootstrap")
    args = ap.parse_args()

    T = TEMPERATURES[args.temperature]
    kappa = args.kappa_kT * KB * T
    r, z = load_profile(args.profile, args.fit_range)
    model = LensModel(r, z, free_baseline=args.baseline == "free")

    res = fit(model)
    theta = res.x
    resid = res.fun
    N, k = len(z), len(theta)
    s2 = float(resid @ resid) / (N - k)
    best = derived(model, theta, kappa)

    unc = {"method": args.uncertainty}
    if args.uncertainty == "jacobian":
        Jt = res.jac
        cov_theta = s2 * np.linalg.inv(Jt.T @ Jt)
        _, C = delta_method(model, theta, cov_theta, kappa, KEYS)
        sd = np.sqrt(np.diag(C))
        unc["sigma"] = dict(zip(KEYS, sd.tolist()))
        corr = C / np.outer(sd, sd)
        unc["correlation"] = {"keys": KEYS, "matrix": corr.round(4).tolist()}
        unc["note"] = ("Gaussian errors from the Jacobian with sigma_z^2 = SSR/(N-k); "
                       "residuals of a profile are correlated, so treat as a lower bound.")
    else:
        # Resamples are drawn up front, in the same order a serial loop would draw
        # them, so the result does not depend on --procs.
        rng = np.random.default_rng(args.seed)
        zfit = z + resid
        zbs = [zfit + rng.choice(resid, size=N, replace=True) for _ in range(args.n_boot)]
        job = partial(_refit, r=r, theta=theta, kappa=kappa, free_baseline=model.free_baseline)
        if args.procs > 1:
            with ProcessPoolExecutor(max_workers=args.procs) as ex:
                results = list(ex.map(job, zbs))
        else:
            results = [job(zb) for zb in zbs]
        S = np.array([row for row in results if row is not None])
        unc["n_boot_ok"] = len(S)
        unc["sigma"] = dict(zip(KEYS, S.std(axis=0, ddof=1).tolist()))
        unc["percentile_16_84"] = {kk: np.percentile(S[:, i], [16, 84]).tolist() for i, kk in enumerate(KEYS)}
        unc["correlation"] = {"keys": KEYS, "matrix": np.corrcoef(S.T).round(4).tolist()}
        unc["note"] = "iid residual bootstrap; ignores residual correlation along r."

    out = dict(
        model="Deslandes, Thiam & Foret, Biophys. J. 113:15 (2017), SI Sect. 5: symmetric lens",
        assumptions=dict(kappa_m_kT=args.kappa_kT, temperature_K=T, kappa_m_J=kappa,
                         profile="z(r) = upper oil/monolayer interface above a flat bilayer plane, up-down symmetric lens",
                         tau=0.0, C0=0.0),
        decisions=dict(temperature=args.temperature, baseline=args.baseline,
                       fit_range=args.fit_range, uncertainty=args.uncertainty),
        best_fit=best,
        theta=dict(ln_lambda_nm=float(theta[0]), g=float(theta[1]), ln_a=float(theta[2]),
                   **({"z0_nm": float(theta[3])} if model.free_baseline else {})),
        uncertainty=unc,
        fit_quality=dict(N=N, n_params=k, rms_nm=float(np.sqrt(np.mean(resid**2))),
                         max_abs_resid_nm=float(np.abs(resid).max()), sigma_z_nm=float(np.sqrt(s2)),
                         optimizer_status=int(res.status), optimizer_message=str(res.message),
                         n_model_evals=model.n_eval),
    )
    Path(args.output).write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
