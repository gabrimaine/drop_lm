"""Refit the profile with an independent implementation and compare with fit_params.

The independent code is data/raw/compare_fit.py (its own shooting on the apex
curvature, its own parametrisation (lambda, delta = 2 - gamma_b/gamma_m, p_bar)
and its own starting guesses). Its fit is run from its five data-driven
"smart" starts, in parallel, with nothing taken from our fit. Its 125-point
grid search is skipped: on this profile all eleven of its starts, grid ones
included, converge to the same optimum.

compare_fit.py has no bilayer-plane offset. With a free baseline our fitted z0
is subtracted from the data first, so both codes fit the same shape.

Writes one row per quantity: ours, ours_std, other, other_std, diff_sigma =
(other - ours) / ours_std, unit.
"""
import argparse
import csv
import importlib.util
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from fitting import TEMPERATURES, load_profile  # noqa: E402
from lens import KB  # noqa: E402

C = None  # the independent module, loaded in each process by _load


def _load(path):
    global C
    spec = importlib.util.spec_from_file_location("compare_fit", path)
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)
    return C


def _one_start(args):
    path, r, z, th0 = args
    if C is None:
        _load(path)
    res = C.fit_profile(r, z, [(th0, None)], verbose=False)
    return float(res.cost), res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--fit-params", required=True)
    ap.add_argument("--other-code", required=True, help="path to compare_fit.py")
    ap.add_argument("--kappa-kT", type=float, default=6.0)
    ap.add_argument("--temperature", required=True, choices=sorted(TEMPERATURES))
    ap.add_argument("--baseline", required=True, choices=["fixed", "free"])
    ap.add_argument("--fit-range", required=True, choices=["all", "r_le_40"])
    ap.add_argument("--procs", type=int, default=5)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    _load(args.other_code)
    kappa = args.kappa_kT * KB * TEMPERATURES[args.temperature]
    r, z = load_profile(args.profile, args.fit_range)
    ours = json.loads(Path(args.fit_params).read_text())
    z0 = ours["best_fit"].get("z0_nm", 0.0) if args.baseline == "free" else 0.0
    zc = z - z0

    starts, _ = C.smart_starts(r, zc, kappa)
    jobs = [(args.other_code, r, zc, th) for th in starts]
    with ProcessPoolExecutor(max_workers=min(args.procs, len(jobs))) as ex:
        results = list(ex.map(_one_start, jobs))
    res = min(results, key=lambda cr: cr[0])[1]

    # Rebuild at the optimum as compare_fit.main does: lowest-SSR branch.
    lam, delta, pbar = C.theta_to_pars(res.x)
    prof = min(C.solve_shape_all(lam, delta, pbar),
               key=lambda p: float(np.sum((C.model_z(p, r) - zc) ** 2)))
    resid = zc - C.model_z(prof, r)
    out = C.derived_outputs(prof, kappa)
    sig = C.parameter_errors(res, r, zc, kappa, prof)

    b, s = ours["best_fit"], ours["uncertainty"]["sigma"]
    q = ours["fit_quality"]
    # quantity, label, unit, ours, ours_std, other, other_std
    rows = [
        ("V", "Droplet volume", "nm^3", b["V_nm3"], s["V_nm3"], out["V"], sig["sig_V"]),
        ("gamma_ratio", "Tension ratio gamma_b/gamma_m", "", b["gamma_ratio"], s["gamma_ratio"],
         out["ratio"], sig["sig_ratio"]),
        ("gamma_m", "Monolayer tension gamma_m", "mN/m", b["gamma_m_mN_per_m"], s["gamma_m_mN_per_m"],
         out["gamma_m"], sig["sig_gm"]),
        ("gamma_b", "Bilayer tension gamma_b", "mN/m", b["gamma_b_mN_per_m"], s["gamma_b_mN_per_m"],
         out["gamma_b"], sig["sig_gb"]),
        ("delta_gamma", "Spreading tension 2 gamma_m - gamma_b", "mN/m", b["delta_gamma_mN_per_m"],
         s["delta_gamma_mN_per_m"], out["dgamma"], None),
        ("lambda_m", "Monolayer bending length lambda_m", "nm", b["lambda_nm"], s["lambda_nm"],
         out["lam"], sig["sig_lam"]),
        ("R", "Contact-line radius R", "nm", b["R_nm"], s["R_nm"], out["R"], sig["sig_R"]),
        # compare_fit's h is the apex height z(0); ours is the thickness 2 z(0)
        ("h", "Lens thickness h = 2 z(0)", "nm", b["h_nm"], s["h_nm"], 2 * out["h"], 2 * sig["sig_h"]),
        ("psi_max", "Maximum monolayer slope angle", "deg", b["psi_max_deg"], None,
         math.degrees(out["psi_max"]), None),
        ("rms", "RMS residual", "pm", 1e3 * q["rms_nm"], None,
         1e3 * float(np.sqrt(np.mean(resid ** 2))), None),
        ("z0_subtracted", "Bilayer-plane offset taken from our fit", "pm", 1e3 * z0, None, 1e3 * z0, None),
    ]
    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "label", "ours", "ours_std", "other", "other_std", "diff_sigma", "unit"])
        for qid, label, unit, vo, so, vt, st in rows:
            d = (vt - vo) / so if so else None
            w.writerow([qid, label, repr(float(vo)), "" if so is None else repr(float(so)),
                        repr(float(vt)), "" if st is None else repr(float(st)),
                        "" if d is None else repr(float(d)), unit])


if __name__ == "__main__":
    main()
