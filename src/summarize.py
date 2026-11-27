"""Flatten fit_params.json into a long table the report can read value by value.

One row per quantity: quantity, label, value, value_std, unit. value_std is
empty where the fit gives no uncertainty for that quantity.
"""
import argparse
import csv
import json
from pathlib import Path

# (quantity id, label, unit, where to read the value, scale to display unit)
ROWS = [
    ("V", "Droplet volume", "nm^3", ("best_fit", "V_nm3"), 1.0),
    ("gamma_ratio", "Tension ratio gamma_b/gamma_m", "", ("best_fit", "gamma_ratio"), 1.0),
    ("gamma_m", "Monolayer tension gamma_m", "mN/m", ("best_fit", "gamma_m_mN_per_m"), 1.0),
    ("gamma_b", "Bilayer tension gamma_b", "mN/m", ("best_fit", "gamma_b_mN_per_m"), 1.0),
    ("delta_gamma", "Spreading tension 2 gamma_m - gamma_b", "mN/m", ("best_fit", "delta_gamma_mN_per_m"), 1.0),
    ("lambda_m", "Monolayer bending length lambda_m", "nm", ("best_fit", "lambda_nm"), 1.0),
    ("R", "Contact-line radius R", "nm", ("best_fit", "R_nm"), 1.0),
    ("h", "Lens thickness h = 2 z(0)", "nm", ("best_fit", "h_nm"), 1.0),
    ("R_LD", "Equal-volume sphere radius", "nm", ("best_fit", "R_LD_nm"), 1.0),
    ("V_over_lambda3", "Reduced volume V/lambda_m^3", "", ("best_fit", "V_over_lambda3"), 1.0),
    ("p_bar", "Reduced pressure p lambda_m/gamma_m", "", ("best_fit", "p_bar"), 1.0),
    ("psi_max", "Maximum monolayer slope angle", "deg", ("best_fit", "psi_max_deg"), 1.0),
    ("psi_large_volume", "Large-volume contact angle arccos(gamma_b/2gamma_m)", "deg",
     ("best_fit", "large_volume_contact_angle_deg"), 1.0),
    ("z0", "Bilayer-plane offset z0", "pm", ("best_fit", "z0_nm"), 1e3),
    ("rms", "RMS residual", "pm", ("fit_quality", "rms_nm"), 1e3),
    ("max_resid", "Largest absolute residual", "pm", ("fit_quality", "max_abs_resid_nm"), 1e3),
    ("n_points", "Profile points fitted", "", ("fit_quality", "N"), 1.0),
    ("temperature", "Temperature for kappa_m", "K", ("assumptions", "temperature_K"), 1.0),
    ("kappa_m", "Monolayer bending rigidity kappa_m", "kBT", ("assumptions", "kappa_m_kT"), 1.0),
]

# best_fit key -> key in uncertainty.sigma (only these carry an uncertainty)
SIGMA = {"V_nm3", "gamma_ratio", "gamma_m_mN_per_m", "gamma_b_mN_per_m",
         "delta_gamma_mN_per_m", "R_nm", "h_nm", "lambda_nm"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-params", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    d = json.loads(Path(args.fit_params).read_text())
    sigma = d["uncertainty"]["sigma"]
    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "label", "value", "value_std", "unit"])
        for qid, label, unit, (sect, key), scale in ROWS:
            if key not in d[sect]:
                continue  # e.g. z0 exists only with a free baseline
            v = d[sect][key] * scale
            s = sigma[key] * scale if sect == "best_fit" and key in SIGMA else ""
            w.writerow([qid, label, repr(float(v)), "" if s == "" else repr(float(s)), unit])


if __name__ == "__main__":
    main()
