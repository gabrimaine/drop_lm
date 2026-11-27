"""Data vs best-fit symmetric-lens model, with residuals."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--fit-profile", required=True)
    ap.add_argument("--fit-params", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    r, z = np.loadtxt(args.profile, comments=("#", "@"), usecols=(0, 1), unpack=True)
    m = np.genfromtxt(args.fit_profile, delimiter=",", names=True)
    d = json.loads(Path(args.fit_params).read_text())
    b, u = d["best_fit"], d["uncertainty"]["sigma"]
    zm = np.interp(r, m["r_nm"], m["z_nm"])
    res = z - zm

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
    ax.plot(r, z, "o", ms=3, color="k", label="data (td_reset2.xvg)")
    ax.plot(r, -z, "o", ms=3, color="0.6", label="mirror (assumed symmetric)")
    mono = m["is_bilayer"] == 0
    ax.plot(m["r_nm"][mono], m["z_nm"][mono], "C3-", lw=2, label="model: monolayer")
    ax.plot(m["r_nm"][mono], -m["z_nm"][mono], "C3-", lw=2, alpha=0.4)
    ax.plot(m["r_nm"][~mono], m["z_nm"][~mono], "C0-", lw=2, label="model: flat bilayer")
    ax.axvline(b["R_nm"], ls=":", color="0.5")
    txt = (f"$V$ = {b['V_nm3']:.0f} ± {u['V_nm3']:.0f} nm$^3$\n"
           f"$\\gamma_b/\\gamma_m$ = {b['gamma_ratio']:.3f} ± {u['gamma_ratio']:.3f}\n"
           f"$\\gamma_m$ = {b['gamma_m_mN_per_m']:.2f} ± {u['gamma_m_mN_per_m']:.2f} mN/m\n"
           f"$\\gamma_b$ = {b['gamma_b_mN_per_m']:.2f} ± {u['gamma_b_mN_per_m']:.2f} mN/m\n"
           f"$\\lambda_m$ = {b['lambda_nm']:.2f} nm, $R$ = {b['R_nm']:.2f} nm\n"
           f"$\\kappa_m$ = {d['assumptions']['kappa_m_kT']:g} $k_BT$ at {d['assumptions']['temperature_K']:g} K\n"
           f"rms = {d['fit_quality']['rms_nm']*1e3:.1f} pm  ({d['uncertainty']['method']})")
    ax.text(0.98, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(fc="white", ec="0.8"))
    ax.set(ylabel="z (nm)", aspect="equal", title="Symmetric LD lens: Deslandes, Thiam & Foret (2017) SI model")
    ax.legend(loc="lower right", fontsize=8)
    axr.plot(r, res * 1e3, "o-", ms=3, color="k")
    axr.axhline(0, color="C3")
    axr.axvline(b["R_nm"], ls=":", color="0.5")
    axr.set(xlabel="r (nm)", ylabel="data − model (pm)")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()
