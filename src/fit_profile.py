"""Best-fit model contour (full arclength parametrisation) as CSV."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from lens import lens  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-params", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--r-max", type=float, default=60.0, help="extend flat bilayer to this r [nm]")
    args = ap.parse_args()
    d = json.loads(Path(args.fit_params).read_text())
    th = d["theta"]
    lam, g, a = np.exp(th["ln_lambda_nm"]), th["g"], np.exp(th["ln_a"])
    z0 = th.get("z0_nm", 0.0)
    sh, _ = lens(a, g)
    s, r, z = sh.s * lam, sh.r * lam, sh.z * lam + z0
    psi = np.degrees(sh.psi)
    curv = sh.u / lam
    # flat bilayer beyond the contact line (z_b = 0, SI Sect. 5)
    rb = np.linspace(r[-1], args.r_max, 50)[1:]
    sb = s[-1] + (rb - r[-1])
    tab = np.vstack([
        np.column_stack([s, r, z, psi, curv, np.zeros_like(s)]),
        np.column_stack([sb, rb, np.full_like(rb, z0), np.zeros_like(rb), np.zeros_like(rb), np.ones_like(rb)]),
    ])
    np.savetxt(args.output, tab, delimiter=",", fmt="%.8g",
               header="s_nm,r_nm,z_nm,psi_deg,dpsi_ds_inv_nm,is_bilayer", comments="")


if __name__ == "__main__":
    main()
