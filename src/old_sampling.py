"""The contour sampling used before fae3bc2, kept to measure what it cost.

Until fae3bc2, lens.shoot(dense=True) capped the DOP853 step at 0.05 lambda and
returned the integrator's own steps; profile_on_grid then interpolated z(r)
linearly between them. The integration itself was as accurate as now (same
rtol/atol); only the sampling was coarse. shoot_fixed_step reproduces that
function line for line, and use_old_sampling() swaps it into lens.lens so that
fitting.fit runs exactly as it did then.
"""
import numpy as np
from scipy.integrate import solve_ivp

import lens as L

OLD_MAX_STEP = 0.05  # lambda


def shoot_fixed_step(a, P, s_max=200.0, dense=False, rtol=1e-10, atol=1e-12,
                     max_step=OLD_MAX_STEP, **_):
    """lens.shoot as it was before fae3bc2 (git show fae3bc2^:src/lens.py)."""

    def contact(s, y, P):
        return y[2]

    contact.terminal = True
    contact.direction = -1

    def overturn(s, y, P):
        return y[2] - np.pi

    overturn.terminal = True

    def axis_hit(s, y, P):
        return y[0] - 1e-6

    axis_hit.terminal = True
    axis_hit.direction = -1

    y0 = L.axis_start(a, P)
    sol = solve_ivp(L.rhs, (y0[0], s_max), y0, args=(P,), method="DOP853",
                    rtol=rtol, atol=atol, events=(contact, overturn, axis_hit),
                    dense_output=False, max_step=max_step if dense else np.inf)
    ok = sol.status == 1 and len(sol.t_events[0]) > 0
    s, Y = sol.t, sol.y
    if ok:
        s = np.append(s[:-1], sol.t_events[0][0])
        Y = np.column_stack([Y[:, :-1], sol.y_events[0][0]])
    z = Y[1] - Y[1, -1]
    return L.Shape(a, P, ok, s, Y[0], z, Y[2], Y[3], Y[4])


def use_old_sampling(max_step=OLD_MAX_STEP):
    """Make lens.lens (and so fitting.fit) use the old sampler in this process."""
    def shoot(a, P, s_max=200.0, dense=False, rtol=1e-10, atol=1e-12, **_):
        return shoot_fixed_step(a, P, s_max=s_max, dense=dense, rtol=rtol, atol=atol,
                                max_step=max_step)
    L.shoot = shoot
