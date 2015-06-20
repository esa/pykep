# -*- coding: iso-8859-1 -*-
"""
This module contains all the classes that allow to construct efficiently
low-thrust tajectories using our own flavour of the Sims-Flanagan model: a trajectory
transcription method that forms the basis for MALTO, the software in use in JPL
for preliminary interplanetary trajectory design.
"""
from PyKEP.sims_flanagan._sims_flanagan import *


def _leg_get_states(self):
    """
    Returns the spacecraft states (t,r,v,m) at the leg grid points

    Examples::

      times,r,v,m = PyKEP.sims_flanagan.leg.get_states()
    """
    from PyKEP import propagate_lagrangian, AU, DAY2SEC, G0, propagate_taylor
    import numpy as np
    from scipy.linalg import norm
    from math import exp

    # We compute the number of segments for forward and backward propagation
    n_seg = len(self.get_throttles())
    fwd_seg = (n_seg + 1) / 2
    back_seg = n_seg / 2

    # We extract information on the spacecraft
    sc = self.get_spacecraft()
    isp = sc.isp
    max_thrust = sc.thrust

    # And on the leg
    throttles = self.get_throttles()
    mu = self.get_mu()

    # time grid
    t_grid = [0.0] * (n_seg * 2 + 2)

    # Forward propagation

    # x,y,z contain the cartesian components of all points (grid+midpints)
    x = [0.0] * (fwd_seg * 2 + 1)
    y = [0.0] * (fwd_seg * 2 + 1)
    z = [0.0] * (fwd_seg * 2 + 1)
    vx = [0.0] * (fwd_seg * 2 + 1)
    vy = [0.0] * (fwd_seg * 2 + 1)
    vz = [0.0] * (fwd_seg * 2 + 1)
    mass = [0.0] * (fwd_seg * 2 + 1)

    state = self.get_xi()

    # Initial conditions
    r = state.r
    v = state.v
    m = state.m
    x[0], y[0], z[0] = r
    vx[0], vy[0], vz[0] = v
    mass[0] = m

    # We compute all points by propagation
    for i, t in enumerate(throttles[:fwd_seg]):
        t_grid[2 * i] = t.start.mjd2000
        t_grid[2 * i + 1] = t.start.mjd2000 + \
            (t.end.mjd2000 - t.start.mjd2000) / 2.
        dt = (t.end.mjd - t.start.mjd) * DAY2SEC
        alpha = min(norm(t.value), 1.0)
        # Keplerian propagation and dV application
        if self.high_fidelity is False:
            dV = [max_thrust / m * dt * dumb for dumb in t.value]

            r, v = propagate_lagrangian(r, v, dt / 2, mu)
            x[2 * i + 1], y[2 * i + 1], z[2 * i + 1] = r
            vx[2 * i + 1], vy[2 * i + 1], vz[2 * i + 1] = v
            mass[2 * i + 1] = m
            # v= v+dV
            v = [a + b for a, b in zip(v, dV)]
            r, v = propagate_lagrangian(r, v, dt / 2, mu)
            m *= exp(-norm(dV) / isp / G0)

            x[2 * i + 2], y[2 * i + 2], z[2 * i + 2] = r
            vx[2 * i + 2], vy[2 * i + 2], vz[2 * i + 2] = v
            mass[2 * i + 2] = m

        # Taylor propagation of constant thrust u
        else:
            u = [max_thrust * dumb for dumb in t.value]

            r, v, m = propagate_taylor(
                r, v, m, u, dt / 2, mu, isp * G0, -10, -10)
            x[2 * i + 1], y[2 * i + 1], z[2 * i + 1] = r
            vx[2 * i + 1], vy[2 * i + 1], vz[2 * i + 1] = v
            mass[2 * i + 1] = m

            r, v, m = propagate_taylor(
                r, v, m, u, dt / 2, mu, isp * G0, -10, -10)
            x[2 * i + 2], y[2 * i + 2], z[2 * i + 2] = r
            vx[2 * i + 2], vy[2 * i + 2], vz[2 * i + 2] = v
            mass[2 * i + 2] = m

    t_grid[2 * i + 2] = t.end.mjd2000

    # Backward propagation

    # x,y,z will contain the cartesian components of
    x_back = [0.123] * (back_seg * 2 + 1)
    y_back = [0.123] * (back_seg * 2 + 1)
    z_back = [0.123] * (back_seg * 2 + 1)
    vx_back = [0.0] * (back_seg * 2 + 1)
    vy_back = [0.0] * (back_seg * 2 + 1)
    vz_back = [0.0] * (back_seg * 2 + 1)
    mass_back = [0.0] * (back_seg * 2 + 1)

    state = self.get_xf()

    # Final conditions
    r = state.r
    v = state.v
    m = state.m
    x_back[-1], y_back[-1], z_back[-1] = r
    vx_back[-1], vy_back[-1], vz_back[-1] = v
    mass_back[-1] = m

    for i, t in enumerate(throttles[-1:-back_seg - 1:-1]):
        t_grid[-2 * i - 2] = t.end.mjd2000 - \
            (t.end  .mjd2000 - t.start.mjd2000) / 2.
        t_grid[-2 * i - 1] = t.end.mjd2000
        dt = (t.end.mjd - t.start.mjd) * DAY2SEC
        alpha = min(norm(t.value), 1.0)
        if self.high_fidelity is False:
            dV = [max_thrust / m * dt * dumb for dumb in t.value]
            r, v = propagate_lagrangian(r, v, -dt / 2, mu)
            x_back[-2 * i - 2], y_back[-2 * i - 2], z_back[-2 * i - 2] = r
            vx_back[-2 * i - 2], vy_back[-2 * i - 2], vz_back[-2 * i - 2] = v
            mass_back[-2 * i - 2] = m
            # v= v+dV
            v = [a - b for a, b in zip(v, dV)]
            r, v = propagate_lagrangian(r, v, -dt / 2, mu)
            m *= exp(norm(dV) / isp / G0)

            x_back[-2 * i - 3], y_back[-2 * i - 3], z_back[-2 * i - 3] = r
            vx_back[-2 * i - 3], vy_back[-2 * i - 3], vz_back[-2 * i - 3] = v
            mass_back[-2 * i - 3] = m

        else:
            u = [max_thrust * dumb for dumb in t.value]
            r, v, m = propagate_taylor(
                r, v, m, u, -dt / 2, mu, isp * G0, -10, -10)
            x_back[-2 * i - 2], y_back[-2 * i - 2], z_back[-2 * i - 2] = r
            vx_back[-2 * i - 2], vy_back[-2 * i - 2], vz_back[-2 * i - 2] = v
            mass_back[-2 * i - 2] = m

            r, v, m = propagate_taylor(
                r, v, m, u, -dt / 2, mu, isp * G0, -10, -10)
            x_back[-2 * i - 3], y_back[-2 * i - 3], z_back[-2 * i - 3] = r
            vx_back[-2 * i - 3], vy_back[-2 * i - 3], vz_back[-2 * i - 3] = v
            mass_back[-2 * i - 3] = m

    t_grid[-2 * i - 3] = t.start.mjd2000
    x = x + x_back
    y = y + y_back
    z = z + z_back
    vx = vx + vx_back
    vy = vy + vy_back
    vz = vz + vz_back
    mass = mass + mass_back

    return t_grid, list(zip(x, y, z)), list(zip(vx, vy, vz)), mass

leg.get_states = _leg_get_states


def _leg_eph(self, t):
    """
    Computes the ephemerides (r, v) along the leg.  Should only be called on high_fidelity legs having
    no state mismatch. Otherwise the values returned will not correspond to physical quantities.

     - t: epoch (either a pykep epoch, or assumes mjd2000). This value must be between self.get_ti() and self.get_tf()
    """
    from PyKEP import epoch, propagate_taylor, G0, DAY2SEC
    from bisect import bisect

    if isinstance(t, epoch):
        t0 = t.mjd2000
    else:
        t0 = t

    # We check that requested epoch is valid
    if (t0 < self.get_ti().mjd2000) or (t0 > self.get_tf().mjd2000):
        raise ValueError("The requested epoch is out of bounds")

    mu = self.get_mu()
    sc = self.get_spacecraft()
    isp = sc.isp
    max_thrust = sc.thrust
    t_grid, r, v, m = self.get_states()

    # We move the starting epoch to 0 for convenience
    T = max(0, t0 - self.get_ti().mjd2000)
    t_grid = [it - self.get_ti().mjd2000 for it in t_grid]
    t_grid[0] = 0.

    # If by chance its in the grid node, we are done
    if T in t_grid:
        idx = t_grid.index(T)
        return r[idx], v[idx], m[idx]

    # Otherwise we repropagate
    idx = bisect(t_grid, T) - 1
    r0 = r[idx]
    v0 = v[idx]
    m0 = m[idx]
    midt = t_grid[len(t_grid) / 2]

    if T < midt:
        idx_thrust = idx / 2
    else:
        idx_thrust = idx / 2 - 1

    dt_int = (T - t_grid[idx]) * DAY2SEC
    th = self.get_throttles()[idx_thrust].value
    dd
    return propagate_taylor(r0, v0, m0, [d * max_thrust for d in th], dt_int, mu, isp * G0, -12, -12)

leg.eph = _leg_eph
