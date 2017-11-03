from PyKEP.sims_flanagan import spacecraft
from PyKEP.core import MU_SUN, EARTH_VELOCITY, G0, AU
import numpy as np


class _dynamics(object):

    def __init__(self, sc=spacecraft(1000, 0.3, 2500), mu=MU_SUN, alpha=1, bound=True):

        # check spacecraft
        if isinstance(sc, spacecraft):
            self.spacecraft = sc
        else:
            raise TypeError(
                "Spacecraft should be instance of pontryagin.spacecraft class.")

        # check gravitational parametre
        if not (isinstance(mu, float) or isinstance(mu, int)):
            raise TypeError(
                "Gravitational parametre, mu, must be supplied as either int or float.")
        elif not mu > 0:
            raise TypeError(
                "Gravitational parametre, mu, must be a positive number.")
        else:
            self.mu = float(mu)

        # check homotopy
        if not (isinstance(alpha, float) or isinstance(alpha, int)):
            raise TypeError(
                "Homotopy parametre, alpha, must be supplied as float or int.")
        elif not (alpha >= 0 and alpha <= 1):
            raise ValueError(
                "Homotopy parametre, alpha, must be between 0 and 1.")
        else:
            self.alpha = float(alpha)

        # check bound
        if not isinstance(bound, bool):
            raise TypeError(
                "Control bounding parametre, bound, supplied as boolean.")
        else:
            self.bound = bool(bound)

        # check control validity
        if (self.alpha == 1 and self.bound == False):
            raise ValueError(
                "Control can only be unbounded with quadratic control, i.e. bound == True if alpha == 1.")
        else:
            pass

        # spacecraft parametres
        self.c1 = self.spacecraft.thrust
        self.c2 = self.spacecraft.thrust / (self.spacecraft.isp * G0)

        # define nondimenional units
        self.L = AU
        self.V = EARTH_VELOCITY
        self.M = self.spacecraft.mass
        self.A = (self.V * self.V) / self.L
        self.F = self.M * self.A
        self.T = self.L / self.V
        self.Q = self.F / self.V

        # nondimensionalise parametres
        self.c1 /= self.F
        self.c2 /= self.Q
        self.mu /= MU_SUN

    def _eom_fullstate(self, fullstate):

        # extract state and control
        x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lvm, obj = fullstate
        u, ix, iy, iz = self._pontryagin(fullstate)

        # common subexpression elimination
        x0 = self.c1 * u / m
        x1 = x**2
        x2 = y**2
        x3 = z**2
        x4 = x1 + x2 + x3
        x5 = self.mu / x4**(3 / 2)
        x6 = x4**(-5 / 2)
        x7 = 3 * lvy * self.mu * x6 * y
        x8 = 3 * lvz * self.mu * x6 * z
        x9 = -x5
        x10 = 3 * self.mu * x6
        x11 = 3 * lvx * self.mu * x * x6
        x12 = self.c1 * u / m**2

        # fullstate transition
        dfs = np.array([
            vx,
            vy,
            vz,
            ix * x0 - x * x5,
            iy * x0 - x5 * y,
            iz * x0 - x5 * z,
            -self.c2 * u,
            -lvx * (x1 * x10 + x9) - x * x7 - x * x8,
            -lvy * (x10 * x2 + x9) - x11 * y - x8 * y,
            -lvz * (x10 * x3 + x9) - x11 * z - x7 * z,
            -lx,
            -ly,
            -lz,
            ix * lvx * x12 + iy * lvy * x12 + iz * lvz * x12,
            self.alpha*u + (1 - self.alpha)*u**2
        ])

        return dfs

    def _eom_fullstate_jac(self, fullstate):

        # extract state and control
        x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lvm = fullstate
        u, ix, iy, iz = self._pontryagin(fullstate)

        # common subexpression elimination
        x0 = x**2
        x1 = y**2
        x2 = z**2
        x3 = x0 + x1 + x2
        x4 = self.mu / x3**(3 / 2)
        x5 = -x4
        x6 = x3**(-5 / 2)
        x7 = 3 * self.mu * x6
        x8 = x0 * x7
        x9 = self.mu * x * x6
        x10 = 3 * x9
        x11 = x10 * y
        x12 = x10 * z
        x13 = self.c1 * u / m**2
        x14 = ix * x13
        x15 = x1 * x7
        x16 = self.mu * x6 * y
        x17 = 3 * x16
        x18 = x17 * z
        x19 = iy * x13
        x20 = x2 * x7
        x21 = iz * x13
        x22 = -lvy * x17
        x23 = self.mu * x6 * z
        x24 = 3 * x23
        x25 = -lvz * x24
        x26 = x3**(-7 / 2)
        x27 = 15 * self.mu * x0 * x26
        x28 = x27 * y
        x29 = x27 * z
        x30 = 15 * self.mu * x26
        x31 = 15 * self.mu * x * x26 * y * z
        x32 = lvz * x31
        x33 = 15 * self.mu * x * x26
        x34 = x1 * x33
        x35 = lvy * x31
        x36 = x2 * x33
        x37 = -x11
        x38 = -x12
        x39 = -lvx * x10
        x40 = x1 * x30 * z
        x41 = lvx * x31
        x42 = x2 * x30 * y
        x43 = -x18
        x44 = 2 * self.c1 * u / m**3

        # fullstate transition jacobian
        dfsj = np.array([
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [x5 + x8, x11, x12, 0, 0, 0, -x14, 0, 0, 0, 0, 0, 0, 0],
            [x11, x15 + x5, x18, 0, 0, 0, -x19, 0, 0, 0, 0, 0, 0, 0],
            [x12, x18, x20 + x5, 0, 0, 0, -x21, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [-lvx * (-x**3 * x30 + 9 * x9) + lvy * x28 + lvz * x29 + x22 + x25, -lvx * (x17 - x28) - lvy * x10 + lvy *
             x34 + x32, -lvx * (x24 - x29) - lvz * x10 + lvz * x36 + x35, 0, 0, 0, 0, 0, 0, 0, x4 - x8, x37, x38, 0],
            [-lvx * x17 + lvx * x28 - lvy * (x10 - x34) + x32, lvx * x34 - lvy * (9 * x16 - x30 * y**3) + lvz * x40 + x25 + x39, -lvy * (
                x24 - x40) - lvz * x17 + lvz * x42 + x41, 0, 0, 0, 0, 0, 0, 0, x37, -x15 + x4, x43, 0],
            [-lvx * x24 + lvx * x29 - lvz * (x10 - x36) + x35, -lvy * x24 + lvy * x40 - lvz * (x17 - x42) + x41, lvx *
             x36 + lvy * x42 - lvz * (9 * x23 - x30 * z**3) + x22 + x39, 0, 0, 0, 0, 0, 0, 0, x38, x43, -x20 + x4, 0],
            [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, -lvx * ix * x44 - lvy * iy *
                x44 - lvz * iz * x44, 0, 0, 0, x14, x19, x21, 0]
        ])

        return dfsj

    def _hamiltonian(self, fullstate):

        # extract fullstate and control
        x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm, obj = fullstate
        u, ix, iy, iz = self._pontryagin(fullstate)

        # common subexpression elimination
        x0 = self.c1 * u / m
        x1 = self.mu / (x**2 + y**2 + z**2)**(3 / 2)

        # Hamiltonian
        H = -lm * self.c2 * u + lx * vx + ly * vy + lz * vz + lvx * (ix * x0 - x * x1) + lvy * (
            iy * x0 - x1 * y) + lvz * (iz * x0 - x1 * z) + self.alpha * u + u**2 * (-self.alpha + 1)

        return H

    def _pontryagin(self, fullstate):

        # extract fullstate
        x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm, obj = fullstate

        # magnitude of lv
        lv = (lvx**2 + lvy**2 + lvz**2)**0.5

        # if mass-optimal control
        if self.alpha == 1:

            # switching function
            s = 1 - self.c1 * lv / m - self.c2 * lm

            # bang-bang control
            if s >= 0:
                u = 0.

            elif s < 0:
                u = 1.

        # if quadratic control
        else:
            u = (self.c1 * lv + m * (self.c2 * lm - self.alpha)) / \
                (2 * m * (1 - self.alpha))

            # if throttle is bounded
            if self.bound:
                u = max(u, 0.0)
                u = min(u, 1.0)

            elif not self.bound:
                pass

        # throttle direction
        ix = -lvx / lv
        iy = -lvy / lv
        iz = -lvz / lv

        # assemble control decision
        control = np.array([u, ix, iy, iz])

        return control
