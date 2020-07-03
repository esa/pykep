from math import acos, asin, cos, log, pi, sin

import numpy as np
from numpy.linalg import norm
from pykep import DAY2SEC, epoch, ic2par
from pykep.core import fb_prop, lambert_problem, propagate_lagrangian
from pykep.planet import jpl_lp
from pykep.trajopt import launchers


class _solar_orbiter_udp:
    def __init__(
        self, t0=[epoch(0), epoch(1000)], multi_objective=False, tof_encoding="direct"
    ):
        """
        Args:
            - multi_objective (``bool``): when True the problem fitness will return also the time of flight as an added objective
            - tof_encoding (``str``): one of 'direct', 'eta' or 'alpha'. Selects the encoding for the time of flights
            - tof (``list`` or ``list`` of ``list``): time of flight bounds. As documented in ``pykep.mga_1dsm``

        """

        vinf = [0.5, 2.5]
        multi_objective = False
        eta_lb = 0.1
        eta_ub = 0.9
        rp_ub = 30

        # Redefining the planets as to change their safe radius. 350km was given as safe distance.
        earth = jpl_lp("earth")
        earth.safe_radius = (earth.radius + 350000) / earth.radius
        venus = jpl_lp("venus")
        venus.safe_radius = (venus.radius + 350000) / venus.radius
        seq = [
            earth,
            venus,
            earth,
            earth,
            venus,
        ]  # alternative: Launch-Venus-Venus-Earth-Venus
        tof = [[10, 400]] * (len(seq) - 1)

        # Sanity checks
        # 1 - Planets need to have the same mu_central_body
        if [r.mu_central_body for r in seq].count(seq[0].mu_central_body) != len(seq):
            raise ValueError(
                "All planets in the sequence need to have identical mu_central_body"
            )
        # 2 - tof encoding needs to be one of 'alpha', 'eta', 'direct'
        if tof_encoding not in ["alpha", "eta", "direct"]:
            raise TypeError("tof encoding must be one of 'alpha', 'eta', 'direct'")
        # 3 - tof is expected to have different content depending on the tof_encoding
        if tof_encoding == "direct":
            if np.shape(np.array(tof)) != (len(seq) - 1, 2):
                raise TypeError(
                    "tof_encoding is "
                    + tof_encoding
                    + " and tof must be a list of two dimensional lists and with length equal to the number of legs"
                )
        if tof_encoding == "alpha":
            if np.shape(np.array(tof)) != (2,):
                raise TypeError(
                    "tof_encoding is "
                    + tof_encoding
                    + " and tof must be a list of two floats"
                )
        if tof_encoding == "eta":
            if np.shape(np.array(tof)) != ():
                raise TypeError(
                    "tof_encoding is " + tof_encoding + " and tof must be a float"
                )
        # 4 - Check launch window t0. If defined in terms of floats transform into epochs
        if len(t0) != 2:
            raise TypeError(
                "t0 is " + t0 + " while should be a list of two floats or epochs"
            )
        if type(t0[0]) is not epoch:
            t0[0] = epoch(t0[0])
        if type(t0[1]) is not epoch:
            t0[1] = epoch(t0[1])

        self._seq = seq
        self._t0 = t0
        self._tof = tof
        self._vinf = vinf
        self._tof_encoding = tof_encoding
        self._multi_objective = multi_objective
        self._eta_lb = eta_lb
        self._eta_ub = eta_ub
        self._rp_ub = rp_ub

        self.n_legs = len(seq) - 1
        self.common_mu = seq[0].mu_central_body

    def _decode_times_and_vinf(self, x):
        # 1 - we decode the times of flight
        if self._tof_encoding == "alpha":
            # decision vector is  [t0] + [u, v, Vinf, eta1, a1] + [beta, rp/rV, eta2, a2] + ... + [T]
            T = list([0] * (self.n_legs))
            for i in range(len(T)):
                T[i] = -log(x[5 + 4 * i])
            alpha_sum = sum(T)
            retval_T = [x[-1] * time / alpha_sum for time in T]
        elif self._tof_encoding == "direct":
            # decision vector is  [t0] + [u, v, Vinf, eta1, T1] + [beta, rp/rV, eta2, T2] + ...
            retval_T = x[5::4]
        elif self._tof_encoding == "eta":
            # decision vector is  [t0] + [u, v, Vinf, eta1, n1] + [beta, rp/rV, eta2, n2] + ...
            dt = self._tof
            T = [0] * self.n_legs
            for i in range(self.n_legs):
                T[i] = (dt - sum(T[:i])) * x[5 + 4 * i]
            retval_T = T

        # 2 - we decode the hyperbolic velocity at departure
        theta = 2 * pi * x[1]
        try:
            phi = acos(2 * x[2] - 1) - pi / 2
        except ValueError as e:
            print("x[2]:" + str(x[2]) + " is invalid value.")
            raise (e)

        Vinfx = x[3] * cos(phi) * cos(theta)
        Vinfy = x[3] * cos(phi) * sin(theta)
        Vinfz = x[3] * sin(phi)

        return (retval_T, Vinfx, Vinfy, Vinfz)

    def fitness(self, x):
        T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)
        # We transform it (only the needed component) to an equatorial system rotating along x
        # (this is an approximation, assuming vernal equinox is roughly x and the ecliptic plane is roughly xy)
        earth_axis_inclination = 0.409072975
        # This is different from the GTOP tanmEM problem, I think it was bugged there as the rotation was in the wrong direction.
        Vinfz = -Vinfy * sin(earth_axis_inclination) + Vinfz * cos(
            earth_axis_inclination
        )
        # And we find the vinf declination (in degrees)
        sindelta = Vinfz / x[3]
        declination = asin(sindelta) / np.pi * 180.0
        # We now have the initial mass of the spacecraft
        m_initial = launchers.atlas501(x[3] / 1000.0, declination)

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.n_legs + 1))
        r_P = list([None] * (self.n_legs + 1))
        v_P = list([None] * (self.n_legs + 1))
        DV = list([0.0] * (self.n_legs + 1))
        for i in range(len(self._seq)):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self._seq[i].eph(t_P[i])

        # 3 - We start with the first leg
        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        r, v = propagate_lagrangian(r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu)

        # Lambert arc to reach seq[1]
        dt = (1 - x[4]) * T[0] * DAY2SEC
        l = lambert_problem(r, r_P[1], dt, self.common_mu, cw=False, max_revs=0)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

        # 4 - And we proceed with each successive leg
        for i in range(1, self.n_legs):
            # Fly-by
            v_out = fb_prop(
                v_end_l,
                v_P[i],
                x[7 + (i - 1) * 4] * self._seq[i].radius,
                x[6 + (i - 1) * 4],
                self._seq[i].mu_self,
            )
            # s/c propagation before the DSM
            r, v = propagate_lagrangian(
                r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC, self.common_mu
            )
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[8 + (i - 1) * 4]) * T[i] * DAY2SEC
            l = lambert_problem(r, r_P[i + 1], dt, self.common_mu, cw=False, max_revs=0)
            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])

        # Last fly-by
        i = self.n_legs
        v_out = fb_prop(
            v_end_l,
            v_P[self.n_legs],
            x[7 + (i - 1) * 4] * self._seq[i].radius,
            x[6 + (i - 1) * 4],
            self._seq[i].mu_self,
        )

        a, e, i, W, w, E = ic2par(r_P[i], v_out, self.common_mu)

        if not self._multi_objective:
            return (-i, sum(DV) - 10, 209 - m_initial)
        else:
            return (-i, sum(T), sum(DV) - 10, 209 - m_initial)

    def get_nobj(self):
        return self._multi_objective + 1

    def get_bounds(self):
        t0 = self._t0
        tof = self._tof
        vinf = self._vinf
        seq = self._seq
        # Base for all possiblities (eta encoding)
        lb = (
            [t0[0].mjd2000]
            + [0.0, 0.0, vinf[0] * 1000, self._eta_lb, 1e-3]
            + [-2 * pi, np.nan, self._eta_lb, 1e-3] * (self.n_legs - 1)
        )
        ub = (
            [t0[1].mjd2000]
            + [1.0, 1.0, vinf[1] * 1000, self._eta_ub, 1.0 - 1e-3]
            + [2 * pi, self._rp_ub, self._eta_ub, 1.0 - 1e-3] * (self.n_legs - 1)
        )
        # Distinguishing among cases (only direct and alpha)
        if self._tof_encoding == "alpha":
            lb = lb + [tof[0]]
            ub = ub + [tof[1]]
        elif self._tof_encoding == "direct":
            for i in range(self.n_legs):
                lb[5 + 4 * i] = tof[i][0]
                ub[5 + 4 * i] = tof[i][1]

        # Setting the minimum rp/rP using the planet safe radius
        for i, pl in enumerate(seq[1:-1]):
            lb[7 + 4 * i] = pl.safe_radius / pl.radius

        # Setting final flyby
        lb = lb + [0, seq[-1].safe_radius / seq[-1].radius]
        ub = ub + [2 * pi, self._rp_ub]

        return (lb, ub)

    def get_nic(self):
        return 2
