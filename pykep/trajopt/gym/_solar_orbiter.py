from math import pi, sin, sqrt
from typing import Any, List, Tuple

import numpy as np

from pykep import AU, DAY2SEC, RAD2DEG, SEC2DAY, epoch, ic2par
from pykep.core import fb_prop, fb_vel, lambert_problem, propagate_lagrangian
from pykep.planet import jpl_lp
from pykep.trajopt import launchers
from pykep.trajopt._lambert import lambert_problem_multirev


class _solar_orbiter_udp:

    earth = jpl_lp("earth")
    venus = jpl_lp("venus")

    def __init__(
        self,
        t0=[epoch(0), epoch(1000)],
        multi_objective=False,
        tof_encoding="direct",
        max_revs=0,
        dummy_DSMs=False,
        seq = [
            earth,
            venus,
            venus,
            earth,
            venus,
            venus,
            venus,
            venus,
            venus,
        ],
    ) -> None:
        """
        Args:
            - multi_objective (``bool``): when True the problem fitness will return also the time of flight as an added objective
            - tof_encoding (``str``): one of 'direct', 'eta' or 'alpha'. Selects the encoding for the time of flights
            - tof (``list`` or ``list`` of ``list``): time of flight bounds. As documented in ``pykep.mga_1dsm``
            - max_revs (``int``): maximal number of revolutions for lambert transfer

        """

        multi_objective = False
        eta_lb = 0.1
        eta_ub = 0.9
        rp_ub = 30
        safe_distance = 350000
        min_start_mass = 1800

        # Redefining the planets as to change their safe radius. 350km was given as safe distance.

        # alternative: Launch-Venus-Venus-Earth-Venus
        for i in range(len(seq)):
            seq[i].safe_radius = (seq[i].radius + safe_distance) / seq[i].radius
        
        tof = [[10, 600]] * (len(seq) - 1)

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

        # 5 - Resonant Flybys mess up the lambert arcs. Watch out for planet sequences that could have them:
        possibly_resonant: List[bool] = [
            dummy_DSMs and seq[i] == seq[i + 1] for i in range(len(seq) - 1)
        ]
        if possibly_resonant[0]:
            raise NotImplementedError("Resonant flybys not yet supported at launch.")

        self._seq = seq
        self._t0 = t0
        self._tof = tof
        self._tof_encoding = tof_encoding
        self._multi_objective = multi_objective
        self.max_revs = max_revs
        self._eta_lb = eta_lb
        self._eta_ub = eta_ub
        self._rp_ub = rp_ub
        self._safe_distance = safe_distance
        self._min_start_mass = min_start_mass
        self._dummy_DSM = possibly_resonant

        self._n_legs = len(seq) - 1
        self._common_mu = seq[0].mu_central_body

    def _decode_tofs(self, x: List[float]) -> List[float]:
        tail = 2 + 3 * sum(self._dummy_DSM)
        if self._tof_encoding == "alpha":
            # decision vector is  [t0, T, a1, a2, ....]
            T = np.log(x[2:-tail])
            return T / sum(T) * x[1]
        elif self._tof_encoding == "direct":
            # decision vector is  [t0, T1, T2, T3, ... ]
            return x[1:-tail]
        elif self._tof_encoding == "eta":
            # decision vector is  [t0, n1, n2, n3, ... ]
            dt = self._tof
            T = [0] * self._n_legs
            T[0] = dt * x[1]
            for i in range(1, len(T)):
                T[i] = (dt - sum(T[:i])) * x[i + 1]
            return T
        else:
            raise TypeError("tof encoding must be one of 'alpha', 'eta', 'direct'")


    def _compute_dvs(self, x: List[float]) -> Tuple[List[float], List[Any], List[float], List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]], List[float]]:

        # 1 -  we 'decode' the times of flights and compute epochs (mjd2000)
        T = self._decode_tofs(x)  # [T1, T2 ...]
        ep = np.insert(T, 0, x[0])  # [t0, T1, T2 ...]
        ep = np.cumsum(ep)  # [t0, t1, t2, ...]

        if not len(ep) == len(self._seq):
            print(len(ep), len(self._seq))
        # 2 - we compute the ephemerides
        r = [(0,0,0)] * len(self._seq)
        v = [(0,0,0)] * len(self._seq)
        for i in range(len(self._seq)):
            r[i], v[i] = self._seq[i].eph(ep[i])

        rf_index = np.cumsum(self._dummy_DSM)
        # 3 - we solve the lambert problems
        l = list()
        DVdsm = list()
        ballistic_legs: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = list()
        ballistic_ep: List[float] = list()

        vi = v[0]
        for i in range(self._n_legs):
            # start leg at planet i, before flyby
            ri = r[i]
            Ti = T[i]

            if np.any(np.isnan(vi)):
                return ([np.nan], [], ep, [], [])

            if self._dummy_DSM[i]:
                # insert dummy DSM
                flyby_param_index = (
                    len(x) - 3 * sum(self._dummy_DSM) + (rf_index[i] - 1) * 3
                )
                tof_ratio, beta, r_p = x[
                    flyby_param_index : flyby_param_index + 3
                ]  # TODO: adapt this to different encodings
                if tof_ratio < 0 or tof_ratio > 1:
                    raise ValueError(
                        "Time of flight ratio of " + str(tof_ratio) + " is invalid."
                    )

                # perform unpowered flyby
                vi = fb_prop(
                    vi, ri, r_p * self._seq[i].radius, beta, self._seq[i].mu_self
                )  # TODO: is seq[i] the correct planet?
                ballistic_legs.append((ri, vi))
                ballistic_ep.append(ep[i])

                if np.any(np.isnan(vi)):
                    return ([np.nan], [], ep, [], [])

                if tof_ratio > 0:
                    # propagate after flyby
                    ri, vi = propagate_lagrangian(
                        ri, vi, T[i] * DAY2SEC * tof_ratio, self._common_mu
                    )

                # adapt the remaining time for the lambert leg
                Ti = Ti * (1 - tof_ratio)
                if Ti == 0:
                    Ti = SEC2DAY

            # call lambert solver on remaining leg - either after flyby or DSM
            lp = lambert_problem_multirev(
                vi,
                lambert_problem(
                    ri, r[i + 1], Ti * DAY2SEC, self._common_mu, False, self.max_revs,
                ),
            )
            l.append(lp)

            # add delta v of DSM
            if self._dummy_DSM[i]:
                DVdsm.append(
                    np.linalg.norm([a - b for a, b in zip(vi, lp.get_v1()[0])])
                )

            vi = lp.get_v2()[0]
            ep_start_lambert = ep[i+1] - Ti
            # ri is now either the position of planet i or the position of the DSM
            ballistic_legs.append((ri, lp.get_v1()[0]))
            ballistic_ep.append(ep_start_lambert)

        # TODO: add ballistic leg of final flyby

        # 4 - we compute the various dVs needed at fly-bys to match incoming
        # and outcoming, also delta v of deep space maneuvers
        assert len(DVdsm) == sum(self._dummy_DSM)
        DV = list()
        j = 0  # index of DSM
        for i in range(len(l) - 1):
            if self._dummy_DSM[i + 1]:
                # use delta v of deep space maneuver
                DV.append(DVdsm[j])
                j += 1
            else:
                # use delta v of flyby
                vin = [a - b for a, b in zip(l[i].get_v2()[0], v[i + 1])]
                vout = [a - b for a, b in zip(l[i + 1].get_v1()[0], v[i + 1])]
                DV.append(fb_vel(vin, vout, self._seq[i + 1]))
        assert(j == len(DVdsm))
        assert(len(ballistic_legs) == sum(self._dummy_DSM) + len(l))
        assert(len(ballistic_ep) == len(ballistic_legs))
        return (DV, l, ep, ballistic_legs, ballistic_ep)

    # Objective function
    def fitness(self, x):
        if len(x) != len(self.get_bounds()[0]):
            raise ValueError(
                "Expected "
                + str(len(self.get_bounds()[0]))
                + " parameters but got "
                + str(len(x))
            )

        DV, lamberts, ep, _, _ = self._compute_dvs(x)
        T = self._decode_tofs(x)
        rf_index = np.cumsum(self._dummy_DSM)

        if np.any(np.isnan(DV)):
            return [np.inf] + [T]*self._multi_objective + [np.nan]*4

        # compute launch velocity and declination
        Vinfx, Vinfy, Vinfz = [
            a - b for a, b in zip(lamberts[0].get_v1()[0], self._seq[0].eph(ep[0])[1])
        ]
        Vinf_launch = np.linalg.norm([Vinfx, Vinfy, Vinfz])

        # We now have the initial mass of the spacecraft
        m_initial = launchers.atlas551(Vinf_launch / 1000.0)

        # compute final flyby and resulting trajectory
        eph = self._seq[-1].eph(ep[-1])
        v_out = fb_prop(
            lamberts[-1].get_v2()[0],
            eph[1],
            x[-1] * self._seq[-1].radius,
            x[-2],
            self._seq[-1].mu_self,
        )
        a, e, i, W, w, E = ic2par(eph[0], v_out, self._common_mu)
        final_perihelion = a * (1 - e)
        # orbit should be as polar as possible, but we do not care about prograde/retrograde
        corrected_inclination = abs(abs(i) % pi - pi / 2)

        # check perihelion and aphelion bounds during the flight
        min_sun_distance = final_perihelion
        max_sun_distance = AU

        for l_i in range(self._n_legs):
            # project lambert leg, compute perihelion and aphelion
            eph = self._seq[l_i].eph(ep[l_i])
            ri = eph[0]
            if self._dummy_DSM[l_i]:
                fl_i = (
                    len(x) - 3 * sum(self._dummy_DSM) + (rf_index[l_i] - 1) * 3
                )
                tof_ratio, beta, r_p = x[
                    fl_i : fl_i + 3
                ]  # TODO: adapt this to different encodings
                if tof_ratio < 0 or tof_ratio > 1:
                    raise ValueError(
                        "Time of flight ratio of " + str(tof_ratio) + " is invalid."
                    )

                assert l_i > 0
                vi = lamberts[l_i - 1].get_v2()[0]
                # perform flyby, then propagate lagrangian
                vi = fb_prop(
                    vi, ri, r_p * self._seq[l_i].radius, beta, self._seq[l_i].mu_self
                )  # TODO: is seq[l_i] the correct planet?

                # then propagate after flyby
                if tof_ratio > 0:
                    ri, vi = propagate_lagrangian(
                        ri, vi, T[l_i] * DAY2SEC * tof_ratio, self._common_mu
                    )

            transfer_v = lamberts[l_i].get_v1()[0]
            transfer_a, transfer_e, _, _, _, E = ic2par(ri, transfer_v, self._common_mu)
            transfer_period = 2 * pi * sqrt(transfer_a ** 3 / self._common_mu)

            # check whether extremum happens during the transfer
            M = E - transfer_e * sin(E)
            mean_angle_to_apoapsis = pi - M
            if mean_angle_to_apoapsis < 0:
                mean_angle_to_apoapsis += 2 * pi
            mean_angle_to_periapsis = 2 * pi - M

            # update min and max sun distance
            if lamberts[l_i].get_tof() > mean_angle_to_apoapsis * transfer_period:
                max_sun_distance = max(max_sun_distance, transfer_a * (1 + transfer_e))

            if lamberts[l_i].get_tof() > mean_angle_to_periapsis * transfer_period:
                min_sun_distance = min(min_sun_distance, transfer_a * (1 - transfer_e))

        if self._multi_objective:
            return [
                corrected_inclination + 2 * min_sun_distance / AU,
                T,
                np.sum(DV) - 10,
                self._min_start_mass - m_initial,
                0.28 - min_sun_distance / AU,
                max_sun_distance / AU - 1.2,
            ]
        else:
            return [
                corrected_inclination + 2 * min_sun_distance / AU,
                np.sum(DV) - 10,
                self._min_start_mass - m_initial,
                0.28 - min_sun_distance / AU,
                max_sun_distance / AU - 1.2,
            ]

    def get_nobj(self):
        return self._multi_objective + 1

    def get_bounds(self):
        t0 = self._t0
        tof = self._tof
        n_legs = self._n_legs

        if self._tof_encoding == "alpha":
            # decision vector is  [t0, T, a1, a2, ....]
            lb = [t0[0].mjd2000, tof[0]] + [1e-3] * (n_legs)
            ub = [t0[1].mjd2000, tof[1]] + [1.0 - 1e-3] * (n_legs)
        elif self._tof_encoding == "direct":
            # decision vector is  [t0, T1, T2, T3, ... ]
            lb = [t0[0].mjd2000] + [it[0] for it in self._tof]
            ub = [t0[1].mjd2000] + [it[1] for it in self._tof]
        elif self._tof_encoding == "eta":
            # decision vector is  [t0, n1, n2, ....]
            lb = [t0[0].mjd2000] + [1e-3] * (n_legs)
            ub = [t0[1].mjd2000] + [1.0 - 1e-3] * (n_legs)

        # add final flyby
        pl = self._seq[-1]
        lb = lb + [-2 * pi, (pl.radius + self._safe_distance) / pl.radius]
        ub = ub + [2 * pi, 30]

        # something of a hack: parameters for dummy DSMs of possibly resonant flybys
        for i_fl in range(len(self._seq) - 1):
            if self._dummy_DSM[i_fl]:
                pl = self._seq[i_fl]
                lb = lb + [0, -2 * pi, (pl.radius + self._safe_distance) / pl.radius]
                ub = ub + [1, 2 * pi, 30]

        if self._tof_encoding == "alpha":
            assert len(lb) == 2 + n_legs + 2 + 3 * sum(self._dummy_DSM)
        elif self._tof_encoding in ["direct", "eta"]:
            assert len(lb) == 1 + n_legs + 2 + 3 * sum(self._dummy_DSM)
        elif self._tof_encoding == "eta":
            assert False

        assert len(ub) == len(lb)

        return (lb, ub)

    def get_nic(self):
        return 4

    def eph(self, x, t):
        if len(x) != len(self.get_bounds()[0]):
            raise ValueError(
                "Expected chromosome of length "
                + str(len(self.get_bounds()[0]))
                + " but got length "
                + str(len(x))
            )

        DVfb, lamberts, ep, _, _ = self._compute_dvs(x)
        rf_index = np.cumsum(self._dummy_DSM)
        if t <= ep[0]:
            raise ValueError(
                "Given epoch " + str(t) + " is at or before launch date " + str(ep[0])
            )

        i = 0
        while i < len(ep) and t >= ep[i]:
            # lambert leg i goes from planet i to planet i+1
            i += 1

        assert i >= 1 and i <= len(ep)
        if i < len(ep):
            assert t < ep[i]

        r_P, v_P = self._seq[i - 1].eph(ep[i - 1])

        elapsed_seconds = (t - ep[i - 1]) * DAY2SEC
        assert elapsed_seconds >= 0

        if i < len(ep):
            # get velocity from start of lambert leg i
            vel = lamberts[i - 1].get_v1()[0]
            if self._dummy_DSM[i - 1]:
                # adjust for DSM
                fl_i = (
                    len(x)
                    - 3 * sum(self._dummy_DSM)
                    + (rf_index[i - 1] - 1) * 3
                )
                tof_ratio, beta, r_p = x[fl_i : fl_i + 3]
                if tof_ratio < 0 or tof_ratio > 1:
                    raise ValueError(
                        "Time of flight ratio of " + str(tof_ratio) + " is invalid."
                    )
                seconds_until_DSM = tof_ratio * (ep[i] - ep[i - 1]) * DAY2SEC

                if elapsed_seconds < seconds_until_DSM:
                    # we are still in the part before the DSM, not in the lambert arc
                    assert i > 1  # resonances not yet supported at launch
                    vel = lamberts[i - 2].get_v2()[0]
                    # perform flyby, then propagate lagrangian
                    vel = fb_prop(
                        vel,
                        v_P,
                        r_p * self._seq[i - 1].radius,
                        beta,
                        self._seq[i - 1].mu_self,
                    )
        else:
            # get velocity after last flyby
            vel = fb_prop(
                lamberts[-1].get_v2()[0],
                v_P,
                x[-1] * self._seq[-1].radius,
                x[-2],
                self._seq[-1].mu_self,
            )

        # propagate the lagrangian
        r, v = propagate_lagrangian(r_P, vel, elapsed_seconds, self._common_mu)

        return r, v

    def pretty(self, x):
        """pretty(x)

        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).

        Prints human readable information on the trajectory represented by the decision vector x
        """
        T = self._decode_tofs(x)
        ep = np.insert(T, 0, x[0])  # [t0, T1, T2 ...]
        ep = np.cumsum(ep)  # [t0, t1, t2, ...]
        DVfb, l, _, _, _ = self._compute_dvs(x)  # TODO: reduce redundant computation of ep
        Vinfx, Vinfy, Vinfz = [
            a - b for a, b in zip(l[0].get_v1()[0], self._seq[0].eph(ep[0])[1])
        ]

        print("Multiple Gravity Assist (MGA) problem: ")
        print("Planet sequence: ", [pl.name for pl in self._seq])

        print("Departure: ", self._seq[0].name)
        print("\tEpoch: ", ep[0], " [mjd2000]")
        print("\tSpacecraft velocity: ", l[0].get_v1()[0], "[m/s]")
        print("\tLaunch velocity: ", [Vinfx, Vinfy, Vinfz], "[m/s]")
        _, _, transfer_i, _, _, _ = ic2par(
            self._seq[0].eph(ep[0])[0], l[0].get_v1()[0], self._common_mu
        )
        print("\tOutgoing Inclination:", transfer_i * RAD2DEG, "[deg]")

        for pl, e, dv, leg in zip(self._seq[1:-1], ep[1:-1], DVfb, l[1:]):
            print("Fly-by: ", pl.name)
            print("\tEpoch: ", e, " [mjd2000]")
            print("\tDV: ", dv, "[m/s]")
            eph = pl.eph(e)
            transfer_v = leg.get_v1()[0]
            _, _, transfer_i, _, _, _ = ic2par(eph[0], transfer_v, self._common_mu)
            print("\tOutgoing Inclination:", transfer_i * RAD2DEG, "[deg]")

        print("Final Fly-by: ", self._seq[-1].name)
        print("\tEpoch: ", ep[-1], " [mjd2000]")
        print("\tSpacecraft velocity: ", l[-1].get_v2()[0], "[m/s]")
        print("\tBeta: ", x[-2])
        print("\tr_p: ", x[-1])

        print("Resulting Solar orbit:")
        r_P, v_P = self._seq[-1].eph(ep[-1])
        v_out = fb_prop(
            l[-1].get_v2()[0],
            v_P,
            x[-1] * self._seq[-1].radius,
            x[-2],
            self._seq[-1].mu_self,
        )
        a, e, i, W, w, E = ic2par(r_P, v_out, self._common_mu)
        print("Perihelion: ", (a * (1 - e)) / AU, " AU")
        print("Aphelion: ", (a * (1 + e)) / AU, " AU")
        print("Inclination: ", i * RAD2DEG, " degrees")

        print("Time of flights: ", T, "[days]")

    def plot(self, x, axes=None, units=AU, N=60):
        """plot(self, x, axes=None, units=pk.AU, N=60)

        Plots the spacecraft trajectory.

        Args:
            - x (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - axes (``matplotlib.axes._subplots.Axes3DSubplot``): 3D axes to use for the plot
            - units (``float``, ``int``): Length unit by which to normalise data.
            - N (``float``): Number of points to plot per leg
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        from pykep.orbit_plots import plot_planet, plot_lambert, plot_kepler

        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams["legend.fontsize"] = 10
            fig = plt.figure()
            axes = fig.gca(projection="3d")

        T = self._decode_tofs(x)
        ep = np.insert(T, 0, x[0])  # [t0, T1, T2 ...]
        ep = np.cumsum(ep)  # [t0, t1, t2, ...]
        _, l, _, _, _ = self._compute_dvs(x)
        for pl, e in zip(self._seq, ep):
            plot_planet(
                pl, epoch(e), units=units, legend=True, color=(0.7, 0.7, 1), axes=axes
            )
        for lamb in l:
            plot_lambert(
                lamb,
                N=N,
                sol=0,
                units=units,
                color="k",
                legend=False,
                axes=axes,
                alpha=0.8,
            )

        # compute final flyby
        r_P, v_P = self._seq[-1].eph(ep[-1])
        v_out = fb_prop(
            l[-1].get_v2()[0],
            v_P,
            x[-1] * self._seq[-1].radius,
            x[-2],
            self._seq[-1].mu_self,
        )

        a, e, i, W, w, E = ic2par(r_P, v_out, self._common_mu)

        # final trajectory
        plot_kepler(
            r_P,
            v_out,
            365 * DAY2SEC,
            self._common_mu,
            N=100,
            color="r",
            units=units,
            axes=axes,
            label="Final Orbit",
        )

        return axes


solar_orbiter = _solar_orbiter_udp()
