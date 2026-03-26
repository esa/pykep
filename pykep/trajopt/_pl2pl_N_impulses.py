import pykep as pk
import numpy as _np
import matplotlib.pyplot as _plt

from math import pi, cos, sin, log, acos, sqrt


# Avoiding scipy dependency
def norm(x):
    return sqrt(sum([it * it for it in x]))


class pl2pl_N_impulses:
    """
    This class is a pygmo (http://esa.github.io/pygmo/) UDP representing a single leg transfer
    between two planets allowing up to a maximum number of impulsive Deep Space Maneuvers.

    The decision vector is::

      [t0,T] + [alpha,u,v,V_inf] * (N-2) + [alpha] + ([tf])

    ... in the units: [mjd2000, days] + [nd, nd, m/sec, nd] + [nd] + [mjd2000]

    Each time-of-flight can be decoded as follows, T_n = T log(alpha_n) / \\sum_i(log(alpha_i))

    .. note::

       The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.
    """

    def __init__(
        self,
        start=pk.planet(pk.udpla.jpl_lp("earth")),
        target=pk.planet(pk.udpla.jpl_lp("venus")),
        N_max=3,
        tof_bounds=[20.0, 400.0],
        DV_max_bounds=[0.0, 4.0],
        phase_free=True,
        multi_objective=False,
        t0_bounds=None,
    ):
        """pykep.trajopt.pl2pl_N_impulses(start='earth', target='venus', N_max=3, tof=[20., 400.], vinf=[0., 4.], phase_free=True, multi_objective=False, t0=None)

        Args:
            *start* (:class:`~pykep.planet`): initial body.

            *target* (:class:`~pykep.planet`): target body.

            *N_max* (:class:`int`): maximum number of impulses.

            *tof_bounds* (:class:`list`): the box bounds [lower,upper] for the total time of flight (days).

            *DV_max_bounds* (:class:`list`): the box bounds [lower,upper] for each DV magnitude (km/sec).

            *phase_free* (:class:`bool`): when True, no rendezvous condition are enforced and start and arrival anomalies will be free.

            *multi_objective* (:class:`bool`):  when True, a multi-objective problem is constructed with DV and time of flight as objectives.

            *t0_bounds* (:class:`list`):  the box bounds on the launch window containing two pykep.epoch. This is not needed if phase_free is True.
        """

        # Sanity checks
        # 0) This is not working for only two impulses
        if N_max <= 2:
            raise ValueError(
                "This UDP is not wsuitable for only two impulse trajectories. Lambert multiple revolutions should be allowed for that (to be implemented)."
            )
        # 1) all planets need to have the same mu_central_body
        if start.mu_central_body != target.mu_central_body:
            raise ValueError(
                "Starting and ending pykep.planet must have the same mu_central_body"
            )
        # 2) Number of impulses must be at least 2
        if N_max < 2:
            raise ValueError("Number of impulses N is less than 2")
        # 3) If phase_free is True, t0 does not make sense
        if t0_bounds is None and not phase_free:
            t0 = [pk.epoch(0), pk.epoch(1000)]
        if t0_bounds is not None and phase_free:
            raise ValueError("When phase_free is True no t0 can be specified")
        if not phase_free:
            if type(t0_bounds[0]) != type(pk.epoch(0)):
                t0_bounds[0] = pk.epoch(t0_bounds[0])
            if type(t0_bounds[1]) != type(pk.epoch(0)):
                t0_bounds[1] = pk.epoch(t0_bounds[1])

        self.obj_dim = multi_objective + 1
        # We then define all class data members
        self.start = start
        self.target = target
        self.N_max = N_max
        self.phase_free = phase_free
        self.multi_objective = multi_objective
        self.DV_max = [s * 1000 for s in DV_max_bounds]

        self._common_mu = start.mu_central_body

        # And we compute the bounds
        if phase_free:
            self._lb = (
                [0, tof_bounds[0]]
                + [0.0, 0.0, 0.0, DV_max_bounds[0] * 1000] * (N_max - 2)
                + [0]
                + [0]
            )
            self._ub = (
                [2 * start.period() * pk.SEC2DAY, tof_bounds[1]]
                + [1.0, 1.0, 1.0, DV_max_bounds[1] * 1000] * (N_max - 2)
                + [1.0]
                + [2 * target.period() * pk.SEC2DAY]
            )
        else:
            self._lb = (
                [t0_bounds[0].mjd2000, tof_bounds[0]]
                + [1e-3, 0.0, 0.0, DV_max_bounds[0] * 1000] * (N_max - 2)
                + [1e-3]
            )
            self._ub = (
                [t0_bounds[1].mjd2000, tof_bounds[1]]
                + [1.0 - 1e-3, 1.0, 1.0, DV_max_bounds[1] * 1000] * (N_max - 2)
                + [1.0 - 1e-3]
            )

    def get_nobj(self):
        return self.obj_dim

    def get_bounds(self):
        return (self._lb, self._ub)

    def decode(self, x):
        """
        Args:
            *x* (:class:`list`): The decision vector in the correct tof encoding.
        """
        retval = []
        # 1 We decode the tofs
        T = pk.alpha2direct(x[2::4], x[1])

        # 2 - We compute the starting and ending position
        r_start, v_start = self.start.eph(pk.epoch(x[0]))
        if self.phase_free:
            r_target, v_target = self.target.eph(pk.epoch(x[-1]))
        else:
            r_target, v_target = self.target.eph(pk.epoch(x[0] + x[1]))

        # 3 - We loop across inner impulses
        rsc = r_start
        vsc = v_start
        for i, time in enumerate(T[:-1]):
            DV = pk.utils.uvV2cartesian(x[3 + 4 * i : 6 + 4 * i])
            retval.append([[rsc, vsc], DV, T[i] * pk.DAY2SEC])

            # We apply the (i+1)-th impulse
            vsc = [a + b for a, b in zip(vsc, DV)]
            rsc, vsc = pk.propagate_lagrangian(
                [rsc, vsc], T[i] * pk.DAY2SEC, self._common_mu
            )
        cw = pk.ic2par([rsc, vsc], self.start.mu_central_body)[2] > pi / 2

        # We now compute the remaining two final impulses
        # Lambert arc to reach seq[1]
        dt = T[-1] * pk.DAY2SEC
        l = pk.lambert_problem(rsc, r_target, dt, self._common_mu, cw, False)
        v_beg_l = l.v0[0]
        v_end_l = l.v1[0]

        DV1 = [a - b for a, b in zip(v_beg_l, vsc)]
        retval.append([[rsc, vsc], DV1, T[-1] * pk.DAY2SEC])

        DV2 = [a - b for a, b in zip(v_target, v_end_l)]
        retval.append([[r_target, v_end_l], DV2, 0])
        return retval

    def fitness(self, x):
        mit = self.decode(x)

        DVs = [norm(node[1]) for node in mit]
        if self.obj_dim == 1:
            return (sum(DVs),)
        else:
            return (sum(DVs), x[1])

    def plot(
        self,
        x,
        ax=None,
        units=pk.AU,
        N=60,
        c_orbit="dimgray",
        c_segments=["royalblue", "indianred"],
        figsize=(5, 5),
        **kwargs
    ):
        """
        Plots the trajectory encoded into *x* in 3D axes.

        Args:
            *x* (:class:`list`): The decision vector in the correct tof encoding.

            *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`, optional): The 3D axis to plot on. Defaults to None.

            *units* (:class:`float`, optional): The unit scale for the plot. Defaults to pk.AU.

            *N* (:class:`int`, optional): The number of points to use when plotting the trajectory. Defaults to 60.

            *c_orbit* (:class:`str`, optional): The color of the planet orbits. Defaults to 'dimgray'.

            *c_segments* (:class:`list`, optional): The colors to alternate the various trajectory segments (inbetween DSMs). Defaults to ["royalblue", "indianred"].

            *figsize* (:class:`tuple`): The figure size (only used if a*ax* is None and axis have to be created.), Defaults to (5, 5).

            *\\*\\*kwargs*: Additional keyword arguments to pass to the trajectory plot (common to Lambert arcs and ballistic arcs)

        Returns:
            :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The 3D axis where the trajectory was plotted.
        """
        if ax is None:
            ax = pk.plot.make_3Daxis(figsize=figsize)

        # Adding the main central body (Sun-like)
        ax = pk.plot.add_sun(ax=ax)

        ax = pk.plot.add_planet_orbit(
            pla=self.start, ax=ax, units=units, N=N, c=c_orbit, label=self.start.name
        )
        ax = pk.plot.add_planet_orbit(
            pla=self.target, ax=ax, units=units, N=N, c=c_orbit, label=self.target.name
        )

        # We decode the chromosome
        mit = self.decode(x)  # [[r,v], DV, DT]

        ax = pk.plot.add_mit(
            ax, mit, self._common_mu, units=units, c_segments=c_segments, N=N, **kwargs
        )
        return ax

    def pretty(self, x):
        # We decode the chromosome
        mit = self.decode(x)  # [[r,v], DV, DT]

        DVs = [norm(node[1]) for node in mit]
        T = [node[2] for node in mit]

        print("Total DV (m/s): ", sum(DVs))
        print("Dvs (m/s): ", [float(it) for it in DVs])
        print("Total DT (m/s): ", sum(T))
        print(
            "Tofs (days): ", [float(it) for it in T[:-1]]
        )  # last node has a zero TOF by convention

    def plot_primer_vector(self, x, N=200, ax=None):
        """Plots the primer vector magnitude along the trajectory encoded in *x*.

        Args:
            *x* (:class:`list`): The decision vector in the correct tof encoding.

            *N* (:class:`int`, optional): The number of points to use when plotting the primer vector. Defaults to 200.

            *ax* (:class:`matplotlib.axes.Axes`, optional): The axis to plot on. Defaults to None.

        Returns:
            :class:`matplotlib.axes.Axes`: The axis where the primer vector was plotted.
            :class:`tuple`: A tuple containing the grid and the primer vector magnitude.
        """
        # We start by decoding the chromosome into the structure [[r,v], DV, DT]
        decoded = self.decode(x)

        # We explicitly extract the encoded information
        dts = [it[2] for it in decoded]
        DVs = [it[1] for it in decoded]
        norm_DVs = [norm(it) for it in DVs]

        posvels = [it[0] for it in decoded]

        if min(norm_DVs) < 1e-3:
            raise ValueError(
                "Impulse magnitude too small, primer vector computation is not possible. Decrease the number of impulses."
            )

        # We create one grid per segment (e.g. part of the trajectory between two impulses)
        # (this is not guaranteed to have the requested size N, nor has uniform spacing, since all impulses
        # must belong to the grid points)
        N = N + len(
            DVs
        )  # heuristic to make sure we are close to the requested number of points
        tgrids = [
            _np.linspace(
                sum(dts[:i]),
                sum(dts[: i + 1]),
                max(
                    int(dts[i] // (sum(dts) / (N - 1))), 5
                ),  # we force a minimum 5 points per segment
            )
            for i in range(len(dts) - 1)
        ]
        # We assemble all the grids into one single final_grid
        final_grid = _np.array([0])
        for i in range(len(dts) - 1):
            final_grid = _np.concatenate((final_grid, tgrids[i][1:]))

        # These are the indices of the final_grid where the impulses are given.
        idxs = [0] + [len(tgrids[0]) - 1]
        for grid in tgrids[1:]:
            idxs += [idxs[-1] + len(grid) - 1]

        # We now compute the various STMs.
        retvals = []
        for posvel, DV, tgrid in zip(posvels, DVs, tgrids):
            ic = [posvel[0], [a + b for a, b in zip(posvel[1], DV)]]
            retvals.append(
                pk.propagate_lagrangian_grid(ic, tgrid, mu=pk.MU_SUN, stm=True)
            )

        # And now assemble them in correspondance to the final_grid and in the Mn0 form.
        posvels = [item[0] for item in retvals[0]]
        stms = [item[1] for item in retvals[0]]

        M = stms[-1]
        for retval in retvals[1:]:
            posvels = posvels + [item[0] for item in retval[1:]]
            stms = stms + [item[1] @ M for item in retval[1:]]
            M = stms[-1]

        res = []
        # When computing the primer vector we must choose which impulses to use.
        # We choose the first and last impulse. But we could choose any pair of impulses,
        # and if the trajectory is optimal (locally) the primer vector would not change.
        idx_i = idxs[0]
        idx_j = idxs[-1]
        DVi = DVs[0]
        DVj = DVs[-1]
        for idx_k, _ in enumerate(final_grid):
            Mji = stms[idx_j] @ _np.linalg.inv(stms[idx_i])
            Mjk = stms[idx_j] @ _np.linalg.inv(stms[idx_k])
            res.append(
                _np.linalg.norm(pk.trajopt.primer_vector(DVi, DVj, Mji, Mjk)[0])
            )

        if ax is None:
            ax = _plt.figure().add_subplot()
        ax.plot(res, label="primer vector magnitude")
        ax.vlines(0, 0.8, 1.2, "k", linestyles="dashed", label="impulse")
        for idx in idxs:
            ax.vlines(idx, 0.8, 1.2, "k", linestyles="dashed")

        ax.hlines(1, 0, len(final_grid), "r")
        ax.legend(loc="lower right")
        return ax, (final_grid, res)
