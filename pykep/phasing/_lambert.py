from pygmo.problem._base import base # pylint: disable=import-error
from pygmo import hypervolume
from pykep.planet import gtoc7
from pykep.orbit_plots import plot_planet, plot_lambert
from pykep.core import lambert_problem, DAY2SEC, epoch, AU, damon
from numpy.linalg import norm
from math import exp


class lambert_metric(base):

    """
    This class defines the Lambert phasing metric as published in the paper:

    Hennes, Izzo, Landau: "Fast approximators for optimal low-thrust hops between main belt asteroids" - IEEE SSCI 2016

    The result is a pygmo multi-objective problem that can be solved efficiently by MO optimization algorithms
    """

    def __init__(self, epoch_bounds=[0, 1000], A1=gtoc7(1), A2=gtoc7(2), single_objective=False, Tmax=0.3, Isp=3000, ms=1500):
        """
        pykep.phasing.lambert_metric(epoch_bounds,A1, A2, max_acc, multi_objective)

        - epoch_bounds: a list containing the lower and upper bounds in mjd2000 for the launch and arrival epochs
        - A1: a planet
        - A2: a planet
        - Tmax: maximum spacecraft thrust [N]
        - Isp: specific impulse of the spacecarft propulsion system [s]
        - ms: spacecraft mass at dparture [kg]
        - single_objective: if True defines a single objectiove problem (only DV)

        Example::

        lm = planet(epoch_bounds=[0, 1000], A1=gtoc7(1), A2=gtoc7(2), single_objective=False, Tmax = 0.3, Isp = 3000, ms = 1500)
        """

        # First we call the constructor of the base class telling
        # essentially to pygmo what kind of problem to expect (2 objective, 0
        # contraints etc.)
        super().__init__(
            2, 0, 1 + (not single_objective), 0, 0, 0)

        # then we set the problem bounds (in this case equal for all
        # components)
        self.set_bounds(epoch_bounds[0], epoch_bounds[1])
        self._ast1 = A1
        self._ast2 = A2
        self._Tmax = Tmax
        self._Isp = Isp
        self._ms = ms
        self._UNFEASIBLE = 1e+20

    # We reimplement the virtual method that defines the objective function.
    def _objfun_impl(self, x):
        from math import sqrt
        # 1 - We check that the transfer time is positive
        if x[0] >= x[1]:
            if self.f_dimension == 1:
                return (self._UNFEASIBLE, )
            else:
                return (self._UNFEASIBLE, self._UNFEASIBLE)

        # 2 - We compute the asteroid positions
        r1, v1 = self._ast1.eph(x[0])
        r2, v2 = self._ast2.eph(x[1])

        # 3 - We compute Lambert arc
        l = lambert_problem(
            r1, r2, (x[1] - x[0]) * DAY2SEC, self._ast1.mu_central_body, False, 0)

        # 4 - We compute the two impulses
        v1l = l.get_v1()[0]
        v2l = l.get_v2()[0]
        DV1 = [a - b for a, b in zip(v1, v1l)]
        DV2 = [a - b for a, b in zip(v2, v2l)]

        a1, a2, tau, dv = damon(DV1, DV2, (x[1] - x[0]) * DAY2SEC)

        Isp = self._Isp
        g0 = 9.80665
        Tmax = self._Tmax
        ms = self._ms
        MIMA = 2 * Tmax / norm(a1) / (1. + exp(-norm(a1)
                                               * (x[1] - x[0]) * DAY2SEC / Isp / g0))

        DV1 = sum([l * l for l in DV1])
        DV2 = sum([l * l for l in DV2])
        totDV = sqrt(DV1) + sqrt(DV2)
        totDT = x[1] - x[0]

        if ms > MIMA:
            if self.f_dimension == 1:
                return (self._UNFEASIBLE, )
            else:
                return (self._UNFEASIBLE, self._UNFEASIBLE)

        if self.f_dimension == 1:
            return (totDV, )
        else:
            return (totDV, x[1])

    def plot_orbits(self, pop, ax=None):
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D

        A1, A2 = self._ast1, self._ast2

        if ax is None:
            fig = plt.figure()
            axis = fig.add_subplot(111, projection='3d')
        else:
            axis = ax

        plot_planet(A1, axes=axis, s=10, t0=epoch(self.lb[0]))
        plot_planet(A2, axes=axis, s=10, t0=epoch(self.ub[0]))
        for ind in pop:
            if ind.cur_f[0] == self._UNFEASIBLE:
                continue
            dep, arr = ind.cur_x
            rdep, vdep = A1.eph(epoch(dep))
            rarr, varr = A2.eph(epoch(arr))
            l = lambert_problem(rdep, rarr, (arr - dep) *
                                DAY2SEC, A1.mu_central_body, False, 1)
            axis = plot_lambert(l, axes=axis, alpha=0.8, color='k')

        if ax is None:
            plt.show()

        return axis

    def plot_pareto_front(self, pop):
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D

        if pop.champion.f[0] == self._UNFEASIBLE:
            raise Exception(
                'Input population contains only unfeasible individuals')

        rx, ry = self._compute_ref_point()
        axis = pop.plot_pareto_fronts()
        axis.set_xlim([0, rx])
        axis.set_ylim([self.lb[0], ry])
        plt.xlabel("[m/s]")
        plt.ylabel("[MJD2000]")
        plt.draw()

        return axis

    def compute_hypervolume(self, pop):
        if pop.champion.f[0] == self._UNFEASIBLE:
            raise Exception(
                'Input population contains only unfeasible individuals')
        hv = hypervolume(pop)
        return (hv.compute(self._compute_ref_point()) * DAY2SEC / AU)

    def _compute_ref_point(self):
        rx = self._Tmax / self._ms * DAY2SEC * (self.ub[0] - self.lb[0])
        ry = self.ub[0]
        return (rx, ry)

    def human_readable_extra(self):
        retval = "\n\tAsteroid 1: " + self._ast1.name
        retval = retval + "\n\tAsteroid 2: " + self._ast2.name
        return retval
