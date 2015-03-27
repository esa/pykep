from PyGMO.problem._base import base
from PyKEP.planets import gtoc7
from PyKEP.core import lambert_problem, DAY2SEC


class lambert_metric(base):

    """
    This class defines the Lambert phasing metric as introduced by the ESA/ACT team during GTOC7.
    The result is a PyGMO problem that can be solved efficiently by optimization algorithms
    """

    def __init__(self, epoch_bounds=[0, 1000], A1=gtoc7(1), A2=gtoc7(2), multi_objective=False, max_acc=5e-4):
        """
PyKEP.phasing.lambert_metric(epoch_bounds,A1, A2)

- epoch_bounds: a list containing the lower and upper bounds in mjd2000 for the launch and arrival epochs
- A1: a planet
- A2: a planet
- max_acc: maximum acceleration from the thrust [m/s^2]
- multi_objective: if True creates a multi-objective problem

Example::

  lm = planet(epoch_bounds=[0, 1000], A1=gtoc7(1), A2=gtoc7(2))
        """

        # First we call the constructor of the base class telling
        # essentially to PyGMO what kind of problem to expect (1 objective, 0
        # contraints etc.)
        super(lambert_metric, self).__init__(2, 0, 1 + multi_objective, 0, 0, 0)

        # then we set the problem bounds (in this case equal for all
        # components)
        self.set_bounds(epoch_bounds[0], epoch_bounds[1])
        self._ast1 = A1
        self._ast2 = A2
        self._max_acc = max_acc

    # We reimplement the virtual method that defines the objective function.
    def _objfun_impl(self, x):
        from math import sqrt
        # 1 - We check that the transfer time is positive
        if x[0] >= x[1]:
            if self.f_dimension == 1:
                return (1e20, )
            else:
                return (1e20, 1e20)

        # 2 - We compute the asteroid positions
        r1, v1 = self._ast1.eph(x[0])
        r2, v2 = self._ast2.eph(x[1])

        # 3 - We compute Lambert arc
        l = lambert_problem(r1, r2, (x[1] - x[0]) * DAY2SEC, self._ast1.mu_central_body, False, 0)

        # 4 - We compute the two impulses
        v1l = l.get_v1()[0]
        v2l = l.get_v2()[0]
        DV1 = [a - b for a, b in zip(v1, v1l)]
        DV2 = [a - b for a, b in zip(v2, v2l)]
        DV1 = sum([l * l for l in DV1])
        DV2 = sum([l * l for l in DV2])
        totDV = sqrt(DV1) + sqrt(DV2)
        totDT = x[1] - x[0]

        if totDV >= self._max_acc * totDT * DAY2SEC:
            if self.f_dimension == 1:
                return (1e20, )
            else:
                return (1e20, 1e20)

        if self.f_dimension == 1:
            return (totDV, )
        else:
            return (totDV, totDT)

    def human_readable_extra(self):
        retval = "\n\tAsetroid 1: " + self._ast1.name
        retval = retval + "\n\tAsetroid 2: " + self._ast2.name
        return retval
