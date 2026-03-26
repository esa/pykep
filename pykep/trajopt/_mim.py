import pykep as pk
from copy import deepcopy as _deepcopy
import heyoka as _hy
import pygmo as _pg


class _mim:
    def __init__(
        self,
        source,
        target,
        t0=pk.epoch(0),
        tof_days=250.0,
        T_max=0.6,
        veff=4000.0 * pk.G0,
        L=pk.AU,
        TIME=pk.YEAR2DAY * pk.DAY2SEC,
        MASS=1500,
        with_gradient=False,
    ):
        # Non-dimensional units
        VEL = L / TIME  # Unit for velocity (1 AU/year)
        ACC = VEL / TIME  # Unit for acceleration (1 AU/year^2)

        # Store user inputs
        self.veff = veff
        self.t0 = t0

        # We redefine the user inputs in non dimensional units.
        self.mu = source.get_mu_central_body() / (L**3 / TIME**2)
        self.c1 = T_max / (MASS * ACC)
        self.c2 = (veff) / VEL
        self.tof = tof_days * pk.DAY2SEC / TIME

        # Initial position is computed once only upon construction.
        r0, v0 = source.eph(t0)
        r0 = [it / L for it in r0]
        v0 = [it / VEL for it in v0]
        self.posvel0 = [r0, v0]

        # Final position is computed once only upon construction.
        rf, vf = target.eph(t0 + tof_days)
        rf = [it / L for it in rf]
        vf = [it / VEL for it in vf]
        self.posvelf = [rf, vf]

        # Storing the planets
        self.source = source
        self.target = target

        # And the Taylor integrators
        self.ta = pk.ta.get_pc(1e-16, pk.optimality_type.TIME)
        dyn = pk.ta.pc_dyn(pk.optimality_type.TIME)
        args = _hy.make_vars("lx", "ly", "lz", "lvx", "lvy", "lvz", "lm", "m")
        aug_dyn = _hy.var_ode_sys(sys=dyn, args=args, order=1)
        self.ta_var = _hy.taylor_adaptive(aug_dyn, tol=1e-8, compact_mode=True)
        self.ic_var = _deepcopy(self.ta_var.state[14:])
        # Compiled functions
        self.dyn_func = pk.ta.get_pc_dyn_cfunc(pk.optimality_type.TIME)

        # Non dimensional units
        self.MASS = MASS
        self.L = L
        self.TIME = TIME
        self.VEL = VEL
        self.ACC = ACC

        # Do we provide the gradient
        self.with_gradient = with_gradient

    def get_bounds(self):
        lb = [-1.0] * 6 + [0.0] * 2
        ub = [1.0] * 7 + [20.0]
        return [lb, ub]

    def set_ta_state(self, x):
        # Preparing the numerical integration parameters
        self.ta.pars[0] = self.mu
        self.ta.pars[1] = self.c1
        self.ta.pars[2] = self.c2
        self.ta.time = 0.0

        # And initial conditions
        self.ta.state[:3] = self.posvel0[0]
        self.ta.state[3:6] = self.posvel0[1]
        self.ta.state[6] = x[7]
        self.ta.state[7:14] = x[:7]

    def set_ta_var_state(self, x):
        # Preparing the numerical integration parameters
        self.ta_var.pars[0] = self.mu
        self.ta_var.pars[1] = self.c1
        self.ta_var.pars[2] = self.c2
        self.ta_var.time = 0.0

        # And initial conditions
        self.ta_var.state[:3] = self.posvel0[0]
        self.ta_var.state[3:6] = self.posvel0[1]
        self.ta_var.state[6] = x[7]
        self.ta_var.state[7:14] = x[:7]
        self.ta_var.state[14:] = self.ic_var

    def fitness(self, x):
        # x = [lx, ly, lz, lvx, lvy, lvz, lm, m]
        # Single Shooting
        self.set_ta_state(x)
        self.ta.propagate_until(self.tof)

        # Assembling the constraints (no need for lm ....)
        ceq = [self.ta.state[0] - self.posvelf[0][0]]
        ceq += [self.ta.state[1] - self.posvelf[0][1]]
        ceq += [self.ta.state[2] - self.posvelf[0][2]]
        ceq += [self.ta.state[3] - self.posvelf[1][0]]
        ceq += [self.ta.state[4] - self.posvelf[1][1]]
        ceq += [self.ta.state[5] - self.posvelf[1][2]]
        ceq += [self.ta.state[13]]  # Not necessary but squares the problem lm = 0
        ceq += [sum([it * it for it in x[:7]]) - 1.0]  # |lambdas|^2 = 1
        return [1.0] + ceq

    def gradient(self, x):
        # Single Shooting of variational equations
        self.set_ta_var_state(x)
        self.ta_var.propagate_until(self.tof)

        # Assembling the gradient
        # Fitness (sparsity takes care of the fitness since it does not depend on the decision vector)
        retval = []
        # State Constraints
        for i in range(6):
            sl = self.ta_var.get_vslice(order=1, component=i)
            retval.extend(list(self.ta_var.state[sl]))
        # Costate Constraint
        sl = self.ta_var.get_vslice(order=1, component=13)
        retval.extend(list(self.ta_var.state[sl]))
        # Norm constraint
        for item in x[:7]:
            retval.extend([2 * item])
        return retval

    def gradient_sparsity(self):
        # x = [lx,ly,lz,lvx,lvy,lvz,lm,m]
        retval = []
        # The final constraints
        for i in range(1, 8):
            for j in range(8):
                retval.append((i, j))
        # The norm constraint does not depend on m
        for i in range(8, 9):
            for j in range(7):
                retval.append((i, j))
        return retval

    def has_gradient(self):
        return self.with_gradient

    def get_nec(self):
        return 8

def mim_from_hop(pl_s, pl_t, when_s, when_f, Tmax, veff, MASS=1000.0, max_trials=10):
    """mim_from_hop(pl_s, pl_t, when_s, when_t, Tmax, veff, MASS=1000.0, max_trials=10)

    The Maximum Initial Mass from a hop transfer. A hop transfer is fixed time randevouz between two planets.
    This function is slower then its approximations :class:`~pykep.mima_from_hop` and :class:`~pykep.mima2_from_hop`,
    if succesfull, returns the exact solution to the problem. Behind the scene Potryagin Maximum Principle
    is applied to the problem and the TPBVP is solved using ipopt / heyoka.

    Args:

        *pl_s* (:class:`~pykep.planet`): the source planet

        *pl_t* (:class:`~pykep.planet`): the target planet

        *when_s* (:class:`~pykep.epoch`): the start epoch

        *when_f* (:class:`~pykep.epoch`): the final epoch

        *Tmax* (:class:`float`, optional): Maximum spacecraft thrust.

        *veff* (:class:`float`, optional): Isp*G0.

        *MASS* (:class:`float`, optional): Units to use to scale the mass intrnally. It also influences the bounds for the MIM in the
        optimization problem.           

        *max_trials* (:class:`int`, optional): Maximum number of trials to find the MIM. If the algorithm does not converge, it will
        return None.                        

    Returns:
        :class:`float`: The Maximum Initial Mass

    Examples:   

        >>> import pykep as pk
        >>> ... # assuming to have two planets pl_s and pl_t
        >>> when_s = pk.epoch(0.0)
        >>> when_f = pk.epoch(100.0)
        >>> Tmax = 0.6
        >>> veff = 3000*pk.G0
        >>> mima = pk.mim_from_hop(pl_s, pl_t, when_s, when_f, Tmax, veff)
        >>> print("Maximum initial mass:", mima, "kg")  
    """
    # We setup ipopt as a solver
    ip = _pg.ipopt()
    ip.set_numeric_option("tol", 1e-5)  # Change the relative convergence tolerance
    ip.set_integer_option("max_iter", 50)  # Change the maximum iterations
    ip.set_integer_option("print_level", 0)  # Makes Ipopt unverbose
    ip.set_string_option(
        "nlp_scaling_method", "none"
    )  # Removes any scaling made in auto mode
    ip.set_string_option(
        "mu_strategy", "adaptive"
    )  # Alternative is to tune the initial mu value
    ipopt = _pg.algorithm(ip)

    # We set up the mim problem
    tof_in_days = when_f.mjd2000 - when_s.mjd2000
    L = pl_s.elements()[0]
    TIME = pl_s.period()
    udp = _mim(
        pl_s,
        pl_t,
        when_s,
        tof_in_days,
        Tmax,
        veff,
        L,
        TIME,
        MASS,
        with_gradient=True,
    )
    prob = _pg.problem(udp)
    prob.c_tol = 1e-6

    # And solve
    masses = []
    xs = []
    for i in range(max_trials):
        pop = _pg.population(prob, 1)
        pop = ipopt.evolve(pop)
        if prob.feasibility_f(pop.champion_f):
            return pop.champion_x[7] * udp.MASS
    print(
        "Could not converge to find the MIM. tof might be too long, mass bounds reached, or the algo may need more iterations (to be manually set)"
    )
