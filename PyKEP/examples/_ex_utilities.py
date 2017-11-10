# We define a meta-problem augmenting a generic UDP with a numercial gradient so that algorithms
# needing the gradient can be used. Note that algorithm such as snopt7 have an internal system to compute gradients
# in case the problem is fed to it without. Hence the kwarg with_grad is provided to deactivate the gradient.

import pygmo as pg
class add_gradient:
    def __init__(self, udp, with_grad=False):
        self.udp_inner = udp
        self.prob = pg.problem(udp)
        self.with_grad = with_grad

    def fitness(self, x):
        return self.prob.fitness(x)

    def get_bounds(self):
        return self.prob.get_bounds()

    def get_nec(self):
        return self.prob.get_nec()

    def get_nic(self):
        return self.prob.get_nic()

    def get_nobj(self):
        return self.prob.get_nobj()

    def gradient(self, x):
        # we here use the low precision gradient
        return pg.estimate_gradient(lambda x: self.fitness(x), x, 1e-8)

    def has_gradient(self):
        return self.with_grad

# We define a utility class that sets the algorithm to be used (in pagmo's terminology the UDA or user defined algorithm).
def algo_factory(name):
    if name is "slsqp":
        uda = pg.nlopt('slsqp')
        uda.xtol_rel = 1e-6
        uda.ftol_rel = 1e-10
        algo = pg.algorithm(uda)
        algo.set_verbosity(1)
        return algo
    elif name is "ipopt":
        uda = pg.ipopt()
        uda.set_integer_option("print_level", 5)
        uda.set_integer_option("acceptable_iter", 4)
        uda.set_integer_option("max_iter", 150)

        uda.set_numeric_option("tol", 1e-8)
        uda.set_numeric_option("dual_inf_tol", 1e-8)
        uda.set_numeric_option("constr_viol_tol", 1e-8)
        uda.set_numeric_option("compl_inf_tol", 1e-8)

        uda.set_numeric_option("acceptable_tol", 1e-3)
        uda.set_numeric_option("acceptable_dual_inf_tol", 1e-2)
        uda.set_numeric_option("acceptable_constr_viol_tol", 1e-6)
        uda.set_numeric_option("acceptable_compl_inf_tol", 1e-6)

        algo = pg.algorithm(uda)
        return algo
    elif name is "snopt7":
        import pygmo_plugins_nonfree as pg7
        uda = pg7.snopt7(True, "/usr/local/lib/libsnopt7_c.so")
        uda.set_integer_option("Major iterations limit", 2000)
        uda.set_integer_option("Iterations limit", 200000)
        uda.set_numeric_option("Major optimality tolerance", 1e-2)
        uda.set_numeric_option("Major feasibility tolerance", 1e-9)

        algo = pg.algorithm(uda)
        return algo