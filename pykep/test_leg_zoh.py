import pykep as pk
import numpy as np
import pygmo as pg

from copy import deepcopy
import unittest as _ut

# Helper to compute numerical gradients of the mc.
def compute_mismatch_constraints_n(leg, state0, controls, state1, tgrid):
   leg.tgrid = tgrid
   leg.state0 = state0
   leg.state1 = state1
   leg.controls = controls
   leg.state1 = state1
   return leg.compute_mismatch_constraints()

def compute_throttle_constraints_n(leg, controls):
    leg.controls = controls
    return leg.compute_throttle_constraints()

# The actual tests.
class leg_zoh_test(_ut.TestCase):
    def test_zoh_mc_ballistic(self):
        import numpy as np
        import pykep as pk

        # The integrators (Keplerian Propagation in Cartesian)
        tol = 1e-12
        ta = pk.ta.get_zoh_kep(tol)

        # We create a ballistic leg
        t0 = 10000  # MJD2000
        t1 = 10400  # MJD2000
        pl0 = pk.planet(pk.udpla.jpl_lp("Earth"))
        pl1 = pk.planet(pk.udpla.jpl_lp("Mars"))
        r0, _ = pl0.eph(t0)
        r1, _ = pl1.eph(t1)
        l = pk.lambert_problem(r0=r0, r1=r1, tof=(t1 - t0) * pk.DAY2SEC, mu=pk.MU_SUN)
        m0 = 1000
        m1 = 1000

        # nd units (in these units mu must be one as to use the tas)
        L = pk.AU
        MU = pk.MU_SUN
        TIME = np.sqrt(L**3 / MU)
        V = L / TIME
        ACC = V / TIME
        MASS = 1000

        # We instantiate a ballistic leg and test for small mismatches
        n_trials = 50
        for i in range(n_trials):
            # leg data
            nseg = int(np.random.uniform(4, 20))
            veff = np.random.uniform(2000, 7000) * pk.G0

            # nd construction data
            state0 = [it / L for it in r0] + [it / V for it in l.v0[0]] + [m0 / MASS]
            state1 = [it / L for it in r1] + [it / V for it in l.v1[0]] + [m1 / MASS]
            veff_nd = veff / V
            tgrid = np.linspace(
                t0 * pk.DAY2SEC / TIME, t1 * pk.DAY2SEC / TIME, nseg + 1
            )
            controls = np.random.uniform(-1, 1, (4 * nseg,)) * 1e-3
            controls[0::4] = 0  # zeroing the thrust magnitude
            cut = np.random.uniform(0, 1)
            ta.pars[4] = 1.0 / veff_nd
            # construct the leg
            leg = pk.leg.zoh(
                state0, controls.tolist(), state1, tgrid, cut=cut, tas=[ta, None]
            )
            # test
            mc = leg.compute_mismatch_constraints()
            self.assertTrue(np.array([np.max(i) < 1e-12 for i in mc]).all())

    def test_zoh_mc_thrust(self):
        # The integrator (Keplerian Propagation in Cartesian)
        tol = 1e-12
        ta = pk.ta.get_zoh_kep(tol)

        # We create some boundary conditions
        t0 = 1234  # MJD2000
        t1 = 3456  # MJD2000
        pl0 = pk.planet(pk.udpla.jpl_lp("Venus"))
        pl1 = pk.planet(pk.udpla.jpl_lp("Earth"))
        r0, _ = pl0.eph(t0)
        r1, _ = pl1.eph(t1)
        l = pk.lambert_problem(r0=r0, r1=r1, tof=(t1 - t0) * pk.DAY2SEC, mu=pk.MU_SUN)
        v0 = l.v0[0]
        v1 = l.v1[0]
        m0 = 1000
        m1 = 1000

        # nd units (in these units mu must be one as to use the tas)
        L = pk.AU
        MU = pk.MU_SUN
        TIME = np.sqrt(L**3 / MU)
        V = L / TIME
        ACC = V / TIME
        MASS = 1000
        F = MASS * ACC

        # leg data
        nseg = 5
        veff = 6000 * pk.G0  # Isp G0

        # nd construction data
        state0 = [it / L for it in r0] + [it / V for it in v0] + [m0 / MASS]
        state1 = [it / L for it in r1] + [it / V for it in v1] + [m1 / MASS]
        veff_nd = veff / V
        tgrid_nd = np.linspace(t0 * pk.DAY2SEC / TIME, t1 * pk.DAY2SEC / TIME, nseg + 1)
        controls_nd = np.array(
            [0.03 / F, 1.0, 0.0, 0]
            + [0.02 / F, 1.0, 0.0, 0]
            + [0.003 / F, 0.0, 1.0, 0]
            + [0.03 / F, 1.0, 0.0, 0]
            + [0.001 / F, 0.0, 0.0, 1.0]
        )
        cut = 0.5
        ta.pars[4] = 1.0 / veff_nd
        # construct the leg
        leg = pk.leg.zoh(
            state0, controls_nd.tolist(), state1, tgrid_nd, cut=cut, tas=[ta, None]
        )
        # test
        mc = leg.compute_mismatch_constraints()
        mc_gt = [
            0.021765250746026865,
            -0.10717348832415041,
            -0.014802731334149222,
            0.14380114520179108,
            -0.04229446021097845,
            -0.004488343342377676,
            -0.054814461615332655,
        ]
        diff = [abs(a - b) for a, b in zip(mc, mc_gt)]
        self.assertTrue(np.array([np.max(i) < 1e-12 for i in diff]).all())
        
    def test_compute_mc_grad(self):
        # The integrators (Keplerian Propagation in Cartesian)
        tol=1e-14
        tol_var = 1e-10

        ta = pk.ta.get_zoh_kep(tol)
        ta_var = pk.ta.get_zoh_kep_var(tol_var)
        
        # We create a ballistic leg 
        t0 = 10000
        t1 = 10400
        pl0 = pk.planet(pk.udpla.jpl_lp("Earth"))
        pl1 = pk.planet(pk.udpla.jpl_lp("Mars"))
        r0, _ = pl0.eph(t0)
        r1, _ = pl1.eph(t1)
        # We create some starting conditions from a Lambert arc
        l = pk.lambert_problem(r0=r0, r1=r1, tof = (t1-t0) * pk.DAY2SEC, mu = pk.MU_SUN)
        m0 = 1000
        m1 = 1000

        # nd units
        L = pk.AU
        MU = pk.MU_SUN
        TIME = np.sqrt(L**3/MU)
        V =  L/TIME
        ACC = V/TIME
        MASS = 1000
        F = MASS*ACC
        
        # nd data
        state0 = [it/L for it in r0] + [it/V for it in l.v0[0]] + [m0/MASS]
        state1 = [it/L for it in r1] + [it/V for it in l.v1[0]] + [m1/MASS]

        for i in range(100):
            veff = np.random.uniform(4000, 8000) * pk.G0
            veff_nd = veff / V
            ta.pars[4] = 1. / veff_nd
            ta_var.pars[4] = 1. / veff_nd
            nseg = int(np.random.uniform(4, 20))
            tgrid = np.linspace(t0*pk.DAY2SEC/TIME, t1*pk.DAY2SEC/TIME, nseg+1)
            controls = np.random.uniform(-1,1, (4*nseg,))
            controls[0::4] /= (F * 100.)  # force will be in [-0.01, 0.01] N
            controls[0::4] = np.abs(controls[0::4]) # force will be in [0.s, 0.01] N
            cut = np.random.uniform(0,1)

            leg = pk.leg.zoh(state0, controls.tolist(), state1, tgrid, cut = cut, tas = [ta, ta_var])
            leg_copy = deepcopy(leg)
            
            dmcdx0, dmcdxf, dmcdcon, dmcdtgrid = leg.compute_mc_grad()
            dtcdcon = leg.compute_tc_grad()
            
            # Check on dmc/dx0
            dmcdx0_n = pg.estimate_gradient(lambda x: compute_mismatch_constraints_n(leg_copy, x, leg_copy.controls, leg_copy.state1, leg_copy.tgrid), leg.state0).reshape(7,-1)
            self.assertTrue(np.linalg.norm(dmcdx0_n-dmcdx0) < 1e-4)
            
            # Check on dmc/dxf
            dmcdxf_n = pg.estimate_gradient(lambda x: compute_mismatch_constraints_n(leg_copy, leg_copy.state0, leg_copy.controls, x, leg_copy.tgrid), leg.state1).reshape(7,-1)
            self.assertTrue(np.linalg.norm(dmcdxf_n-dmcdxf) < 1e-4)

            # Check on dmc/dcontrols
            dmcdcon_n = pg.estimate_gradient(lambda x: compute_mismatch_constraints_n(leg_copy, leg_copy.state0, x, leg_copy.state1, leg_copy.tgrid), leg.controls, dx=1E-8).reshape(7,-1)
            self.assertTrue(np.linalg.norm(dmcdcon_n-dmcdcon) < 1e-4)
            
            # Check on dmc/dtgrid
            dmcdtgrid_n = pg.estimate_gradient(lambda x: compute_mismatch_constraints_n(leg_copy, leg_copy.state0, leg_copy.controls, leg_copy.state1, x), leg.tgrid).reshape(7,-1)
            self.assertTrue(np.linalg.norm(dmcdtgrid_n-dmcdtgrid) < 1e-4)
            
            # Check on dtc/dtcdcon
            dtcdtgrid_n = pg.estimate_gradient(lambda x: compute_throttle_constraints_n(leg_copy, x), leg.controls).reshape(leg.nseg,-1)
            self.assertTrue(np.linalg.norm(dtcdtgrid_n-dtcdcon) < 1e-4)


