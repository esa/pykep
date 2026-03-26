import pykep as pk
import numpy as np

import unittest as _ut

def compute_numerical_gradient(sf_leg, sf_leg_type = 'lf'):
    import numpy as np
    import pykep as pk
    import pygmo as pg

    state_length = np.array(sf_leg.rvs).flatten().size + 1
    throttle_length = np.array(sf_leg.throttles).size
    chromosome = np.zeros((state_length * 2 + throttle_length + 1))
    chromosome[0:state_length] = np.append(np.array(sf_leg.rvs).flatten(), sf_leg.ms)
    chromosome[state_length:state_length+throttle_length] = np.array(sf_leg.throttles)
    chromosome[state_length+throttle_length:state_length*2+throttle_length] = np.append(np.array(sf_leg.rvf).flatten(), sf_leg.mf)
    chromosome[-1] = sf_leg.tof

    def set_and_compute_constraints(chromosome, sf_leg_type = 'lf'):

        if sf_leg_type == 'hf' or sf_leg_type == 'high-fidelity':
            sf_leg_constraint = pk.leg.sims_flanagan_hf()
        else:
            sf_leg_constraint = pk.leg.sims_flanagan()
        sf_leg_constraint.cut = 0.5
        sf_leg_constraint.max_thrust = 1
        sf_leg_constraint.mu = 1
        sf_leg_constraint.veff = 1
        sf_leg_constraint.rvs = [chromosome[0:3],chromosome[3:6]]
        sf_leg_constraint.ms = chromosome[6]
        sf_leg_constraint.throttles = chromosome[state_length:state_length+throttle_length]
        sf_leg_constraint.rvf = [chromosome[state_length+throttle_length:state_length+throttle_length+3],chromosome[state_length+throttle_length+3:state_length+throttle_length+6]]
        sf_leg_constraint.mf = chromosome[2*state_length+throttle_length-1]
        sf_leg_constraint.tof = chromosome[2*state_length+throttle_length]
        eq_con = sf_leg_constraint.compute_mismatch_constraints()
        ineq_con = sf_leg_constraint.compute_throttle_constraints()
        return np.concatenate((eq_con, ineq_con))

    return pg.estimate_gradient_h(callable = lambda x : set_and_compute_constraints(x, sf_leg_type), x=chromosome)


class leg_sims_flanagan_test(_ut.TestCase):

    def test_sims_flanagan(self):
        import numpy as np
        import pykep as pk

        udpla_e = pk.udpla.vsop2013("earth_moon", 1e-2)
        udpla_j = pk.udpla.vsop2013("jupiter", 1e-2)
        earth = pk.planet(udpla_e)
        jupiter = pk.planet(udpla_j)
        dt_days = 1000
        dt = dt_days * pk.DAY2SEC
        t0 = 1233.3
        rv0 = earth.eph(t0)
        rv1 = jupiter.eph(t0 + dt_days)
        lp = pk.lambert_problem(rv0[0], rv1[0], dt, pk.MU_SUN)
        rv0[1] = lp.v0[0]
        rv1[1] = lp.v1[0]

        cut_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        mc_list = []
        for i in range(1, 34):
            for cut in cut_values:
                throttles = [0.0] * i * 3
                sf_hf_leg = pk.leg.sims_flanagan(rv0, 1.0, throttles, rv1, 1.0, dt, 1.0, 1.0, pk.MU_SUN, cut)
                mc = sf_hf_leg.compute_mismatch_constraints()
                mc[0] /= pk.AU
                mc[1] /= pk.AU
                mc[2] /= pk.AU
                mc[3] /= pk.EARTH_VELOCITY
                mc[4] /= pk.EARTH_VELOCITY
                mc[5] /= pk.EARTH_VELOCITY
                mc[6] /= 1000
                mc_list.append(mc)
                
        self.assertTrue(np.array([np.max(i) < 1e-8 for i in mc_list]).all())

    def test_mc_grad(self):
        import numpy as np
        import pykep as pk

        sf_leg = pk.leg.sims_flanagan()
        sf_leg.cut = 0.5
        sf_leg.throttles = np.array([0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24,
        0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34])
        sf_leg.rvs = np.array([[1, 0.1, -0.1], [0.2, 1.0, -0.2]])
        sf_leg.ms = 1
        sf_leg.rvf = np.array([[1.2, -0.1, 0.1], [-0.2, 1.023, -0.44]])
        sf_leg.mf = 13 / 15
        sf_leg.max_thrust = 1
        sf_leg.mu = 1
        sf_leg.veff = 1
        sf_leg.tof = 1
        #sf_leg.tol = 1e-16
        state_length = np.array(sf_leg.rvs).flatten().size + 1
        throttle_length = np.array(sf_leg.throttles).size

        num_grad = compute_numerical_gradient(sf_leg, sf_leg_type = 'lf')
        num_grad = num_grad.reshape((17, 45), order='C')
        grad_rvm, grad_rvm_bck, grad_final = sf_leg.compute_mc_grad()
        a_tc_grad = sf_leg.compute_tc_grad()
        a_grad = np.zeros((state_length+throttle_length // 3, 2 * state_length + throttle_length + 1))
        a_grad[0:state_length, 0:state_length] = grad_rvm
        a_grad[0:state_length, state_length:state_length + throttle_length] = grad_final[:,0:throttle_length] 
        a_grad[0:state_length, state_length+throttle_length:state_length*2+throttle_length] = grad_rvm_bck
        a_grad[0:state_length, state_length*2+throttle_length] = grad_final[:, throttle_length:throttle_length + 1].reshape(7,)
        a_grad[state_length:, state_length:state_length+throttle_length] = a_tc_grad
        self.assertTrue(np.allclose(num_grad, a_grad, atol=1e-8))

    def test_pickling(self):
        import pickle
        import io
        import pykep as pk

        # An example object
        data = pk.leg.sims_flanagan()
        data.cut = 0.5
        data.throttles = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        data.rvs = np.array([[1, 0.1, -0.1], [0.2, 1.0, -0.2]])
        data.ms = 1
        data.rvf = np.array([[1.2, -0.1, 0.1], [-0.2, 1.023, -0.44]])
        data.mf = 13 / 15
        data.max_thrust = 1
        data.mu = 1
        data.veff = 1
        data.tof = 1

        # Create an in-memory bytes buffer
        buffer = io.BytesIO()

        # Pickle the object into the buffer
        pickle.dump(data, buffer)

        # Reset buffer position to the start for reading
        buffer.seek(0)

        # Unpickle the data from the buffer
        loaded_data = pickle.load(buffer)

        # Verify that the original and loaded data are the same
        self.assertTrue(loaded_data.__repr__()==data.__repr__())

class sims_flanagan_hf_test(_ut.TestCase):
    def test_comparison_sf_and_sf_hf(self):
        import pykep as pk
        import numpy as np

        sf_leg = pk.leg.sims_flanagan()
        sf_leg.cut = 0.5
        sf_leg.throttles = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        sf_leg.rvs = np.array([[1, 0.1, -0.1], [0.2, 1.0, -0.2]])
        sf_leg.ms = 1
        sf_leg.rvf = np.array([[1.2, -0.1, 0.1], [-0.2, 1.023, -0.44]])
        sf_leg.mf = 13 / 15
        sf_leg.max_thrust = 1
        sf_leg.mu = 1
        sf_leg.veff = 1
        sf_leg.tof = 1
        rvm_mc_sf = sf_leg.compute_mismatch_constraints()

        sf_hf_leg = pk.leg.sims_flanagan_hf()
        sf_hf_leg.cut = 0.5
        sf_hf_leg.throttles = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        sf_hf_leg.rvms = np.array([1, 0.1, -0.1, 0.2, 1.0, -0.2, 1])
        sf_hf_leg.rvmf = np.array([1.2, -0.1, 0.1, -0.2, 1.023, -0.44, 13 / 15])
        sf_hf_leg.max_thrust = 1
        sf_hf_leg.mu = 1
        sf_hf_leg.veff = 1
        sf_hf_leg.tof = 1
        rvm_mc_sf_hf = sf_hf_leg.compute_mismatch_constraints()
        self.assertTrue(np.allclose(rvm_mc_sf, rvm_mc_sf_hf, atol=1e-13))

    def test_sims_flanagan_hf(self):
        import numpy as np
        import pykep as pk

        udpla_e = pk.udpla.vsop2013("earth_moon", 1e-2)
        udpla_j = pk.udpla.vsop2013("jupiter", 1e-2)
        earth = pk.planet(udpla_e)
        jupiter = pk.planet(udpla_j)
        dt_days = 1000
        dt = dt_days * pk.DAY2SEC
        t0 = 1233.3
        rv0 = earth.eph(t0)
        rv1 = jupiter.eph(t0 + dt_days)
        lp = pk.lambert_problem(rv0[0], rv1[0], dt, pk.MU_SUN)
        rv0[1] = lp.v0[0]
        rv1[1] = lp.v1[0]


        cut_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        mc_list = []
        for i in range(1, 34):
            for cut in cut_values:
                throttles = [0.0] * i * 3
                sf_hf_leg = pk.leg.sims_flanagan_hf(rv0, 1.0, throttles, rv1, 1.0, dt, 1.0, 1.0, pk.MU_SUN, cut)
                mc = sf_hf_leg.compute_mismatch_constraints()
                mc[0] /= pk.AU
                mc[1] /= pk.AU
                mc[2] /= pk.AU
                mc[3] /= pk.EARTH_VELOCITY
                mc[4] /= pk.EARTH_VELOCITY
                mc[5] /= pk.EARTH_VELOCITY
                mc[6] /= 1000
                mc_list.append(mc)
                
        self.assertTrue(np.array([np.max(i) < 1e-8 for i in mc_list]).all())
        
    def test_mc_grad_hf(self):
        import numpy as np
        import pykep as pk

        sf_leg = pk.leg.sims_flanagan_hf()
        sf_leg.cut = 0.5
        sf_leg.throttles = np.array([0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24,
        0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34])
        sf_leg.rvs = np.array([[1, 0.1, -0.1], [0.2, 1.0, -0.2]])
        sf_leg.ms = 1
        sf_leg.rvf = np.array([[1.2, -0.1, 0.1], [-0.2, 1.023, -0.44]])
        sf_leg.mf = 13 / 15
        sf_leg.max_thrust = 1
        sf_leg.mu = 1
        sf_leg.veff = 1
        sf_leg.tof = 1
        state_length = np.array(sf_leg.rvs).flatten().size + 1
        throttle_length = np.array(sf_leg.throttles).size

        num_grad = compute_numerical_gradient(sf_leg, sf_leg_type = 'hf')
        num_grad = num_grad.reshape((17, 45), order='C')
        grad_rvm, grad_rvm_bck, grad_final = sf_leg.compute_mc_grad()
        a_tc_grad = sf_leg.compute_tc_grad()
        a_grad = np.zeros((state_length+throttle_length // 3, 2 * state_length + throttle_length + 1))
        a_grad[0:state_length, 0:state_length] = grad_rvm
        a_grad[0:state_length, state_length:state_length + throttle_length] = grad_final[:,0:throttle_length] 
        a_grad[0:state_length, state_length+throttle_length:state_length*2+throttle_length] = grad_rvm_bck
        a_grad[0:state_length, state_length*2+throttle_length] = grad_final[:, throttle_length:throttle_length + 1].reshape(7,)
        a_grad[state_length:, state_length:state_length+throttle_length] = a_tc_grad
        self.assertTrue(np.allclose(num_grad, a_grad, atol=1e-8))

    def test_pickling(self):
        import pickle
        import io
        import pykep as pk

        # An example object
        data = pk.leg.sims_flanagan_hf()
        data.cut = 0.5
        data.throttles = np.array([0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24,
        0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34])
        data.rvs = np.array([[1, 0.1, -0.1], [0.2, 1.0, -0.2]])
        data.ms = 1
        data.rvf = np.array([[1.2, -0.1, 0.1], [-0.2, 1.023, -0.44]])
        data.mf = 13 / 15
        data.max_thrust = 1
        data.mu = 1
        data.veff = 1
        data.tof = 1

        # Create an in-memory bytes buffer
        buffer = io.BytesIO()

        # Pickle the object into the buffer
        pickle.dump(data, buffer)

        # Reset buffer position to the start for reading
        buffer.seek(0)

        # Unpickle the data from the buffer
        loaded_data = pickle.load(buffer)

        # Verify that the original and loaded data are the same
        self.assertTrue(loaded_data.__repr__()==data.__repr__())
        