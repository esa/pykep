from PyGMO.problem import base

class traj_tutorial_01(base):
	def __init__(self):
		super(traj_tutorial_01,self).__init__(36,0,1,18,11,1e-7)
		import PyKEP
		self.__earth = PyKEP.planet_ss('earth')
		self.__mars = PyKEP.planet_ss('mars')
		self.__sc = PyKEP.sims_flanagan.spacecraft(1000,0.05,2500)
		self.set_bounds([0,60,500,-4000,-4000,-4000] + [-1] * 30,[3000,1000,self.__sc.mass,4000,4000,4000] + [1] * 30)
	def _objfun_impl(self,x):
		return (-x[2],)
	def _compute_constraints_impl(self,x):
		import PyKEP
		start = PyKEP.epoch(x[0])
		end = PyKEP.epoch(x[0] + x[1])
		r,v = self.__earth.eph(start)
		v_list = list(v)
		v_list[0] += x[3]
		v_list[1] += x[4]
		v_list[2] += x[5]
		x0 = PyKEP.sims_flanagan.sc_state(r,v_list,self.__sc.mass)
		r,v = self.__mars.eph(end)
		xe = PyKEP.sims_flanagan.sc_state(r, v ,x[2])
		l = PyKEP.sims_flanagan.leg(start,x0,x[-30:],end,xe,self.__sc,PyKEP.MU_SUN)
		v_inf_con = (x[3] * x[3] + x[4] * x[4] + x[5] * x[5] - 1.6E7) / (PyKEP.EARTH_VELOCITY * PyKEP.EARTH_VELOCITY)
		retval = list(l.mismatch_constraints() + l.throttles_constraints()) + [v_inf_con]
		retval[0] /= PyKEP.AU
		retval[1] /= PyKEP.AU
		retval[2] /= PyKEP.AU
		retval[3] /= PyKEP.EARTH_VELOCITY
		retval[4] /= PyKEP.EARTH_VELOCITY
		retval[5] /= PyKEP.EARTH_VELOCITY
		retval[6] /= self.__sc.mass
		return retval
