try:
	from PyGMO.problem import base as PyGMO_problem

	"""
	This example on the use of PyKEP constructs, using PyGMO for optimization, an interplanetary low-thrust optimization
	problem that can be then solved using one of the available PyGMO solvers. The problem is a non-linear constrained
	problem thas uses the sims-flanagan transcription to model the low-thrust trajectory. PyKEP plots capabilities
	are also demonstrated via the plot method. The interplanetary mission modelled is a lt-mga mission Earth-Venus-Mercury.

	"""
	class mga_lt_EVMe(PyGMO_problem):
		"""
		This constructs a PaGMO.problem object that represents a low-thrust transfer between Earth and Mercury with a Venus
		fly-by. The decision vector contains
		[
		t0,
		T1, mf1, Vxi1, Vyi1, Vzi1, Vxf1, Vyf1, Vzf1,
		T2, mf2, Vxi2, Vyi2, Vzi2, Vxf2, Vyf2, Vzf2,
		[throttles1], [throttles2]
		]
		in the following units: [mjd2000, days, kg, m/s,m/s,m/s, [non-dimensional]]
		"""
		def __init__(self, mass=2000, Tmax=0.5, Isp=3500, Vinf_dep=3, Vinf_arr=2, nseg1=5, nseg2=20):
			#First we call the constructor for the base PyGMO problem
			#(dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
			super(mga_lt_EVMe,self).__init__(17 + 3 * (nseg1+nseg2), 0, 1, 14 + 1 + nseg1+nseg2 + 3, nseg1+nseg2 + 3, 1e-4)

			#We then define some data members (we use the double underscore to indicate they are private)
			from PyKEP import planet_ss, MU_SUN
			from PyKEP.sims_flanagan import spacecraft, leg
			self.__earth = planet_ss('earth')
			self.__venus = planet_ss('venus')
			self.__mercury = planet_ss('mercury')
			self.__sc = spacecraft(mass,Tmax,Isp)
			self.__Vinf_dep = Vinf_dep*1000
			self.__Vinf_arr = Vinf_arr*1000
			self.__leg1 = leg()
			self.__leg2 = leg()
			self.__leg1.set_mu(MU_SUN)
			self.__leg1.set_spacecraft(self.__sc)
			self.__leg2.set_mu(MU_SUN)
			self.__leg2.set_spacecraft(self.__sc)
			self.__nseg1 = nseg1
			self.__nseg2 = nseg2

			#And the box-bouds (launch windows, allowed velocities etc.)
			lb = [3000,100,mass/2] + [-self.__Vinf_dep]*3 + [-6000]*3 + [200,mass/9] + [-6000]*3 + [-self.__Vinf_arr]*3 + [-1,-1,-1] * (nseg1+nseg2)
			ub = [4000,1000,mass]   + [self.__Vinf_dep]*3  + [6000]*3  + [2000,mass]    + [6000]*3  + [ self.__Vinf_arr]*3 + [1,1,1]    * (nseg1+nseg2)
			self.set_bounds(lb,ub)

		#This is the objective function
		def _objfun_impl(self,x):
			return (-x[10],)

		#And these are the constraints
		def _compute_constraints_impl(self,x):
			from PyKEP import epoch, AU, EARTH_VELOCITY, fb_con
			from PyKEP.sims_flanagan import leg, sc_state
			from numpy.linalg import norm
			from math import sqrt, asin, acos

			#Ephemerides
			t_E = epoch(x[0])
			t_V = epoch(x[0] + x[1])
			t_M = epoch(x[0] + x[1] + x[9])
			rE, vE = self.__earth.eph(t_E)
			rV, vV = self.__venus.eph(t_V)
			rM, vM = self.__mercury.eph(t_M)

			#First Leg
			v = [a+b for a,b in zip(vE,x[3:6])]
			x0 = sc_state(rE,v,self.__sc.mass)
			v = [a+b for a,b in zip(vV,x[6:9])]
			xe = sc_state(rV, v ,x[2])
			self.__leg1.set(t_E,x0,x[-3 * (self.__nseg1 + self.__nseg2):-self.__nseg2 * 3],t_V,xe)

			#Second leg
			v = [a+b for a,b in zip(vV,x[11:14])]
			x0 = sc_state(rV,v,x[2])
			v = [a+b for a,b in zip(vM,x[14:17])]
			xe = sc_state(rM, v ,x[10])
			self.__leg2.set(t_E,x0,x[(-3 * self.__nseg2):],t_V,xe)

			#Defining the costraints
			#departure
			v_dep_con = (x[3]  * x[3]  + x[4]  * x[4]  + x[5]  * x[5]  - self.__Vinf_dep * self.__Vinf_dep) / (EARTH_VELOCITY * EARTH_VELOCITY)
			#arrival
			v_arr_con = (x[14] * x[14] + x[15] * x[15] + x[16] * x[16] - self.__Vinf_arr * self.__Vinf_arr) / (EARTH_VELOCITY * EARTH_VELOCITY)
			#fly-by at Venus
			DV_eq, alpha_ineq = fb_con(x[6:9],x[11:14],self.__venus)

			#Assembling the constraints
			retval = list(self.__leg1.mismatch_constraints() + self.__leg2.mismatch_constraints()) + [DV_eq] + list(self.__leg1.throttles_constraints() + self.__leg2.throttles_constraints()) + [v_dep_con] + [v_arr_con] + [alpha_ineq]

			#We then scale all constraints to non dimensional values
			#leg 1
			retval[0] /= AU
			retval[1] /= AU
			retval[2] /= AU
			retval[3] /= EARTH_VELOCITY
			retval[4] /= EARTH_VELOCITY
			retval[5] /= EARTH_VELOCITY
			retval[6] /= self.__sc.mass
			#leg 2
			retval[7] /= AU
			retval[8] /= AU
			retval[9] /= AU
			retval[10] /= EARTH_VELOCITY
			retval[11] /= EARTH_VELOCITY
			retval[12] /= EARTH_VELOCITY
			retval[13] /= self.__sc.mass
			#fly-by at Venus
			retval[14] /= (EARTH_VELOCITY*EARTH_VELOCITY)

			return retval

		#This transforms the leg into a high fidelity one
		def high_fidelity(self,boolean):
			self.__leg.high_fidelity = boolean

		#And this help visualizing the trajectory
		def plot(self,x):
			import matplotlib as mpl
			from mpl_toolkits.mplot3d import Axes3D
			import matplotlib.pyplot as plt
			from PyKEP import epoch, AU
			from PyKEP.sims_flanagan import sc_state
			from PyKEP.orbit_plots import plot_planet, plot_sf_leg

			t_E = epoch(x[0])
			t_V = epoch(x[0] + x[1])
			t_M = epoch(x[0] + x[1] + x[9])
			rE, vE = self.__earth.eph(t_E)
			rV, vV = self.__venus.eph(t_V)
			rM, vM = self.__mercury.eph(t_M)

			#First Leg
			v = [a+b for a,b in zip(vE,x[3:6])]
			x0 = sc_state(rE,v,self.__sc.mass)
			v = [a+b for a,b in zip(vV,x[6:9])]
			xe = sc_state(rV, v ,x[2])
			self.__leg1.set(t_E,x0,x[-3 * (self.__nseg1 + self.__nseg2):-self.__nseg2 * 3],t_V,xe)

			#Second leg
			v = [a+b for a,b in zip(vV,x[11:14])]
			x0 = sc_state(rV,v,x[2])
			v = [a+b for a,b in zip(vM,x[14:17])]
			xe = sc_state(rM, v ,x[10])
			self.__leg2.set(t_E,x0,x[(-3 * self.__nseg2):],t_V,xe)

			fig = plt.figure()
			ax = fig.gca(projection='3d')

			#The Sun
			ax.scatter(0,0,0, color='y')
			#The legs
			plot_sf_leg(ax, self.__leg1, units=AU,N=10)
			plot_sf_leg(ax, self.__leg2, units=AU,N=10)
			#The planets
			plot_planet(ax, self.__earth, t_E, units=AU, legend = True,color=(0.7,0.7,1))
			plot_planet(ax, self.__venus, t_V, units=AU, legend = True,color=(0.7,0.7,1))
			plot_planet(ax, self.__mercury, t_M, units=AU, legend = True,color=(0.7,0.7,1))
			plt.show()

	def run_example3():
		from PyGMO import algorithm, island
		prob = mga_lt_EVMe()
		algo = algorithm.scipy_slsqp(max_iter = 500, acc =1e-5)
		#algo = algorithm.snopt(major_iter=2000, opt_tol=1e-3, feas_tol=1e-9)
		algo2 = algorithm.mbh(algo,5,0.05)
		algo2.screen_output = True
		isl = island(algo2,prob,1)
		print "Running Monotonic Basin Hopping ...."
		isl.evolve(1); isl.join()
		print "Is the solution found a feasible trajectory? " + str(prob.feasibility_x(isl.population.champion.x))
		prob.plot(isl.population.champion.x)

except:
	print "Could not import PyGMO. PyGMO is required for some PyKEP examples"
