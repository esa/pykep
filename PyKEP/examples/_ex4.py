try:
	from PyGMO.problem import base as PyGMO_problem
	"""
	This example demonstarates the use of the interplanetary leg in the Sundmann variable to obtain automated mesh optimization
	"""
	class mga_lt_earth_mars_sundmann(PyGMO_problem):
		"""
		This constructs a PaGMO.problem object that represents a low-thrust transfer between Earth and Mars. 
		The interplantary leg uses the Sundmann variable. The decision vector
		contains [t0,T,sf,mf,Vx,Vy,Vz,[throttles]] in the following units: [mjd2000, days, days, kg, m/s,m/s,m/s, [non-dimensional]]
		"""
		def __init__(self,mass=1000,Tmax=0.1,Isp=2500,Vinf=3.0,nseg=20):
			#First we call the constructor for the base PyGMO problem 
			#(dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
			super(mga_lt_earth_mars_sundmann,self).__init__(7 + nseg*3,0,1,9 + nseg,nseg+1,1e-4)
			
			#We then define some data members (we use the double underscore to indicate they are private)
			from PyKEP import planet_ss, MU_SUN, AU, planet_gtoc5
			from PyKEP.sims_flanagan import spacecraft, leg_s
			self.__earth = planet_ss('earth')
			self.__mars = planet_ss('jupiter')
			self.__sc = spacecraft(mass,Tmax,Isp)
			self.__Vinf = Vinf*1000
			#here we construct the trajectory leg in the Sundmann variable t = (r/10AU)^1.5 s
			self.__leg = leg_s(nseg,1.0/(100*AU)**1.0,1.0)
			self.__leg.set_mu(MU_SUN)
			self.__leg.set_spacecraft(self.__sc)
			self.__nseg = nseg
			#The bounds on the Sundmann variable can be evaluated considering circular orbits at r=1AUand r=0.7AU (for example)
			self.set_bounds([5000,2400,10000, self.__sc.mass/10,-self.__Vinf,-self.__Vinf,-self.__Vinf] + [-1] * 3 *nseg,[8000,2500,150000,self.__sc.mass,self.__Vinf,self.__Vinf,self.__Vinf] + [1] * 3 * nseg)

		#This is the objective function
		def _objfun_impl(self,x):
			return (-x[3],)
			
		#And these are the constraints
		def _compute_constraints_impl(self,x):
			from PyKEP import epoch, AU, EARTH_VELOCITY, DAY2SEC
			from PyKEP.sims_flanagan import sc_state
			
			start = epoch(x[0])
			end = epoch(x[0] + x[1])
			
			r,v = self.__earth.eph(start)
			v = [a+b for a,b in zip(v,x[4:7])]
			x0 = sc_state(r,v,self.__sc.mass)
			
			r,v = self.__mars.eph(end)
			xe = sc_state(r, v ,x[3])
			self.__leg.set(start,x0,x[-3 * self.__nseg:],end,xe, x[2] * DAY2SEC)
			v_inf_con = (x[4] * x[4] + x[5] * x[5] + x[6] * x[6] - self.__Vinf * self.__Vinf) / (EARTH_VELOCITY * EARTH_VELOCITY)
			try: 
				retval = list(self.__leg.mismatch_constraints() + self.__leg.throttles_constraints()) + [v_inf_con]
			except:
				print "warning: CANNOT EVALUATE constraints .... possible problem in the taylor integration in the Sundmann variable"
				return (1e14,)*(8+1+self.__nseg+2)
			#We then scale all constraints to non dimensiona values
			retval[0] /= AU
			retval[1] /= AU
			retval[2] /= AU
			retval[3] /= EARTH_VELOCITY
			retval[4] /= EARTH_VELOCITY
			retval[5] /= EARTH_VELOCITY
			retval[6] /= self.__sc.mass
			retval[7] /= 365.25 * DAY2SEC
			return retval
			
			
		#And this help visualizing the trajectory
		def plot(self,x):
			import matplotlib as mpl
			from mpl_toolkits.mplot3d import Axes3D
			import matplotlib.pyplot as plt
			from PyKEP import epoch, AU, DAY2SEC
			from PyKEP.sims_flanagan import sc_state
			from PyKEP.orbit_plots import plot_planet, plot_leg
			
			start = epoch(x[0])
			end = epoch(x[0] + x[1])
			
			r,v = self.__earth.eph(start)
			v = [a+b for a,b in zip(v,x[4:7])]
			x0 = sc_state(r,v,self.__sc.mass)
			
			r,v = self.__mars.eph(end)
			xe = sc_state(r, v ,x[3])
			self.__leg.set(start,x0,x[-3 * self.__nseg:],end,xe, x[2] * DAY2SEC)
			
			fig = plt.figure()
			ax = fig.gca(projection='3d')
			#The Sun
			ax.scatter([0],[0],[0], color='y')
			#The leg
			plot_leg(ax, self.__leg, units=AU,N=10)
			#The planets
			plot_planet(ax, self.__earth, start, units=AU, legend = True,color=(0.8,0.8,1))
			plot_planet(ax, self.__mars, end, units=AU, legend = True,color=(0.8,0.8,1))
			plt.show()
			
	def run_example4():
		from PyGMO import algorithm, island,population
		N=20
		prob = mga_lt_earth_mars_sundmann(nseg=N)
		algo = algorithm.scipy_slsqp(max_iter = 500, acc=1e-5)
		#algo = algorithm.snopt(major_iter=1000, opt_tol=1e-6, feas_tol=1e-11)
		algo2 = algorithm.mbh(algo,5,0.05)
		algo2.screen_output = True
		isl = island(algo2,prob,1)
		print "Running Monotonic Basin Hopping ...."
		isl.evolve(1); isl.join()
		print "Is the solution found a feasible trajectory? " + str(prob.feasibility_x(isl.population.champion.x))
		print isl.population.champion.x
		prob.plot(isl.population.champion.x)

except:
	print "Could not import PyGMO. PyGMO is required for some PyKEP examples"
