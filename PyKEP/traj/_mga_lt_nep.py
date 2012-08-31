from PyGMO.problem import base as base_problem
from PyKEP import planet_ss,epoch, fb_con, EARTH_VELOCITY, AU, MU_SUN
from PyKEP.sims_flanagan import leg, spacecraft, sc_state

"""
This class represents, as a global optimization problem (linearly constrained, high diemensional), a low-thrust
interplanetary trajectory modelled as a Multiple Gravity Assist trajectory with legs represented by an impulsive transcription...

The decision vector (chromosome) is:
[t0] + 
[T1, mf1, Vxi1, Vyi1, Vzi1, Vxf1, Vyf1, Vzf1] +
[T2, mf2, Vxi2, Vyi2, Vzi2, Vxf2, Vyf2, Vzf2] + .....
[throttles1] + [throttles2] + ....

SEE : Yam, C.H., di Lorenzo, D., and Izzo, D.,	 Low-Thrust Trajectory Design as a Constrained Global Optimization Problem, 
Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering, 225(11), pp.1243-1251, 2011.
"""

class mga_lt_nep(base_problem):
	"""
	Constructs an mga_lt problem
	
	USAGE: opt_prob = mga_lt_nep(seq = [planet_ss('earth'),planet_ss('venus'),planet_ss('earth')], nseg = [10]*2, 
	t0 = [epoch(0),epoch(1000)], T = [[200,500],[200,500]], Vinf_dep=2.5, Vinf_arr=2.0, mass=4000.0, Tmax=1.0, Isp=2000.0,
	multi_objective = False, fb_rel_vel = 6)

	* seq: list of PyKEP.planet defining the encounter sequence for the trajectoty (including the initial planet)
	* nseg: list of integers containing the number of segments to be used for each leg (len(nseg) = len(seq)-1)
	* t0: list of PyKEP epochs defining the launch window
	* T: minimum and maximum time of each leg (days)
	* vinf_dep: maximum launch hyperbolic velocity allowed (in km/sec)
	* vinf_arr: maximum arrival hyperbolic velocity allowed (in km/sec)
	* mass: spacecraft starting mass
	* Tmax: maximum thrust
	* Isp: engine specific impulse
	* multi-objective: when True defines the problem as a multi-objective problem, returning total DV and time of flight
	* fb_rel_vel = determines the bounds on the maximum allowed relative velocity at all fly-bys (in km/sec)
	"""
	def __init__(self, seq = [planet_ss('earth'),planet_ss('venus'),planet_ss('earth')], nseg = [10]*2, 
	t0 = [epoch(0),epoch(1000)], T = [[200,500],[200,500]], Vinf_dep=2.5, Vinf_arr=2.0, mass=4000.0, Tmax=1.0, Isp=2000.0,
	multi_objective = False, fb_rel_vel = 6):
		#1) We comute the problem dimensions .... and call the base problem constructor
		self.__n_legs = len(seq) - 1
		n_fb = self.__n_legs - 1
		# 1a) The decision vector length
		dim = 1 + self.__n_legs * 8 + sum(nseg) * 3
		# 1b) The total number of constraints (mismatch + fly-by + boundary + throttles
		c_dim = self.__n_legs * 7 + n_fb * 2 + 2 + sum(nseg)
		# 1c) The number of inequality constraints (boundary + fly-by angle + throttles)
		c_ineq_dim = 2 + n_fb + sum(nseg)
		# 1d) the number of objectives
		f_dim = multi_objective + 1
		#First we call the constructor for the base PyGMO problem 
		#As our problem is n dimensional, box-bounded (may be multi-objective), we write
		#(dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
		super(mga_lt_nep,self).__init__(dim,0,f_dim,c_dim,c_ineq_dim,1e-4)

		#2) We then define some class data members
		#public:
		self.seq = seq
		#private:
		self.__nseg = nseg
		self.__Vinf_dep = Vinf_dep*1000
		self.__Vinf_arr = Vinf_arr*1000
		self.__sc = spacecraft(mass,Tmax,Isp)
		self.__leg = leg()
		self.__leg.set_mu(MU_SUN)
		self.__leg.set_spacecraft(self.__sc)
		fb_rel_vel*=1000
		#3) We compute the bounds
		lb = [t0[0].mjd2000] + [0, mass / 2, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel] * self.__n_legs + [-1,-1,-1] * sum(self.__nseg)
		ub = [t0[1].mjd2000] + [1, mass, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel] * self.__n_legs + [1,1,1] * sum(self.__nseg)
		#3a ... and account for the bounds on the vinfs......
		lb[3:6] = [-self.__Vinf_dep]*3
		ub[3:6] = [self.__Vinf_dep]*3
		lb[-sum(self.__nseg)*3-3:-sum(self.__nseg)*3] = [-self.__Vinf_arr]*3
		ub[-sum(self.__nseg)*3-3:-sum(self.__nseg)*3] = [self.__Vinf_arr]*3
		# 3b... and for the time of flight
		lb[1:1+8*self.__n_legs:8] = [el[0] for el in T]
		ub[1:1+8*self.__n_legs:8] = [el[1] for el in T]
		
		#4) And we set the bounds
		self.set_bounds(lb,ub)

	#Objective function
	def _objfun_impl(self,x):
		if self.f_dimension == 1:
			return (-x[2 + (self.__n_legs - 1) * 8],)
		else:
			return (-x[2 + (self.__n_legs - 1) * 8], sum(x[1:1+8*self.__n_legs:8]))		

	#Constraints function
	def _compute_constraints_impl(self,x):
		# 1 - We decode the chromosome extracting the time of flights
		T = list([0]*(self.__n_legs))
		for i in range(self.__n_legs):
			T[i] = x[1+i*8]
			
		#2 - We compute the epochs and ephemerides of the planetary encounters
		t_P = list([None] * (self.__n_legs+1))
		r_P = list([None] * (self.__n_legs+1))
		v_P = list([None] * (self.__n_legs+1))
		
		for i,planet in enumerate(self.seq):
			t_P[i] = epoch(x[0] + sum(T[0:i]))
			r_P[i],v_P[i] = self.seq[i].eph(t_P[i])
		
		#3 - We iterate through legs to compute mismatches and throttles constraints
		ceq = list()
		cineq = list()
		m0 = self.__sc.mass
		for i in range(self.__n_legs):
			#First Leg
			v = [a+b for a,b in zip(v_P[i],x[(3 + i * 8):(6 + i * 8)])]
			x0 = sc_state(r_P[i],v,m0)
			v = [a+b for a,b in zip(v_P[i+1],x[(6 + i * 8):(11 + i * 8)])]
			xe = sc_state(r_P[i+1], v ,x[2 + i * 8])
			throttles = x[(1 + 8*self.__n_legs + 3*sum(self.__nseg[:i])):(1 + 8*self.__n_legs + 3*sum(self.__nseg[:i])+3*self.__nseg[i])]
			self.__leg.set(t_P[i],x0,throttles,t_P[i+1],xe)
			#update mass!
			m0 = x[2+8*i]
			ceq.extend(self.__leg.mismatch_constraints())
			cineq.extend(self.__leg.throttles_constraints())

		#Adding the boundary constraints
		#departure
		v_dep_con = (x[3] ** 2  + x[4] ** 2  + x[5] **2  - self.__Vinf_dep ** 2) / (EARTH_VELOCITY**2)
		#arrival
		v_arr_con = (x[6 + (self.__n_legs-1)*8]**2 + x[7+ (self.__n_legs-1)*8]**2 + x[8+ (self.__n_legs-1)*8]**2 - self.__Vinf_arr ** 2) / (EARTH_VELOCITY**2)
		cineq.append(v_dep_con*100)
		cineq.append(v_arr_con*100)

		#We add the fly-by constraints
		for i in range(self.__n_legs-1):
			DV_eq, alpha_ineq = fb_con(x[6 + i*8:9 + i*8],x[11+ i * 8:14+ i * 8],self.seq[i+1])
			ceq.append(DV_eq / (EARTH_VELOCITY**2))
			cineq.append(alpha_ineq)

		#Making the mismatches non dimensional
		for i in range(self.__n_legs):
			ceq[0+i*7] /= AU
			ceq[1+i*7] /= AU
			ceq[2+i*7] /= AU
			ceq[3+i*7] /= EARTH_VELOCITY
			ceq[4+i*7] /= EARTH_VELOCITY
			ceq[5+i*7] /= EARTH_VELOCITY
			ceq[6+i*7] /= 1000

		#We assemble the constraint vector
		retval = list()
		retval.extend(ceq)
		retval.extend(cineq)

		return retval
	
	#And this helps visualizing the trajectory
	def plot(self,x):
		import matplotlib as mpl
		from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
		from PyKEP import epoch, AU
		from PyKEP.sims_flanagan import sc_state
		from PyKEP.orbit_plots import plot_planet, plot_sf_leg
		
		fig = plt.figure()
		ax = fig.gca(projection='3d')

		#Plotting the Sun ........
		ax.scatter([0],[0],[0], color='y')
		
		#Plotting the legs .......
		
		# 1 - We decode the chromosome extracting the time of flights
		T = list([0]*(self.__n_legs))
		for i in range(self.__n_legs):
			T[i] = x[1+i*8]
			
		#2 - We compute the epochs and ephemerides of the planetary encounters
		t_P = list([None] * (self.__n_legs+1))
		r_P = list([None] * (self.__n_legs+1))
		v_P = list([None] * (self.__n_legs+1))
		
		for i,planet in enumerate(self.seq):
			t_P[i] = epoch(x[0] + sum(T[0:i]))
			r_P[i],v_P[i] = self.seq[i].eph(t_P[i])
		
		#3 - We iterate through legs to compute mismatches and throttles constraints
		ceq = list()
		cineq = list()
		m0 = self.__sc.mass
		for i in range(self.__n_legs):
			#First Leg
			v = [a+b for a,b in zip(v_P[i],x[(3 + i * 8):(6 + i * 8)])]
			x0 = sc_state(r_P[i],v,m0)
			v = [a+b for a,b in zip(v_P[i+1],x[(6 + i * 8):(11 + i * 8)])]
			xe = sc_state(r_P[i+1], v ,x[2 + i * 8])
			throttles = x[(1 +8 * self.__n_legs + 3*sum(self.__nseg[:i])):(1 +8 * self.__n_legs + 3*sum(self.__nseg[:i])+3*self.__nseg[i])]
			self.__leg.set(t_P[i],x0,throttles,t_P[i+1],xe)
			#update mass!
			m0 = x[2+8*i]
			plot_sf_leg(ax, self.__leg, units=AU,N=10)
			
		#Plotting planets
		for i,planet in enumerate(self.seq):
			plot_planet(ax, planet, t_P[i], units=AU, legend = True,color=(0.7,0.7,1))

		plt.show()
			
