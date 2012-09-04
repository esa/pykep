from PyGMO.problem import base as base_problem
from PyKEP import epoch,DAY2SEC,planet_ss,MU_SUN,lambert_problem,propagate_lagrangian,fb_prop, AU
from math import pi, cos, sin, acos
from scipy.linalg import norm

class mga_1dsm(base_problem):
	"""
	This class represents a global optimization problem (box-bounded, continuous) relative to an interplanetary trajectory modelled
	as a Multiple Gravity Assist trajectory that allows one only Deep Space Manouvre between each leg.

	SEE : Izzo: "Global Optimization and Space Pruning for Spacecraft Trajectory Design, Spacecraft Trajectory Optimization, Conway, B. (Eds.), Cambridge University Press, pp.178-199, 2010)

	The decision vector is [t0,u,v,Vinf,eta1,T] + [beta, rp/rV, eta2,a2] ..... in the units: [mjd2000,nd,nd,km/s,nd,years] + [rad,nd,nd,nd] + ....
	where Vinf = Vinf_mag*(cos(theta)*cos(phi)i+cos(theta)*sin(phi)j+sin(phi)k) and theta = 2*pi*u and phi = acos(2*v-1)-pi/2

	Each leg time-of-flight can be obtained as Tn = T*an, T(n-1) = (T - Tn)*a(n-1), .... , Ti = (T-T(i+1)-T(i+2)- .... - Tn)*ai

	NOTE: The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.

	"""
	def __init__(self, seq = [planet_ss('earth'),planet_ss('venus'),planet_ss('earth')], t0 = [epoch(0),epoch(1000)], tof = [1.0,5.0], vinf = 2.5, multi_objective = False):
		"""
		prob = mga_1dsm(seq = [planet_ss('earth'),planet_ss('venus'),planet_ss('earth')], t0 = [epoch(0),epoch(1000)], tof = [1.0,5.0], vinf = 2.5, multi_objective = False)

		* seq:  list of PyKEP.planet defining the encounter sequence for the trajectoty (including the initial planet)
		* t0:   list of PyKEP epochs defining the launch window
		* tof:  minimum and maximum time of flight allowed (in years)
		* vinf: maximum launch hyperbolic velocity allowed
		* multi-objective: when True defines the problem as a multi-objective problem, minimizing total DV and time of flight
		"""
		self.__n_legs = len(seq) - 1
		dim = 6 + (self.__n_legs-1) * 4
		obj_dim = multi_objective + 1
		#First we call the constructor for the base PyGMO problem 
		#As our problem is n dimensional, box-bounded (may be multi-objective), we write
		#(dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
		super(mga_1dsm,self).__init__(dim,0,obj_dim,0,0,0)

		#We then define all planets in the sequence as data members
		self.seq = seq
		
		#And we compute the bounds
		lb = [t0[0].mjd2000,0.0,0.0,0.0      ,1e-5    ,tof[0]*365.25] + [0   ,1.1 ,1e-5    ,1e-5]     * (self.__n_legs-1)
		ub = [t0[1].mjd2000,1.0,1.0,vinf*1000,1.0-1e-5,tof[1]*365.25] + [2*pi,30.0,1.0-1e-5,1.0-1e-5] * (self.__n_legs-1)
		
		#Accounting that each planet has a different safe radius......
		for i,pl in enumerate(seq[1:-1]):
			lb[7+4*i] = pl.safe_radius / pl.radius
			
		#And we set them
		self.set_bounds(lb,ub)

	#Objective function
	def _objfun_impl(self,x):
		#1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
		T = list([0]*(self.__n_legs))

		for i in xrange(self.__n_legs-1):	
			j = i+1;
			T[-j] = (x[5] - sum(T[-(j-1):])) * x[-1-(j-1)*4]
		T[0] = x[5] - sum(T)
		
		#2 - We compute the epochs and ephemerides of the planetary encounters
		t_P = list([None] * (self.__n_legs+1))
		r_P = list([None] * (self.__n_legs+1))
		v_P = list([None] * (self.__n_legs+1))
		DV = list([None] * (self.__n_legs+1))
		
		for i,planet in enumerate(self.seq):
			t_P[i] = epoch(x[0] + sum(T[0:i]))
			r_P[i],v_P[i] = self.seq[i].eph(t_P[i])

		#3 - We start with the first leg
		theta = 2*pi*x[1]
		phi = acos(2*x[2]-1)-pi/2

		Vinfx = x[3]*cos(phi)*cos(theta)
		Vinfy =	x[3]*cos(phi)*sin(theta)
		Vinfz = x[3]*sin(phi)

		v0 = [a+b for a,b in zip(v_P[0],[Vinfx,Vinfy,Vinfz])]
		r,v = propagate_lagrangian(r_P[0],v0,x[4]*T[0]*DAY2SEC,MU_SUN)

		#Lambert arc to reach seq[1]
		dt = (1-x[4])*T[0]*DAY2SEC
		l = lambert_problem(r,r_P[1],dt,MU_SUN)
		v_end_l = l.get_v2()[0]
		v_beg_l = l.get_v1()[0]

		#First DSM occuring at time nu1*T1
		DV[0] = norm([a-b for a,b in zip(v_beg_l,v)])

		#4 - And we proceed with each successive leg
		for i in range(1,self.__n_legs):
			#Fly-by 
			v_out = fb_prop(v_end_l,v_P[i],x[7+(i-1)*4]*self.seq[i].radius,x[6+(i-1)*4],self.seq[i].mu_self)
			#s/c propagation before the DSM
			r,v = propagate_lagrangian(r_P[i],v_out,x[8+(i-1)*4]*T[i]*DAY2SEC,MU_SUN)
			#Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
			dt = (1-x[8+(i-1)*4])*T[i]*DAY2SEC
			l = lambert_problem(r,r_P[i+1],dt,MU_SUN)
			v_end_l = l.get_v2()[0]
			v_beg_l = l.get_v1()[0]
			#DSM occuring at time nu2*T2
			DV[i] = norm([a-b for a,b in zip(v_beg_l,v)])

		#Last Delta-v
		DV[-1] = norm([a-b for a,b in zip(v_end_l,v_P[-1])])

		if self.f_dimension == 1:
			return (sum(DV),)
		else:
			return (sum(DV), sum(T))

	def pretty(self,x):
		"""
		Prints human readable information on the trajectory represented by the decision vector x
		
		Example::
		
		  prob.pretty(x)
		"""
		#1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
		T = list([0]*(self.__n_legs))
		#a[-i] = x[-1-(i-1)*4]
		for i in xrange(self.__n_legs-1):	
			j = i+1;
			T[-j] = (x[5] - sum(T[-(j-1):])) * x[-1-(j-1)*4]
		T[0] = x[5] - sum(T)
		
		#2 - We compute the epochs and ephemerides of the planetary encounters
		t_P = list([None] * (self.__n_legs+1))
		r_P = list([None] * (self.__n_legs+1))
		v_P = list([None] * (self.__n_legs+1))
		DV = list([None] * (self.__n_legs+1))
		
		for i,planet in enumerate(self.seq):
			t_P[i] = epoch(x[0] + sum(T[0:i]))
			r_P[i],v_P[i] = self.seq[i].eph(t_P[i])

		#3 - We start with the first leg
		print "First Leg: " + self.seq[0].name + " to " + self.seq[1].name 
		
		theta = 2*pi*x[1]
		phi = acos(2*x[2]-1)-pi/2

		Vinfx = x[3]*cos(phi)*cos(theta)
		Vinfy =	x[3]*cos(phi)*sin(theta)
		Vinfz = x[3]*sin(phi)
		
		print "Departure: " + str(t_P[0]) + " (" + str(t_P[0].mjd2000) + " mjd2000) " 
		print "Duration: " + str(T[0]) + "days"
		print "VINF: " + str(x[3] / 1000) + " km/sec"

		v0 = [a+b for a,b in zip(v_P[0],[Vinfx,Vinfy,Vinfz])]
		r,v = propagate_lagrangian(r_P[0],v0,x[4]*T[0]*DAY2SEC,MU_SUN)
		
		print "DSM after " + str(x[4]*T[0]) + " days"

		#Lambert arc to reach seq[1]
		dt = (1-x[4])*T[0]*DAY2SEC
		l = lambert_problem(r,r_P[1],dt,MU_SUN)
		v_end_l = l.get_v2()[0]
		v_beg_l = l.get_v1()[0]

		#First DSM occuring at time nu1*T1
		DV[0] = norm([a-b for a,b in zip(v_beg_l,v)])
		print "DSM magnitude: " + str(DV[0]) + "m/s"

		#4 - And we proceed with each successive leg
		for i in range(1,self.__n_legs):
			print "\nleg no. " + str(i+1) + ": " + self.seq[i].name + " to " + self.seq[i+1].name 
			print "Duration: " + str(T[i]) + "days"
			#Fly-by 
			v_out = fb_prop(v_end_l,v_P[i],x[7+(i-1)*4]*self.seq[i].radius,x[6+(i-1)*4],self.seq[i].mu_self)
			print "Fly-by epoch: " + str(t_P[i]) + " (" + str(t_P[i].mjd2000) + " mjd2000) " 
			print "Fly-by radius: " + str(x[7+(i-1)*4]) + " planetary radii"
			#s/c propagation before the DSM
			r,v = propagate_lagrangian(r_P[i],v_out,x[8+(i-1)*4]*T[i]*DAY2SEC,MU_SUN)
			print "DSM after " + str(x[8+(i-1)*4]*T[i]) + " days"
			#Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
			dt = (1-x[8+(i-1)*4])*T[i]*DAY2SEC
			l = lambert_problem(r,r_P[i+1],dt,MU_SUN)
			v_end_l = l.get_v2()[0]
			v_beg_l = l.get_v1()[0]
			#DSM occuring at time nu2*T2
			DV[i] = norm([a-b for a,b in zip(v_beg_l,v)])
			print "DSM magnitude: " + str(DV[i]) + "m/s"

		#Last Delta-v
		print "\nArrival at " + self.seq[-1].name
		DV[-1] = norm([a-b for a,b in zip(v_end_l,v_P[-1])])
		print "Arrival epoch: " + str(t_P[-1]) + " (" + str(t_P[-1].mjd2000) + " mjd2000) " 
		print "Arrival Vinf: " + str(DV[-1]) + "m/s"
		print "Total mission time: " + str(sum(T)/365.25) + " years"


	#Plot of the trajectory
	def plot(self,x):
		"""
		Plots the trajectory represented by the decision vector x
		
		Example::
		
		  prob.plot(x)
		"""
		import matplotlib as mpl
		from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
		from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler

		mpl.rcParams['legend.fontsize'] = 10
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		ax.scatter(0,0,0, color='y')
		
		#1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
		T = list([0]*(self.__n_legs))
		#a[-i] = x[-1-(i-1)*4]
		for i in xrange(self.__n_legs-1):	
			j = i+1;
			T[-j] = (x[5] - sum(T[-(j-1):])) * x[-1-(j-1)*4]
		T[0] = x[5] - sum(T)
		
		#2 - We compute the epochs and ephemerides of the planetary encounters
		t_P = list([None] * (self.__n_legs+1))
		r_P = list([None] * (self.__n_legs+1))
		v_P = list([None] * (self.__n_legs+1))
		DV = list([None] * (self.__n_legs+1))
		
		for i,planet in enumerate(self.seq):
			t_P[i] = epoch(x[0] + sum(T[0:i]))
			r_P[i],v_P[i] = planet.eph(t_P[i])
			plot_planet(ax, planet, t0=t_P[i], color=(0.8,0.6,0.8), legend=True, units = AU)

		#3 - We start with the first leg
		theta = 2*pi*x[1]
		phi = acos(2*x[2]-1)-pi/2

		Vinfx = x[3]*cos(phi)*cos(theta)
		Vinfy =	x[3]*cos(phi)*sin(theta)
		Vinfz = x[3]*sin(phi)

		v0 = [a+b for a,b in zip(v_P[0],[Vinfx,Vinfy,Vinfz])]
		r,v = propagate_lagrangian(r_P[0],v0,x[4]*T[0]*DAY2SEC,MU_SUN)
		plot_kepler(ax,r_P[0],v0,x[4]*T[0]*DAY2SEC,MU_SUN,N = 100, color='b', legend=False, units = AU)

		#Lambert arc to reach seq[1]
		dt = (1-x[4])*T[0]*DAY2SEC
		l = lambert_problem(r,r_P[1],dt,MU_SUN)
		plot_lambert(ax,l, sol = 0, color='r', legend=False, units = AU)
		v_end_l = l.get_v2()[0]
		v_beg_l = l.get_v1()[0]

		#First DSM occuring at time nu1*T1
		DV[0] = norm([a-b for a,b in zip(v_beg_l,v)])

		#4 - And we proceed with each successive leg
		for i in range(1,self.__n_legs):
			#Fly-by 
			v_out = fb_prop(v_end_l,v_P[i],x[7+(i-1)*4]*self.seq[i].radius,x[6+(i-1)*4],self.seq[i].mu_self)
			#s/c propagation before the DSM
			r,v = propagate_lagrangian(r_P[i],v_out,x[8+(i-1)*4]*T[i]*DAY2SEC,MU_SUN)
			plot_kepler(ax,r_P[i],v_out,x[8+(i-1)*4]*T[i]*DAY2SEC,MU_SUN,N = 100, color='b', legend=False, units = AU)
			#Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
			dt = (1-x[8+(i-1)*4])*T[i]*DAY2SEC
			l = lambert_problem(r,r_P[i+1],dt,MU_SUN)
			plot_lambert(ax,l, sol = 0, color='r', legend=False, units = AU)
			v_end_l = l.get_v2()[0]
			v_beg_l = l.get_v1()[0]
			#DSM occuring at time nu2*T2
			DV[i] = norm([a-b for a,b in zip(v_beg_l,v)])

		plt.show()
	
	def set_tof(self, minimum, maximum):
		"""
		Sets the minimum and maximum time of flight allowed (in years)
		
		Example::
		  m = 3
		  M = 5
		  prob.set_tof(m,M)
		"""
		lb = list(self.lb)
		ub = list(self.ub)
		lb[5] = minimum*365.25
		ub[5] = maximum*365.25
		self.set_bounds(lb,ub)
		
	def set_launch_window(self, start, end):
		"""
		Sets the launch window allowed in terms of starting and ending epochs
		
		Example::
		
		  start = epoch(0)
		  end = epoch(1000)
		  prob.set_launch_window(start, end)
		"""
		lb = list(self.lb)
		ub = list(self.ub)
		lb[0] = start.mjd2000
		ub[0] = end.mjd2000
		self.set_bounds(lb,ub)
		
	def set_vinf(self, vinf):
		"""
		Sets the allowed launch vinf (in km/s)
		
		Example::
		  
		  M = 5
		  prob.set_vinf(M)
		"""
		lb = list(self.lb)
		ub = list(self.ub)
		lb[3] = 0
		ub[3] = vinf * 1000
		self.set_bounds(lb,ub)