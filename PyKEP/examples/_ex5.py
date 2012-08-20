try:
	from PyGMO.problem import base as PyGMO_problem
	"""
	This constructs a PyGMO.problem object that represents an impulsive transfer Earth-Venus-Earth. The decision 
	vector is [t0,u,v,Vinf,n1,T1] + [beta, rp/rV, n2,T2] in the units: [mjd2000,nd,nd,m/s,nd,days] + [rad,nd,nd,days]
	where Vinf = Vinf_mag*(cos(theta)*cos(phi)i+cos(theta)*sin(phi)j+sin(phi)k) and theta = 2*pi*u and phi = acos(2*v-1)-pi/2
	"""
	class mga_dsm_EVE(PyGMO_problem):

		def __init__(self):
			#First we call the constructor for the base PyGMO problem 
			#As our problem is ten dimensional, box-bounded, single-objective
			#(dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
			super(mga_dsm_EVE,self).__init__(10)
	
			#We then define some data members (we use the double underscore to indicate they are private)
			from PyKEP import planet_ss, MU_SUN
			from math import pi
	   		self.__earth = planet_ss('earth')
			self.__venus = planet_ss('venus')
			self.set_bounds([5844,0,0,0,0.001,60,0,1.1,0.001,60],[6209,1.0,1.0,2.5*1000,0.999,500,2*pi,30,0.999,500])

		#Objective function
		def _objfun_impl(self,x):
			from PyKEP import epoch,DAY2SEC,planet_ss,MU_SUN,lambert_problem,propagate_lagrangian,fb_prop
			from math import cos, sin, pi, acos
			from scipy.linalg import norm

			t_E1 = epoch(x[0])
			t_V = epoch(x[0]+x[5])
			t_E2 = epoch(x[0]+x[5]+x[9])
			rE1, vE1 = self.__earth.eph(t_E1)
			rV, vV = self.__venus.eph(t_V)
			rE2, vE2 = self.__earth.eph(t_E2)

			#FIRST LEG: Earth-Venus
			#s/c propagation during nu1*T1 (first segment)
			theta = 2*pi*x[1]
			phi = acos(2*x[2]-1)-pi/2

			Vinfx = x[3]*cos(phi)*cos(theta)
			Vinfy =	x[3]*cos(phi)*sin(theta)
			Vinfz = x[3]*sin(phi)

			v0 = [a+b for a,b in zip(vE1,[Vinfx,Vinfy,Vinfz])]
			r0 = rE1
			r,v = propagate_lagrangian(r0,v0,x[4]*x[5]*DAY2SEC,MU_SUN)

			#Lambert arc to reach Venus during (1-nu1)*T1 (second segment)
			dt = (1-x[4])*x[5]*DAY2SEC
			l1 = lambert_problem(r,rV,dt,MU_SUN)
			v_end_l1 = l1.get_v2()[0]
			v_beg_l1 = l1.get_v1()[0]

			#First Delta-v occuring at time nu1*T1
			DV1 = [a-b for a,b in zip(v_beg_l1,v)]

			#SECOND LEG: Venus-Earth
			#Fly-by at Venus
			v_out = fb_prop(v_end_l1,vV,x[7]*self.__venus.radius,x[6],self.__venus.mu_self)
			r_out = rV

			#s/c propagation during nu2*T2 (first segment)
			r,v = propagate_lagrangian(r_out,v_out,x[8]*x[9]*DAY2SEC,MU_SUN)

			#Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
			dt = (1-x[8])*x[9]*DAY2SEC
			l2 = lambert_problem(r,rE2,dt,MU_SUN)
			v_end_l2 = l2.get_v2()[0]
			v_beg_l2 = l2.get_v1()[0]

			#Second Delta-v occuring at time nu2*T2
			DV2 = [a-b for a,b in zip(v_beg_l2,v)]

			#Third Delta-v occuring at time (1-nu2)*T2
			DV3 = [a-b for a,b in zip(v_end_l2,vE2)]

			#Total Delta-v
			DV_total = norm(DV1)+norm(DV2)+norm(DV3)
			
			return (DV_total,)

	
		#Plot of the trajectory
		def plot(self,x):
			import matplotlib as mpl
			from mpl_toolkits.mplot3d import Axes3D
			import matplotlib.pyplot as plt
			from math import cos, sin, pi, acos
			from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler
			from PyKEP import epoch, DAY2SEC, planet_ss, AU, MU_SUN, lambert_problem, propagate_lagrangian, fb_prop

			mpl.rcParams['legend.fontsize'] = 10
			fig = plt.figure()
			ax = fig.gca(projection='3d')
			ax.scatter(0,0,0, color='y')
			t_E1 = epoch(x[0])
			t_V = epoch(x[0]+x[5])
			t_E2 = epoch(x[0]+x[5]+x[9])
			rE1, vE1 = self.__earth.eph(t_E1)
			rV, vV = self.__venus.eph(t_V)
			rE2, vE2 = self.__earth.eph(t_E2)

			#FIRST LEG: Earth-Venus
			#s/c propagation during nu1*T1 (first segment)
			theta = 2*pi*x[1]
			phi = acos(2*x[2]-1)-pi/2

			Vinfx = x[3]*cos(phi)*cos(theta)
			Vinfy =	x[3]*cos(phi)*sin(theta)
			Vinfz = x[3]*sin(phi)

			v0 = [a+b for a,b in zip(vE1,[Vinfx,Vinfy,Vinfz])]
			r0 = rE1
			r,v = propagate_lagrangian(r0,v0,x[4]*x[5]*DAY2SEC,MU_SUN)
			plot_kepler(ax,r0,v0,x[4]*x[5]*DAY2SEC,MU_SUN,N = 100, color='b', legend=True, units = AU)

			#Lambert arc to reach Venus during (1-nu1)*T1 (second segment)
			dt = (1-x[4])*x[5]*DAY2SEC
			l1 = lambert_problem(r,rV,dt,MU_SUN)
			v_end_l1 = l1.get_v2()[0]
			v_beg_l1 = l1.get_v1()[0]
			plot_lambert(ax,l1, sol = 0, color='r', legend=True, units = AU)

			#First Delta-v occuring at time nu1*T1
			DV1 = [a-b for a,b in zip(v_beg_l1,v)]

			#SECOND LEG: Venus-Earth
			#Fly-by at Venus
			v_out = fb_prop(v_end_l1,vV,x[7]*self.__venus.radius,x[6],self.__venus.mu_self)
			r_out = rV

			#s/c propagation during nu2*T2 (first segment)
			r,v = propagate_lagrangian(r_out,v_out,x[8]*x[9]*DAY2SEC,MU_SUN)
			plot_kepler(ax,r_out,v_out,x[8]*x[9]*DAY2SEC,MU_SUN,N = 100, color='b', legend=True, units = AU)

			#Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
			dt = (1-x[8])*x[9]*DAY2SEC
			l2 = lambert_problem(r,rE2,dt,MU_SUN)
			v_end_l2 = l2.get_v2()[0]
			v_beg_l2 = l2.get_v1()[0]
			plot_lambert(ax,l2, sol = 0, color='r', legend=True, units = AU)				

			plot_planet(ax, self.__earth, t0=t_E1, color=(0.8,0.6,0.8), legend=True, units = AU)
			plot_planet(ax, self.__venus, t0=t_V, color=(0.8,0.8,0.2), legend=True, units = AU)
			plot_planet(ax, self.__earth, t0=t_E2, color=(0.8,0.6,0.8), legend=True, units = AU)
			plt.show()

	def run_example5():
		from PyGMO import algorithm, archipelago, problem
		prob = mga_dsm_EVE()
		algo = algorithm.jde(100)
		archi = archipelago(algo,prob,8,20)
		print "Running a Self-Adaptive Differential Evolution Algorithm .... on 8 parallel islands"
		archi.evolve(10); archi.join()
		isl = min(archi, key=lambda x:x.population.champion.f[0])
		prob.plot(isl.population.champion.x)

except:
	print 'Could not import PyGMO. PyGMO is required for some PyKEP examples'
	


