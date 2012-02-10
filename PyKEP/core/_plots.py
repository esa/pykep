def plot_planet(ax,plnt,t0='PyKEP.epoch(0)', N=60, units = 'PyKEP.AU', color = 'k', legend = False):
	"""
	Plots the planet position and its orbit
		      
	USAGE: plot_planet(ax,plnt,t0='PyKEP.epoch(0)', N=60, units=AU, legend = False):
	  * ax = 3D axis object created using fig.gca(projection='3d')
	  * plnt = PyKEP.planet object we want to plot
	  * N = number of points to be plotted along one orbit
	  * t0 = PyKEP.epoch object indicating when we want to plot the planet position
	  * units = in meters. The length unit to be used in the plot
	  * legend = when True it plots also the legend with the planet name
	  
	EXAMPLE:
	  from mpl_toolkits.mplot3d import Axes3D
	  import matplotlib.pyplot as plt
	
	  fig = plt.figure()
	  ax = fig.gca(projection='3d')
	  pl = planet_ss('earth')
	  plot_planet(ax,pl)
	  plt.show()
	"""
	from PyKEP import MU_SUN, SEC2DAY, epoch, AU
	from math import pi,sqrt
	import numpy as np
	
	if t0 == 'PyKEP.epoch(0)':
		t0 = epoch(0)
	if units == 'PyKEP.AU':
		units = AU
	
	#orbit semi-major axis
	a = plnt.orbital_elements[0]
	
	#orbital period in days
	T = 2*pi*sqrt(a**3/MU_SUN) * SEC2DAY
	
	#points where the orbit will be plotted
	when = np.linspace(0,T,N)
	
	#Ephemerides Calculation for the given planet
	x = np.array([0.0]*N)
	y = np.array([0.0]*N)
	z = np.array([0.0]*N)
	
	for i,day in enumerate(when):
		r,v = plnt.eph(epoch(t0.mjd2000 + day))
		x[i] = r[0]/units
		y[i] = r[1]/units
		z[i] = r[2]/units

	#Actual plot commands
	ax.plot(x, y, z, label=plnt.name + " " +  t0.__repr__()[0:11], c=color)
	ax.scatter(x[0],y[0],z[0])
	
	if legend:
		ax.legend()
	
def plot_lambert(ax,l, N=60, sol=0, units = 'PyKEP.AU', color = 'b', legend = False):
	"""
	Plots a particular solution to a Lambert's problem
	      
	USAGE: plot_lambert(ax,l, N=60, sol=0, units = 'PyKEP.AU', legend = 'False')
	  * ax = 3D axis object created using fig.gca(projection='3d')
	  * l = PyKEP.lambert_problem object
	  * N = number of points to be plotted along one arc
	  * sol = solution to the Lambert's problem we want to plot (must be in 0..Nmax*2)
	          where Nmax is the maximum number of revolutions for which there exist a solution.
	  * units = in meters. The length unit to be used in the plot
	  * legend = when True it plots also the legend with info on the sol plotted
	  
	EXAMPLE:
	  from mpl_toolkits.mplot3d import Axes3D
	  import matplotlib.pyplot as plt
	  
	  t1 = epoch(0)
	  t2 = epoch(540)
	  dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC
  
	  pl = planet_ss('earth')
	  plot_planet(ax,pl, t0=t1)
	  rE,vE = pl.eph(t1)

	  pl = planet_ss('mars')
	  plot_planet(ax,pl, t0=t2)
	  rM, vM = pl.eph(t2)

	  l = lambert_problem(rE,rM,dt,MU_SUN)
	  plot_lambert(ax,l, sol=0)
	  plot_lambert(ax,l, sol=1)

	  plt.show()
	"""
	from PyKEP import propagate_lagrangian, AU
	import numpy as np
	
	if units == 'PyKEP.AU':
		units = AU
		
	if sol > l.get_Nmax()*2:
		raise ValueError("sol must be in 0 .. NMax*2 \n * Nmax is the maximum number of revolutions for which there exist a solution to the Lambert's problem \n * You can compute Nmax calling the get_Nmax() method of the lambert_problem object")
		
	#We extract the relevant information from the Lambert's problem
	r = l.get_r1()
	v = l.get_v1()[sol]
	T = l.get_tof()
	mu = l.get_mu()
	
	#We define the integration time ...
	dt = T / (N-1)
	
	#... and alocate the cartesian components for r
	x = np.array([0.0]*N)
	y = np.array([0.0]*N)
	z = np.array([0.0]*N)
	
	#We calculate the spacecraft position at each dt
	for i in range(N):
		x[i] = r[0]/units
		y[i] = r[1]/units
		z[i] = r[2]/units
		r,v = propagate_lagrangian(r,v,dt,mu)
	
	#And we plot
	ax.plot(x, y, z, c=color, label='Lambert solution (' + str((sol+1)/2) + ' revs.)')
	
	if legend:
		ax.legend()
