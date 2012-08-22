def plot_planet(ax,plnt,t0='PyKEP.epoch(0)', N=60, units = 1.0, color = 'k', legend = False):
	"""
	Plots the planet position and its orbit
		      
	USAGE: plot_planet(ax,plnt,t0='PyKEP.epoch(0)', N=60, units=AU, legend = False):
	  * ax:		3D axis object created using fig.gca(projection='3d')
	  * plnt:	PyKEP.planet object we want to plot
	  * t0:		PyKEP.epoch object indicating when we want to plot the planet position
	  * units:	the length unit to be used in the plot
	  * color:	matplotlib color to use to plot the line
	  * legend	when True it plots also the legend with the planet name
	  
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
	if legend:
		label=plnt.name + " " +  t0.__repr__()[0:11]
	else:
		label=None
	ax.plot(x, y, z, label=label, c=color)
	ax.scatter([x[0]],[y[0]],[z[0]])
	
	if legend:
		ax.legend()
	
def plot_lambert(ax,l, N=60, sol=0, units = 1.0, color = 'b', legend = False):
	"""
	Plots a particular solution to a Lambert's problem
	      
	USAGE: plot_lambert(ax,l, N=60, sol=0, units = 'PyKEP.AU', legend = 'False')
	  * ax:		3D axis object created using fig.gca(projection='3d')
	  * l:		PyKEP.lambert_problem object
	  * N:		number of points to be plotted along one arc
	  * sol:	solution to the Lambert's problem we want to plot (must be in 0..Nmax*2)
			where Nmax is the maximum number of revolutions for which there exist a solution.
	  * units:	the length unit to be used in the plot
	  * color:	matplotlib color to use to plot the line
	  * legend	when True it plots also the legend with info on the Lambert's solution chosen
	  
	EXAMPLE:
	  from mpl_toolkits.mplot3d import Axes3D
	  import matplotlib.pyplot as plt
	
	  fig = plt.figure()
	  ax = fig.gca(projection='3d')

	  t1 = epoch(0)
	  t2 = epoch(640)
	  dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC

	  pl = planet_ss('earth')
	  plot_planet(ax,pl, t0=t1, color='k')
	  rE,vE = pl.eph(t1)

	  pl = planet_ss('mars')
	  plot_planet(ax,pl, t0=t2, color='r')
	  rM, vM = pl.eph(t2)

	  l = lambert_problem(rE,rM,dt,MU_SUN)
	  plot_lambert(ax,l, color='b')
	  plot_lambert(ax,l,sol=1, color='g')
	  plot_lambert(ax,l,sol=2, color='g')

	  plt.show()
	"""
	from PyKEP import propagate_lagrangian, AU
	import numpy as np
	
		
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
	if legend:
		label = 'Lambert solution (' + str((sol+1)/2) + ' revs.)'
	else:
		label = None
	ax.plot(x, y, z, c=color, label=label)
	
	if legend:
		ax.legend()

def plot_kepler(ax,r,v,t,mu, N=60, units = 1, color = 'b', legend = False):
	"""
	Plots the result of a keplerian propagation
		      
	USAGE: plot_kepler(ax,r,v,t,mu, N=60, units = 1, color = 'b', legend = False):
	  * ax:		3D axis object created using fig.gca(projection='3d')
	  * r:		initial position (cartesian coordinates)
	  * v:		initial velocity (cartesian coordinates)
	  * t:		propagation time
	  * mu:		gravitational parameter 
	  * N:		number of points to be plotted along one arc
	  * units:	the length unit to be used in the plot
	  * color:	matplotlib color to use to plot the line
	  * legend	when True it plots also the legend
	"""
	
	from PyKEP import propagate_lagrangian
  
	#We define the integration time ...
	dt = t / (N-1)
	
	#... and calcuate the cartesian components for r
	x = [0.0]*N
	y = [0.0]*N
	z = [0.0]*N
	
	#We calculate the spacecraft position at each dt
	for i in range(N):
		x[i] = r[0]/units
		y[i] = r[1]/units
		z[i] = r[2]/units
		r,v = propagate_lagrangian(r,v,dt,mu)
	
	#And we plot
	if legend:
		label = 'ballistic arc'
	else:
		label = None
	ax.plot(x, y, z, c=color, label=label)
	
	if legend:
		ax.legend()
		
def plot_taylor(ax,r,v,m,u,t,mu, veff, N=60, units = 1, color = 'b', legend = False):
	"""
	Plots the result of a taylor propagation of constant thrust
		      
	USAGE: plot_taylor(ax,r,v,m,u,t,mu, veff, N=60, units = 1, color = 'b', legend = False):
	  * ax:		3D axis object created using fig.gca(projection='3d')
	  * r:		initial position (cartesian coordinates)
	  * v:		initial velocity (cartesian coordinates)
	  * m: 		initial mass
	  * u:		cartesian components for the constant thrust
	  * t:		propagation time
	  * mu:		gravitational parameter 
	  * veff:	the product Isp * g0
	  * N:		number of points to be plotted along one arc
	  * units:	the length unit to be used in the plot
	  * color:	matplotlib color to use to plot the line
	  * legend:	when True it plots also the legend
	"""
	
	from PyKEP import propagate_taylor
  
	#We define the integration time ...
	dt = t / (N-1)
	
	#... and calcuate the cartesian components for r
	x = [0.0]*N
	y = [0.0]*N
	z = [0.0]*N
	
	#We calculate the spacecraft position at each dt
	for i in range(N):
		x[i] = r[0]/units
		y[i] = r[1]/units
		z[i] = r[2]/units
		r,v,m = propagate_taylor(r,v,m,u,dt,mu,veff,-10,-10)
	
	#And we plot
	if legend:
		label = 'constant thrust arc'
	else:
		label = None
	ax.plot(x, y, z, c=color, label=label)
	
	if legend:
		ax.legend()		

def plot_sf_leg(ax, leg, N=5, units=1, color='b', legend=False, plot_line = True):
	"""
	Plots a Sims-Flanagan leg
		      
	USAGE: plot_sf_leg(ax, leg, N=5, units=1, color='b', legend=False, no_trajectory = False):
	  * ax:		3D axis object created using fig.gca(projection='3d')
	  * leg:	a PyKEP.sims_flanagan.leg
	  * N:		number of points to be plotted along one arc
	  * units:	the length unit to be used in the plot
	  * color:	matplotlib color to use to plot the trajectory and the grid points
	  * legend	when True it plots also the legend
	  * plot_line: 	when True plots also the trajectory (between mid-points and grid points)
	  
	EXAMPLE:
	    from mpl_toolkits.mplot3d import Axes3D
	    import matplotlib.pyplot as plt

	    fig = plt.figure()
	    ax = fig.gca(projection='3d')
	    t1 = epoch(0)
	    pl = planet_ss('earth')
	    rE,vE = pl.eph(t1)
	    plot_planet(ax, pl,t0=t1, units=AU)

	    t2 = epoch(440)
	    pl = planet_ss('mars')
	    rM, vM = pl.eph(t2)
	    plot_planet(ax, pl,t0=t2, units=AU)

	    sc = sims_flanagan.spacecraft(4500,0.5,2500)
	    x0 = sims_flanagan.sc_state(rE,vE,sc.mass)
	    xe = sims_flanagan.sc_state(rM,vM,sc.mass)
	    l = sims_flanagan.leg(t1,x0,[1,0,0]*5,t2,xe,sc,MU_SUN)

	    plot_sf_leg(ax,l,units=AU)
	"""
	from PyKEP import propagate_lagrangian, AU, DAY2SEC, G0, propagate_taylor
	import numpy as np
	from scipy.linalg import norm
	from math import exp
	
	#We compute the number of segments for forward and backward propagation
	n_seg = len(leg.get_throttles())
	fwd_seg =  (n_seg+1)/2
	back_seg = n_seg/2

	#We extract information on the spacecraft
	sc = leg.get_spacecraft()
	isp = sc.isp
	max_thrust = sc.thrust
	
	#And on the leg
	throttles = leg.get_throttles()
	mu = leg.get_mu()
	
	#Forward propagation
	
	#x,y,z contain the cartesian components of all points (grid+midpints)
	x = [0.0]*(fwd_seg*2+1)
	y = [0.0]*(fwd_seg*2+1)
	z = [0.0]*(fwd_seg*2+1)
	
	state = leg.get_xi()
	
	#Initial conditions
	r = state.r; v = state.v; m = state.m
	x[0] = r[0]/units
	y[0] = r[1]/units
	z[0] = r[2]/units
	
	#We compute all points by propagation
	for i,t in enumerate(throttles[:fwd_seg]):
		dt = (t.end.mjd - t.start.mjd)*DAY2SEC
		alpha = min(norm(t.value),1.0)
		#Keplerian propagation and dV application
		if leg.high_fidelity == False:
			dV = [max_thrust / m * dt * dumb for dumb in t.value]
			if plot_line:
				plot_kepler(ax,r,v,dt/2,mu,N=N,units=units,color=(alpha,0,1-alpha))
			r,v = propagate_lagrangian(r,v,dt/2,mu)
			x[2*i+1] = r[0]/units
			y[2*i+1] = r[1]/units
			z[2*i+1] = r[2]/units
			#v= v+dV
			v = [a+b for a,b in zip(v,dV)]
			if plot_line:
				plot_kepler(ax,r,v,dt/2,mu,N=N,units=units,color=(alpha,0,1-alpha))
			r,v = propagate_lagrangian(r,v,dt/2,mu)
			x[2*i+2] = r[0]/units
			y[2*i+2] = r[1]/units
			z[2*i+2] = r[2]/units
			m *= exp( -norm(dV)/isp/G0 )
		#Taylor propagation of constant thrust u
		else:
			u = [max_thrust * dumb for dumb in t.value]
			if plot_line:
				plot_taylor(ax,r,v,m,u,dt/2,mu,isp*G0,N=N,units=units,color=(alpha,0,1-alpha))
			r,v,m = propagate_taylor(r,v,m,u,dt/2,mu,isp*G0,-10,-10)
			x[2*i+1] = r[0]/units
			y[2*i+1] = r[1]/units
			z[2*i+1] = r[2]/units
			if plot_line:
				plot_taylor(ax,r,v,m,u,dt/2,mu,isp*G0,N=N,units=units,color=(alpha,0,1-alpha))
			r,v,m = propagate_taylor(r,v,m,u,dt/2,mu,isp*G0,-10,-10)
			x[2*i+2] = r[0]/units
			y[2*i+2] = r[1]/units
			z[2*i+2] = r[2]/units
		
	x_grid = x[::2]; y_grid = y[::2]; z_grid = z[::2]
	x_midpoint = x[1::2]; y_midpoint = y[1::2]; z_midpoint = z[1::2]
	ax.scatter(x_grid[:-1], y_grid[:-1], z_grid[:-1], label='nodes',marker='o')
	ax.scatter(x_midpoint, y_midpoint, z_midpoint, label='mid-points',marker='x')
	ax.scatter(x_grid[-1],y_grid[-1],z_grid[-1], marker='^', c='y', label='mismatch point')
	
	#Backward propagation
	
	#x,y,z will contain the cartesian components of 
	x = [0.0]*(back_seg*2+1)
	y = [0.0]*(back_seg*2+1)
	z = [0.0]*(back_seg*2+1)
	
	state = leg.get_xf()
	
	#Final conditions
	r = state.r; v = state.v; m = state.m
	x[-1] = r[0]/units
	y[-1] = r[1]/units
	z[-1] = r[2]/units
	
	for i,t in enumerate(throttles[-1:-back_seg-1:-1]):
		dt = (t.end.mjd - t.start.mjd)*DAY2SEC
		alpha = min(norm(t.value),1.0)
		if leg.high_fidelity == False:
			dV = [max_thrust / m * dt * dumb for dumb in t.value]
			if plot_line:
				plot_kepler(ax,r,v,-dt/2,mu,N=N,units=units,color=(alpha,0,1-alpha))
			r,v = propagate_lagrangian(r,v,-dt/2,mu)
			x[-2*i-2] = r[0]/units
			y[-2*i-2] = r[1]/units
			z[-2*i-2] = r[2]/units
			#v= v+dV
			v = [a-b for a,b in zip(v,dV)]
			if plot_line:
				plot_kepler(ax,r,v,-dt/2,mu,N=N,units=units,color=(alpha,0,1-alpha))
			r,v = propagate_lagrangian(r,v,-dt/2,mu)
			x[-2*i-3] = r[0]/units
			y[-2*i-3] = r[1]/units
			z[-2*i-3] = r[2]/units
			m *= exp( norm(dV)/isp/G0 )
		else:
			u = [max_thrust * dumb for dumb in t.value]
			if plot_line:
				plot_taylor(ax,r,v,m,u,-dt/2,mu,isp*G0,N=N,units=units,color=(alpha,0,1-alpha))
			r,v,m = propagate_taylor(r,v,m,u,-dt/2,mu,isp*G0,-10,-10)
			x[-2*i-2] = r[0]/units
			y[-2*i-2] = r[1]/units
			z[-2*i-2] = r[2]/units
			if plot_line:
				plot_taylor(ax,r,v,m,u,-dt/2,mu,isp*G0,N=N,units=units,color=(alpha,0,1-alpha))
			r,v,m = propagate_taylor(r,v,m,u,-dt/2,mu,isp*G0,-10,-10)
			x[-2*i-3] = r[0]/units
			y[-2*i-3] = r[1]/units
			z[-2*i-3] = r[2]/units

	x_grid = x[::2]; y_grid = y[::2]; z_grid = z[::2]
	x_midpoint = x[1::2]; y_midpoint = y[1::2]; z_midpoint = z[1::2]
	
	ax.scatter(x_grid[1:], y_grid[1:], z_grid[1:], marker = 'o')
	ax.scatter(x_midpoint, y_midpoint, z_midpoint, marker = 'x')
	ax.scatter(x_grid[0],y_grid[0],z_grid[0], marker='^', c='y')
	
	if legend:
		ax.legend()
  
def plot_leg(ax, leg, N=5, units=1, color='b', legend=False, plot_line = True):
	"""
	Plots a trajectory leg
		      
	USAGE: plot_leg(ax, leg, N=5, units=1, color='b', legend=False, no_trajectory = False):
	  * ax:		3D axis object created using fig.gca(projection='3d')
	  * leg:	a PyKEP.sims_flanagan.leg
	  * N:		number of points to be plotted along one arc
	  * units:	the length unit to be used in the plot
	  * color:	matplotlib color to use to plot the trajectory and the grid points
	  * legend	when True it plots also the legend
	  * plot_line: 	when True plots also the trajectory (between mid-points and grid points)
	  
	EXAMPLE:
	    from mpl_toolkits.mplot3d import Axes3D
	    import matplotlib.pyplot as plt

	    fig = plt.figure()
	    ax = fig.gca(projection='3d')
	    t1 = epoch(0)
	    pl = planet_ss('earth')
	    rE,vE = pl.eph(t1)
	    plot_planet(ax, pl,t0=t1, units=AU)

	    t2 = epoch(440)
	    pl = planet_ss('mars')
	    rM, vM = pl.eph(t2)
	    plot_planet(ax, pl,t0=t2, units=AU)

	    sc = sims_flanagan.spacecraft(4500,0.5,2500)
	    x0 = sims_flanagan.sc_state(rE,vE,sc.mass)
	    xe = sims_flanagan.sc_state(rM,vM,sc.mass)
	    l = sims_flanagan.leg(t1,x0,[1,0,0]*5,t2,xe,sc,MU_SUN)

	    plot_sf_leg(ax,l,units=AU)
	"""
	from PyKEP import propagate_lagrangian, AU, DAY2SEC, G0, propagate_taylor
	import numpy as np
	from scipy.linalg import norm
	from math import exp
	
	#We compute the number of segments for forward and backward propagation
	n_seg = len(leg.get_throttles())
	fwd_seg =  (n_seg+1)/2
	back_seg = n_seg/2

	#We extract information on the spacecraft

	sc = leg.get_spacecraft()
	isp = sc.isp
	max_thrust = sc.thrust
	
	#And on the leg
	throttles = leg.get_throttles()
	mu = leg.get_mu()
	
	#Forward propagation
	
	#x,y,z contain the cartesian components of all points containd in states
	x = [0.0]*(fwd_seg+1)
	y = [0.0]*(fwd_seg+1)
	z = [0.0]*(fwd_seg+1)
	
	state = leg.get_xi()
	
	#Initial conditions
	r = state.r; v = state.v; m = state.m
	x[0] = r[0]/units
	y[0] = r[1]/units
	z[0] = r[2]/units
	
	#We compute all points by propagation
	for i,t in enumerate(throttles[:fwd_seg]):
		dt = (t.end.mjd - t.start.mjd) * DAY2SEC
		alpha = min(norm(t.value),1.0)
		
		#Taylor propagation of constant thrust u
		u = [max_thrust * dumb for dumb in t.value]
		if plot_line:
			plot_taylor(ax,r,v,m,u,dt,mu,isp*G0,N=N,units=units,color=(alpha,0,1-alpha))
		r,v,m = propagate_taylor(r,v,m,u,dt,mu,isp*G0,-10,-10)
		x[i+1] = r[0]/units
		y[i+1] = r[1]/units
		z[i+1] = r[2]/units
		
	ax.scatter(x[:-1], y[:-1], z[:-1], label='nodes',marker='o')
	ax.scatter(x[-1],y[-1],z[-1], marker='^', c='y', label='mismatch point')
	
	#Backward propagation
	
	#x,y,z will contain the cartesian components of 
	x = [0.0]*(back_seg+1)
	y = [0.0]*(back_seg+1)
	z = [0.0]*(back_seg+1)
	
	state = leg.get_xf()
	
	#Final conditions
	r = state.r; v = state.v; m = state.m
	x[-1] = r[0]/units
	y[-1] = r[1]/units
	z[-1] = r[2]/units
	
	for i,t in enumerate(throttles[-1:-back_seg-1:-1]):
		dt = (t.end.mjd - t.start.mjd)*DAY2SEC
		alpha = min(norm(t.value),1.0)

		u = [max_thrust * dumb for dumb in t.value]
		if plot_line:
			plot_taylor(ax,r,v,m,u,-dt,mu,isp*G0,N=N,units=units,color=(alpha,0,1-alpha))
		r,v,m = propagate_taylor(r,v,m,u,-dt,mu,isp*G0,-10,-10)
		x[-i-2] = r[0]/units
		y[-i-2] = r[1]/units
		z[-i-2] = r[2]/units


	ax.scatter(x[1:], y[1:], z[1:], marker = 'o')
	ax.scatter(x[0],  y[0],  z[0], marker='^', c='y')
	
	if legend:
		ax.legend()
  
