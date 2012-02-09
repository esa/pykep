# -*- coding: iso-8859-1 -*-
from _core import *

"""Extracts astronomical constants defined in the keplerian_toolbox file astro_constants.h"""
AU = _core._get_AU()
DAY2SEC = _core._get_DAY2SEC()
DAY2YEAR = _core._get_DAY2YEAR()
DEG2RAD = _core._get_DEG2RAD()
EARTH_VELOCITY = _core._get_EARTH_VELOCITY()
G0 = _core._get_G0()
MU_SUN = _core._get_MU_SUN()
RAD2DEG = _core._get_RAD2DEG()
SEC2DAY = _core._get_SEC2DAY()

epoch.epoch_type = _core._epoch_type

"""Various 3D plots"""

def plot_planet(ax,plnt,N=60, t0=epoch(0), units = AU):
  	"""
	Plots the planet position and its orbit
	
	NOTE: This function makes use of the axis3D object of matplotlib
	      
	USAGE: plot_planet(ax,plnt,N=60, t0=epoch(0), units=AU)
	  * ax = 3D axis object created using fig.gca(projection='3d')
	  * plnt = PyKEP.planet object we want to plot
	  * N = number of points to be plotted along one orbit
	  * t0 = PyKEP.epoch object indicating when we want to plot the planet position
	  * units = in meters, the unit to be used
	  
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
	ax.plot(x, y, z, label=plnt.name)
	ax.scatter(x[0],y[0],z[0])
	ax.legend()