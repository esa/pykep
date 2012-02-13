import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
from PyKEP import *


mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')

t1 = epoch(0)
t2 = epoch(640)
dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC

ax.scatter(0,0,0, color='y')

pl = planet_ss('earth')
plot_planet(ax,pl, t0=t1, color='k', legend=True, units = AU)
rE,vE = pl.eph(t1)

pl = planet_ss('mars')
plot_planet(ax,pl, t0=t2, color='r', legend=True, units = AU)
rM, vM = pl.eph(t2)

l = lambert_problem(rE,rM,dt,MU_SUN)
plot_lambert(ax,l, color='b', legend=True, units = AU)
plot_lambert(ax,l,sol=1, color='g', legend=True, units = AU)
plot_lambert(ax,l,sol=2, color='g', legend=True, units = AU)

plt.show()

