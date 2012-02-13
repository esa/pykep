def run_example2():
	import matplotlib as mpl
	from mpl_toolkits.mplot3d import Axes3D

	import matplotlib.pyplot as plt
	from PyKEP import epoch, DAY2SEC, planet_ss, AU, MU_SUN, lambert_problem
	from PyKEP.orbit_plots import plot_planet, plot_lambert


	mpl.rcParams['legend.fontsize'] = 10

	fig = plt.figure()
	ax = fig.gca(projection='3d')

	t1 = epoch(0)
	t2 = epoch(640)
	dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC

	ax.scatter(0,0,0, color='y')

	pl = planet_ss('earth')
	plot_planet(ax,pl, t0=t1, color=(0.8,0.8,1), legend=True, units = AU)
	rE,vE = pl.eph(t1)

	pl = planet_ss('mars')
	plot_planet(ax,pl, t0=t2, color=(0.8,0.8,1), legend=True, units = AU)
	rM, vM = pl.eph(t2)

	l = lambert_problem(rE,rM,dt,MU_SUN)
	plot_lambert(ax,l, color='b', legend=True, units = AU)
	plot_lambert(ax,l,sol=1, color='g', legend=True, units = AU)
	plot_lambert(ax,l,sol=2, color='g', legend=True, units = AU)

	plt.show()

