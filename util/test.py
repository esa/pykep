import numpy as np

from gravity_spherical_harmonic import gravity_spherical_harmonic
from load_gravity_model import load_gravity_model


path = "gravity_models/Earth/egm96.txt"
radius, mu, c_coef, s_coef, max_degree, max_order = load_gravity_model(path)

r_orbit = radius + 800
lat = np.pi/6
lon = np.pi/4

x1 = np.array([r_orbit * np.cos(lat) * np.cos(lon),
               r_orbit * np.cos(lat) * np.sin(lon),
               r_orbit * np.sin(lat)])

lat = -np.pi/2
lon = 3*np.pi/4

x2 = np.array([r_orbit * np.cos(lat) * np.cos(lon),
	       r_orbit * np.cos(lat) * np.sin(lon),
	       r_orbit * np.sin(lat)])

x = np.row_stack((x1, x2))

degree = 150
order = 150

acc = gravity_spherical_harmonic(x, radius, mu, c_coef, s_coef, degree, order)

matlab_x = np.array([[1.554201241795927e+03, 1.554201241795926e+03, 1.269e+03],
		     [-1.098898235362841e-13, 1.098898235362841e-13, -2538]])
matlab_acc = np.array([[-5.995780870165296e-03, -5.996501671863385e-03, -4.911882793026260e-03],
		       [9.331106606031310e-08, 4.034312166685789e-08, 9.763696823761261e-03]])

for i in range(len(x)):
    print(f"Acceleration error for position {x[i]}:\n{acc[i] - matlab_acc[i]}")
