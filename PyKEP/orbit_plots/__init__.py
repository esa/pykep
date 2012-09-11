"""
This module contains functions that allow 3D plots of trajectories and orbits.
It makes use of matplotlib mplot 3D extension and thus is fairly in testing phase
"""
from PyKEP import __extensions__

if (__extensions__['mplot3d']):
	from _plots import *

