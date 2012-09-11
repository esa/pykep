"""
This module contains models of interplanetary trajectories
"""
from PyKEP import __extensions__

if (__extensions__['pygmo']):
	from _mga_1dsm import mga_1dsm
	from _mga_lt_nep import mga_lt_nep

