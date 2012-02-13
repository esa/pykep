__doc__ = 'PyKEP is the answer'
__all__ = ['core', 'sims_flanagan', 'orbit_plots', 'kep_examples']

"""Detecting Installed Extensions"""
# Fill up the __extensions__ variable with all detected extensions

__extensions__ = {'matplotlib': False, 'mplot3d': False,'pygmo': False}
try:
	import matplotlib
	__extensions__['matplotlib']=True
	del matplotlib
except:
	pass
try:
	import mpl_toolkits.mplot3d
	__extensions__['mplot3d']=True
	del mpl_toolkits
except:
	pass
try:
	import PyGMO
	__extensions__['pygmo']=True
	del PyGMO
except:
	pass
      
__version__ = '1.1.0'

# For convenience, bring all core classes into the root namespace when importing *.
from core import *
import sims_flanagan
import kep_examples

#Importing dependent modules
if (__extensions__['matplotlib'] == True):
	import orbit_plots

__all__ += filter(lambda name: not name.startswith('_'),dir(core))
