###########################################################
# We check what extensions are available
###########################################################
__extensions__ = {'matplotlib': False, 'mplot3d': False, 'pygmo': False, 'scikit-learn': False, 'scipy': False}
# 1 - matplotlib
try:
    from matplotlib import __version__ as matplotlib_ver
    __extensions__['matplotlib'] = True

    # We detect the version and if more than 1.1.0 mplot3d is there
    mver = matplotlib_ver.split('.')
    mver = int(mver[0]) * 100 + int(mver[1]) * 10 + int(mver[2])
    if mver >= 110:
        __extensions__['mplot3d'] = True
    del mver
except ImportError:
    pass

# 2 - PyGMO
try:
    from PyGMO import __version__ as pygmo_ver
    __extensions__['pygmo'] = True
except ImportError:
    pass

# 3 - scikit-learn is installed
try:
    from sklearn import __version__ as sklearn_ver
    __extensions__['scikit-learn'] = True
except ImportError:
    pass

# 4 - scipy is installed
try:
    from scipy import __version__ as scipy_ver
    __extensions__['scipy'] = True
except ImportError:
    pass

###########################################################
# We import the submodules
###########################################################
from PyKEP import core, sims_flanagan, orbit_plots, examples, trajopt, phasing

# For convenience, bring all core classes into the root namespace when importing *.
from PyKEP.core import *

###########################################################
# We define PyKEP module
###########################################################
__doc__ = 'PyKEP is the answer ... but what was the question?'
__all__ = ['core', 'sims_flanagan', 'orbit_plots', 'examples', 'trajopt', 'phasing']
__version__ = {'major': 1, 'minor': 2, 'bugfix': 0}
__all__ += [name for name in dir(core) if not name.startswith('_')]
