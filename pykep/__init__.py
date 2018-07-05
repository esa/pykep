###########################################################
# We check what extensions are available
###########################################################
__extensions__ = {'matplotlib': False, 'mplot3d': False,
                  'pygmo': False, 'scikit-learn': False, 'scipy': False}
# 1 - matplotlib
try:
    from matplotlib import __version__ as matplotlib_ver
    __extensions__['matplotlib'] = True

    # We detect the version and if more than 1.1.0 mplot3d is there
    mver = matplotlib_ver.split('.')
    mver = int(mver[0]) * 100 + int(mver[1]) * 10
    if mver >= 110:
        __extensions__['mplot3d'] = True
    del mver
except ImportError:
    pass

# 2 - pygmo2
try:
    from pygmo import __version__ as pygmo_ver
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

# For convenience, bring all core classes into the root namespace when
# importing *.
from pykep.core import *
from pykep import core, sims_flanagan, pontryagin, orbit_plots, examples, phasing, util, planet, trajopt


###########################################################
# We define pykep module
###########################################################
version = '2.2'
__doc__ = 'pykep is the answer ... but what was the question?'
__all__ = ['core', 'sims_flanagan', 'pontryagin', 'orbit_plots',
           'examples', 'trajopt', 'phasing', 'util', 'planet']
__version__ = {'major': int(version.split('.')[0]), 'minor': int(version.split('.')[1])}
__all__ += [name for name in dir(core) if not name.startswith('_')]
