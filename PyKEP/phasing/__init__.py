"""
This module contains functions that allow to study the planetary phasing.
That is the relative planetaery position
"""
from PyKEP import __extensions__

if (__extensions__['scipy']):
    from _knn import *

if (__extensions__['scikit-learn']):
    from _dbscan import *
