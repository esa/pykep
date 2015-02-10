"""
This module contains functions that allow to study the planetary phasing.
That is the relative planetary position
"""
from PyKEP import __extensions__

if (__extensions__['scipy']):
    from ._knn import *

if (__extensions__['scikit-learn']):
    from ._dbscan import *


def three_impulses_approximation(pl1, pl2, ep1=None, ep2=None):
    """
    Returns the DV in m/s necessary for an orbit transfer between pl1 and pl2 assuming a perfect phasing.
    The transfer will be made of three impulses. One to match apogees, one to match inclination and RAAN and
    one to match perigees. Two of the three impulses will be merged together at either departure or arrival.
    The argument of perigee is not matched, so that this approximation is only good for near-circular orbits.

    USAGE: DV = three_impulses_approximation(pl1, pl2)
           DV = three_impulses_approximation(pl1, pl2, ep1 = epoch(5500))
           DV = three_impulses_approximation(pl1, pl2, ep1 = epoch(5500), ep2 = epoch(5700))

    * pl1: departure planet
    * pl2: arrival planet
    * ep1: departure epoch (optional and only useful for non keplerian planets)
    * ep1: arrival epoch (default value is ep1).

    * DV: estimated DV cost for the orbital trab=nsfer
    """

    from PyKEP.core._core import _three_impulses_approximation
    if ep2 is None:
        ep2 = ep1
    if ep1 is None:
        return _three_impulses_approximation(pl1, pl2)
    return _three_impulses_approximation(pl1, pl2, ep1, ep2)
