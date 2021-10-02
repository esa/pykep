from math import pi, sqrt, cos, sin
import math

from pykep.core import AU, RAD2DEG, SEC2DAY
from pykep.core.core import propagate_lagrangian, ic2par, epoch, \
    propagate_taylor

import numpy as np

 
class rvt:
    
    """
    Keplerian orbit represented by radius, velocity, time and mu.    
    """
    
    def __init__(self, r, v, time, mu):
        """
        Args:
            - r (``tuple`` of ``float``): cartesian position in m.
            - v: (``tuple`` of ``float``): velocity in m/s.
            - time: (``float``): time in seconds.
            - mu (`float``): gravity parameter of the central body.
        """

        self._r = r
        self._v = v
        self._t = time  # in seconds
        self._mu = mu

    # useful for debugging        
    def __str__(self):
        a, e, i, _, _, _ = self.kepler()
        period = 2 * pi * sqrt(a ** 3 / self._mu)
        apo = a * (1 + e) / AU
        per = a * (1 - e) / AU
        return str(self._r) + " " + str(self._v) + " " + str(self._t * SEC2DAY) + " " + \
                str(apo) + " " + str(per) + " " + \
                str(e) + " " + str(i * RAD2DEG) + " " + str(period * SEC2DAY)
        
    def propagate_lagrangian(self, tof):
        orb = rvt(self._r, self._v, self._t + tof, self._mu)
        orb._r, orb._v = propagate_lagrangian(orb._r, orb._v, tof, self._mu)
        return orb
 
    def propagate_taylor(self, tof, m0, thrust, veff=1, log10tol=-15, log10rtol=-15):
        orb = rvt(self._r, self._v, self._t + tof, self._mu)
        orb._r, orb._v, m = propagate_taylor(orb._r, orb._v, m0, thrust, tof, self._mu,
                                             veff, log10tol, log10rtol)
        return orb, m
    
    # keplarian parameters a, e, i, W, w, E       
    def kepler(self):
        return ic2par(self._r, self._v, self._mu) 
    
    # plots orbit from current time up to time + tof    
    def plot(self, tof, N=60, units=AU, color="b", label=None, axes=None):
        from pykep.orbit_plots import plot_kepler
        plot_kepler(r0=self._r, v0=self._v, tof=tof,
                    mu=self._mu, N=N, units=units, color=color,
                    label=label, axes=axes)

    def period(self):
        kep = ic2par(self._r, self._v, self._mu) 
        a = kep[0]
        meanMotion = sqrt(self._mu / (a ** 3))
        return 2.0 * math.pi / meanMotion;  # in seconds
    
    def rotate(self, k, theta):
        orb = rvt(self._r, self._v, self._t, self._mu)
        orb._r = rotate_vector(self._r, k, theta)
        orb._v = rotate_vector(self._v, k, theta)
        return orb 
    
    def tof(self, rvt2):
        return rvt2._t - self._t


def rvt_planet(pl, time):
    r, v = pl.eph(epoch(time * SEC2DAY))
    return rvt(r, v, time, pl.mu_central_body)


def rotate_vector(v, k, theta):
    dP = np.dot(k, v)
    cosTheta = cos(theta)
    sinTheta = sin(theta)
    # rotate using Rodrigues rotation formula
    r_rot = [
        a * cosTheta + b * sinTheta + c * (1 - cosTheta) * dP
        for a, b, c in zip(v, np.cross(k, v), k)
    ]
    return r_rot

