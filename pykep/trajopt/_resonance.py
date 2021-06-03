import math

from pykep.core import epoch, fb_prop, SEC2DAY
from pykep.trajopt import _resonance
from pykep.trajopt._rvt import rvt

      
class resonance: 
    """
    Determines the best "fitting" resonance orbit.
    """

    def __init__(self, planet, rvt_in, rvt_pl,
                 resonances=[[1, 1], [5, 4], [4, 3], [3, 2], [5, 3]]
        ):
        """
        Args:
            - planet (``planet``): resonance planet. 
            - rvt_in: (``rvt``): incoming orbit.
            - rvt_pl: (``rvt``): planet orbit.
            - resonances (``list`` of ``int``): resonance options. 
        """
        assert rvt_in._t == rvt_pl._t  # timing must be consistent
        
        self._planet = planet
        self._rvt_in = rvt_in
        self._rvt_pl = rvt_pl       
        self._time = rvt_in._t
        self._resonances = resonances
        self._period = planet.compute_period(epoch(self._time * SEC2DAY))
        self._mu = planet.mu_self
        self._timing_error = -1
        self._rvt_out = None
        self._resonance = None
 
    # useful for debugging
    def __str__(self):
        return str(_resonance) + " " + str(self._timing_error * SEC2DAY) + " " + str(self._rvt_out)
    
    # select a resonance option minimizing the timing error  
    def select_resonance(self, beta, safe_distance):
        v_out = fb_prop(self._rvt_in._v, self._rvt_pl._v,
                        self._planet.radius + safe_distance, beta, self._mu)
        self._rvt_out = rvt(self._rvt_in._r, v_out, self._time, self._rvt_in._mu)
        period = self._rvt_out.period()
        self._timing_error = math.inf
        for resonance in self._resonances:
            target = self._period * resonance[1] / resonance[0];
            dt = abs(period - target)
            if dt < self._timing_error:
                self._resonance = resonance
                self._timing_error = dt
        return self._timing_error, self._resonance
    
    # time of flight of the resonance transfer
    def tof(self):
        return self._resonance[1] * self._period
