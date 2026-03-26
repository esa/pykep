## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)## 
## This file is part of the kep3 library.## 
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

from sgp4.api import Satrec
from sgp4 import exporter
import pykep as pk

from sgp4.api import Satrec
from sgp4 import exporter
import pykep as pk
import numpy as np

class tle:
    """__init__(line1, line2)

    This User Defined Planet (UDPLA) represents a satellite orbiting the Earth and defined in the TLE format
    and propagated using the SGP4 propagator.

    Args:
        *line1* (:class:`str`): The first line of a TLE
        
        *line2* (:class:`str`): The second line of a TLE

    .. note::
       The resulting ephemerides will be returned in SI units and in the True Equator Mean Equinox (TEME) reference frame

    Examples:
      >>> import pykep as pk
      >>> line1 = "1 33773U 97051L   23290.57931959  .00002095  00000+0  65841-3 0  9991"
      >>> line2 = "2 33773  86.4068  33.1145 0009956 224.5064 135.5336 14.40043565770064"
      >>> udpla = pk.udpla.tle(line1, line2)
      >>> pla = pk.planet(udpla)
      >>> pla.eph(pk.epoch("2023-10-31"))
    """
    def __init__(self, line1, line2):
        import pykep as pk
        self.satellite = Satrec.twoline2rv(line1, line2)
        self.e = 0
        self.ref_epoch = pk.epoch(self.satellite.jdsatepoch + self.satellite.jdsatepochF, pk.epoch.julian_type.JD)

    def eph(self, mjd2000):
        """Mandatory method of the :class:`~pykep.planet` interface.

        Args:
            *mjd2000* (:class:`float`): Modified Julian Date 2000.

        Returns:
            :class:`list` [:class:`list`, :class:`list`]: the planet ephemerides.
        """
        jd = mjd2000 + 2451544.5
        jd_i = int(jd)
        jd_fr = jd-jd_i
        self.e, r, v = self.satellite.sgp4(jd_i, jd_fr)
        return [[it*1000 for it in r], [it*1000 for it in v]]
    
    def eph_v(self, mjd2000s):
        """Optional method of the :class:`~pykep.planet` interface.

        Args:
            *mjd2000s* (array-like (1,)): Modified Julian Dates 2000.

        Returns:
            :class:`ndarray`: the planet ephemerides as an array of dimension (-1,6)
        """
        jds = [mjd2000 + 2451544.5 for mjd2000 in mjd2000s]
        jd_is = [int(item) for item in jds]
        jd_frs = [a-b for a,b in zip(jds,jd_is)]
        self.e, r, v = self.satellite.sgp4_array(np.array(jd_is), np.array(jd_frs))
        rv = np.hstack((r,v))
        return rv.reshape((-1, 6)) * 1000
    
    def get_name(self):
        """Optional method of the :class:`~pykep.planet` interface.

        Returns:
            :class:`str`: The body name
        """
        return self.satellite.satnum_str + " - SGP4"
    
    def get_extra_info(self):
        """Optional method of the :class:`~pykep.planet` interface.

        Returns:
            :class:`str`: Extra info on the udpla
        """
        line1, line2 = exporter.export_tle(self.satellite)
        return "TLE line1: " + line1 + "\nTLE line2: " + line2 
    
    def get_mu_central_body(self):
        """Optional method of the :class:`~pykep.planet` interface.

        Returns:
            :class:`float`: the graviational parameter of the Earth.
        """
        return pk.MU_EARTH