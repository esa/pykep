# -*- coding: iso-8859-1 -*-
from PyKEP.core._core import *

"""Defining astronomical constants defined in the keplerian_toolbox file astro_constants.h"""
AU = _core._get_AU()
JR = _core._get_JR()
DAY2SEC = _core._get_DAY2SEC()
DAY2YEAR = _core._get_DAY2YEAR()
DEG2RAD = _core._get_DEG2RAD()
EARTH_VELOCITY = _core._get_EARTH_VELOCITY()
G0 = _core._get_G0()
MU_SUN = _core._get_MU_SUN()
RAD2DEG = _core._get_RAD2DEG()
SEC2DAY = _core._get_SEC2DAY()

from PyKEP.core._core import _epoch_type

EPOCH_TYPE = {
    "jd": _epoch_type.JD,
    "mjd": _epoch_type.MJD,
    "mjd2000": _epoch_type.MJD2000
}


def _epoch_ctor(self, julian_date, julian_date_type=None):
    """
PyKEP.epoch(julian_date, julian_date_type="mjd2000")

- julian_date: a julian date
- julian_date_type: julian date type, one of "jd", "mjd" and "mjd2000", defaults to "mjd2000"

Examples::

  e1 = epoch(0)
  e1 = epoch(0,"mjd2000")
  e2 = epoch(54333, "mjd")
    """
    if julian_date_type is None:
        self._orig_init(julian_date)
    else:
        self._orig_init(julian_date, EPOCH_TYPE[julian_date_type])
epoch._orig_init = epoch.__init__
epoch.__init__ = _epoch_ctor


def _planet_ctor(self, *args):
    """
PyKEP.planet(when,orbital_elements, mu_central_body, mu_self,radius, safe_radius [, name = 'unknown'])

PyKEP.planet(when,r,v, mu_central_body, mu_self,radius, safe_radius [, name = 'unknown'])

- when: a :py:class:`PyKEP.epoch` indicating the orbital elements epoch
- orbital_elements: a sequence of six containing a,e,i,W,w,M (SI units, i.e. meters and radiants)
- r,v: position and velocity of an object at when (SI units)
- mu_central_body: gravity parameter of the central body (SI units, i.e. m^2/s^3)
- mu_self: gravity parameter of the planet (SI units, i.e. m^2/s^3)
- radius: body radius (SI units, i.e. meters)
- safe_radius: mimimual radius that is safe during a fly-by of the planet (SI units, i.e. m)
- name: body name

.. note::

   use the derived classes :py:class:`PyKEP.planet_ss`, :py:class:`PyKEP.planet_mpcorb`, :py:class:`PyKEP.planet_js`, :py:class:`PyKEP.planet_tle` to instantiate common planets, asteroids or satellites

Example::

  earth = planet(epoch(54000,"mjd"),(9.99e-01 * AU, 1.67e-02, 8.85e-04 * DEG2RAD, 1.75e+02 * DEG2RAD, 2.87e+02 * DEG2RAD, 2.57e+02 * DEG2RAD), MU_SUN, 398600e9, 6378000, 6900000,  'Earth')"
    """
    self._orig_init(*args)
planet._orig_init = planet.__init__
planet.__init__ = _planet_ctor
