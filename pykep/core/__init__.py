# -*- coding: iso-8859-1 -*-
# We import the proteced symbols we use in this file
from .core import _get_AU, _get_JR, _get_DAY2SEC, _get_DAY2YEAR, _get_DEG2RAD, _get_EARTH_VELOCITY, _get_EARTH_J2, _get_EARTH_RADIUS, _get_MU_EARTH, _get_G0, _get_MU_SUN, _get_RAD2DEG, _get_SEC2DAY
from .core import _epoch_type
# We import symbols we use in this file
from .core import epoch
# We import all symbols in the core namespace (also the ones we do not use
# in this file, but we still want in the namespace core)
from .core import *

"""Defining astronomical constants defined in the keplerian_toolbox file astro_constants.h"""
AU = _get_AU()
JR = _get_JR()
DAY2SEC = _get_DAY2SEC()
DAY2YEAR = _get_DAY2YEAR()
DEG2RAD = _get_DEG2RAD()
EARTH_VELOCITY = _get_EARTH_VELOCITY()
EARTH_J2 = _get_EARTH_J2()
EARTH_RADIUS = _get_EARTH_RADIUS()
G0 = _get_G0()
MU_SUN = _get_MU_SUN()
MU_EARTH = _get_MU_EARTH()
RAD2DEG = _get_RAD2DEG()
SEC2DAY = _get_SEC2DAY()

EPOCH_TYPE = {
    "jd": _epoch_type.JD,
    "mjd": _epoch_type.MJD,
    "mjd2000": _epoch_type.MJD2000
}


def _epoch_ctor(self, julian_date=0, julian_date_type='mjd2000'):
    """
pykep.epoch(julian_date=0, julian_date_type="mjd2000")

- julian_date: a julian date, defaults to 0
- julian_date_type: julian date type, one of "jd", "mjd" and "mjd2000", defaults to "mjd2000"

Examples::

  e1 = epoch(0)
  e1 = epoch(0,"mjd2000")
  e2 = epoch(54333, "mjd")
    """
    self._orig_init(julian_date, EPOCH_TYPE[julian_date_type])

epoch._orig_init = epoch.__init__
epoch.__init__ = _epoch_ctor
