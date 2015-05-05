# -*- coding: iso-8859-1 -*-
from PyKEP.core._core import _get_AU, _get_JR, _get_DAY2SEC, _get_DAY2YEAR, _get_DEG2RAD, _get_EARTH_VELOCITY, _get_G0, _get_MU_SUN, _get_RAD2DEG, _get_SEC2DAY
from PyKEP.core._core import _epoch_type
from PyKEP.core._core import *

"""Defining astronomical constants defined in the keplerian_toolbox file astro_constants.h"""
AU = _get_AU()
JR = _get_JR()
DAY2SEC = _get_DAY2SEC()
DAY2YEAR = _get_DAY2YEAR()
DEG2RAD = _get_DEG2RAD()
EARTH_VELOCITY = _get_EARTH_VELOCITY()
G0 = _get_G0()
MU_SUN = _get_MU_SUN()
RAD2DEG = _get_RAD2DEG()
SEC2DAY = _get_SEC2DAY()

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
