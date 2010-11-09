# -*- coding: iso-8859-1 -*-
from _PyKEP import *

"""Extracts astronomical constants defined in the keplerian_toolbox file astro_constants.h"""
AU = _PyKEP._get_AU()
DAY2SEC = _PyKEP._get_DAY2SEC()
DAY2YEAR = _PyKEP._get_DAY2YEAR()
DEG2RAD = _PyKEP._get_DEG2RAD()
EARTH_VELOCITY = _PyKEP._get_EARTH_VELOCITY()
G0 = _PyKEP._get_G0()
MU_SUN = _PyKEP._get_MU_SUN()
RAD2DEG = _PyKEP._get_RAD2DEG()
SEC2DAY = _PyKEP._get_SEC2DAY()

epoch.epoch_type = _PyKEP._epoch_type