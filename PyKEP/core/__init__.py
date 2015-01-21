# -*- coding: iso-8859-1 -*-
from _core import *

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

epoch.epoch_type = _core._epoch_type


def read_satcat(satcatfilename=None):
    from collections import namedtuple
    satcatentry = namedtuple('satcatentry', 'intdsgn noradn multnameflag payloadflag operationstatus name ownership launchdate launchsite decay period incl apogee perigee radarA orbitstatus')
    satcat = list()

    with open(satcatfilename, 'r') as f:
        for l1 in f:
            satcat.append(satcatentry(l1[0:11], l1[13:18], l1[19:20], l1[20:21], l1[21:22], l1[23:47], l1[49:54], l1[56:66], l1[68:73], l1[75:85], l1[87:94], l1[96:101], l1[103:109], l1[111:117], l1[119:127], l1[129:132]))
    return satcat


def read_tle(tlefilename, verbose=False):
    """
    This function reads a Two-Line-Element file as taken from the NORAD database
    http://www.celestrak.com/NORAD/elements/ and returns a list of PyKEP planet_tle
    objects

    USAGE: planet_list = read_tle(filename, verbose=True):

    * filename: A string containin the file name (assumed to be in the working directory)
    * verbose: Activates some screen output to show the progress.

    * planet_list: a list of PyKEP.planet_tle objects
    """
    from PyKEP import planet_tle
    planet_list = []

    with open(tlefilename, 'r') as f:
        for line in f:
            line1 = f.next()[:69]
            line2 = f.next()[:69]
            planet_list.append(planet_tle(line1, line2))
    return planet_list
