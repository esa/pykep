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

epoch.epoch_type = _core._epoch_type


def read_satcat(satcatfilename=None):
    """
    This function reads the satcat catalogue, as can be downloaded from http://www.celestrak.com/NORAD/elements/ 
    and returns a dictionary keyed with the sattelite international designator (e.g. "1958-002B"). Each entry
    is a named tuple, see http://celestrak.com/satcat/satcat-format.asp for the meaning of all entries. 
    """
    from collections import namedtuple
    satcatentry = namedtuple('satcatentry', 'noradn multnameflag payloadflag operationstatus name ownership launchdate launchsite decay period incl apogee perigee radarA orbitstatus')
    satcat = dict()

    with open(satcatfilename, 'r') as f:
        for l1 in f:
            intdsgn = l1[0:11].strip()
            satcat[intdsgn] = satcatentry(l1[13:18], l1[19:20], l1[20:21], l1[21:22], l1[23:47], l1[49:54], l1[56:66], l1[68:73], l1[75:85], l1[87:94], l1[96:101], l1[103:109], l1[111:117], l1[119:127], l1[129:132])
    return satcat


def read_tle(tlefilename, verbose=False, with_name=False):
    """
    This function reads a Two-Line-Element file as taken from the NORAD database
    http://www.celestrak.com/NORAD/elements/ and returns a list of PyKEP planet_tle
    objects

    USAGE: planet_list = read_tle(filename, verbose=True):

    * tlefilename: A string containin the file name (assumed to be in the working directory)
    * verbose: Activates some screen output to show the progress.
    * with_name: When True the TLE files does contains satellite names (i.e. when downloaded from space-track instead)

    * planet_list: a list of PyKEP tle_planets having as name their international designator (e.g. "1958-002B").
    """
    from PyKEP import planet_tle
    planet_list = []

    with open(tlefilename, 'r') as f:
        for line in f:
            if with_name:
                line1 = f.next()[:69]
            else:
                line1 = line[:69]
            line2 = f.next()[:69]
            planet_list.append(planet_tle(line1, line2))
    return planet_list
