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


def read_satcat(satcatfilename=None):
    """
    PyKEP.read_satcat("my_dir/satcat.txt")

    This function reads the satcat catalogue, as can be downloaded from http://www.celestrak.com/NORAD/elements/
    and returns a dictionary keyed with the satellite international designator (e.g. "1958-002B"). Each entry
    is a named tuple with the following fields::

      noradn, multnameflag, payloadflag, operationstatus, name, ownership, launchdate, launchsite, decay, period, incl, apogee, perigee, radarA, orbitstatus

    see http://celestrak.com/satcat/satcat-format.asp for the meaning of all entries.

    Example::

      satcat = PyKEP.read_satcat("my_dir/satcat.txt")
      cross_section = satcat["1958-002B"].radarA
    """
    from collections import namedtuple
    satcatentry = namedtuple('satcatentry', 'noradn multnameflag payloadflag operationstatus name ownership launchdate launchsite decay period incl apogee perigee radarA orbitstatus')
    satcat = dict()

    with open(satcatfilename, 'r') as f:
        for l1 in f:
            intdsgn = l1[0:11].strip()
            satcat[intdsgn] = satcatentry(l1[13:18], l1[19:20], l1[20:21], l1[21:22], l1[23:47], l1[49:54], l1[56:66], l1[68:73], l1[75:85], l1[87:94], l1[96:101], l1[103:109], l1[111:117], l1[119:127], l1[129:132])
    return satcat


def read_tle(tle_file, verbose=False, with_name=False):
    """
    planet_list = PyKEP.read_tle(tle_file, verbose=False, with_name=False)

    - tle_file: A string containin the file name (assumed to be in the working directory)
    - verbose: Activates some screen output to show the progress.
    - with_name: When True the TLE files contain three lines, the first one containing the satellite common name

    - [out] planet_list: a list of PyKEP tle_planets having as name their international designator (e.g. "1958-002B").

    This function reads a Two-Line-Element file as taken from the NORAD database
    http://www.celestrak.com/NORAD/elements/ or equivalent and returns a list of PyKEP planet_tle
    objects

    .. note::

       The name of each of the instantiated planets will be its international designator and
       thus a valid key to the satcat dictionary

    Example::

      planet_list = PyKEP.read_tle("my_dir/cosmos.tle")
      satcat = PyKEP.read_satcat("my_dir/satcat.txt")
      cross_sections = [satcat[pl.name()].radarA for pl in  planet_list]
    """
    from PyKEP import planet_tle
    planet_list = []

    with open(tle_file, 'r') as f:
        for line in f:
            if with_name:
                line1 = f.next()[:69]
            else:
                line1 = line[:69]
            line2 = f.next()[:69]
            planet_list.append(planet_tle(line1, line2))
    return planet_list
