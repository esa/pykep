# -*- coding: iso-8859-1 -*-
"""
This module contains all the classes that allow to construct efficiently
low-thrust tajectories using our own flavour of the Sims-Flanagan model: a trajectory
transcription method that forms the basis for MALTO, the software in use in JPL
for preliminary interplanetary trajectory design.
"""
from pykep.util._util import *


def read_satcat(satcatfilename=None):
    """
    pykep.read_satcat("my_dir/satcat.txt")

    This function reads the satcat catalogue, as can be downloaded from http://www.celestrak.com/NORAD/elements/
    and returns a dictionary keyed with the satellite international designator (e.g. "1958-002B"). Each entry
    is a named tuple with the following fields::

      noradn, multnameflag, payloadflag, operationstatus, name, ownership, launchdate, launchsite, decay, period, incl, apogee, perigee, radarA, orbitstatus

    see http://celestrak.com/satcat/satcat-format.asp for the meaning of all entries.

    Example::

      satcat = pykep.read_satcat("my_dir/satcat.txt")
      cross_section = satcat["1958-002B"].radarA
    """
    from collections import namedtuple
    satcatentry = namedtuple(
        'satcatentry', 'noradn multnameflag payloadflag operationstatus name ownership launchdate launchsite decay period incl apogee perigee radarA orbitstatus')
    satcat = dict()

    with open(satcatfilename, 'r') as f:
        for l1 in f:
            intdsgn = l1[0:11].strip()
            satcat[intdsgn] = satcatentry(l1[13:18], l1[19:20], l1[20:21], l1[21:22], l1[23:47], l1[49:54], l1[
                                          56:66], l1[68:73], l1[75:85], l1[87:94], l1[96:101], l1[103:109], l1[111:117], l1[119:127], l1[129:132])
    return satcat


def read_tle(tle_file, verbose=False, with_name=True):
    """
    planet_list = pykep.read_tle(tle_file, verbose=False, with_name=True)

    - tle_file: A string containin the file name (assumed to be in the working directory)
    - verbose: Activates some screen output to show the progress.
    - with_name: When True the TLE files contain three lines, the first one containing the satellite common name

    - [out] planet_list: a list of pykep tle_planets having as name their international designator (e.g. "1958-002B").

    This function reads a Two-Line-Element file as taken from the NORAD database
    http://www.celestrak.com/NORAD/elements/ or equivalent and returns a list of pykep planet_tle
    objects

    .. note::

       The name of each of the instantiated planets will be its international designator and
       thus a valid key to the satcat dictionary

    Example::

      planet_list = pykep.read_tle("my_dir/cosmos.tle")
      satcat = pykep.read_satcat("my_dir/satcat.txt")
      cross_sections = [satcat[pl.name()].radarA for pl in  planet_list]
    """
    from pykep.planet import tle
    planet_list = []

    with open(tle_file, 'r') as f:
        for line in f:
            if with_name:
                line1 = next(f)[:69]
            else:
                line1 = line[:69]
            line2 = next(f)[:69]
            planet_list.append(tle(line1, line2))
    return planet_list
