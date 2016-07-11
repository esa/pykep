# -*- coding: iso-8859-1 -*-
"""
This module contains all the classes that allow to construct efficiently
low-thrust tajectories using our own flavour of the Sims-Flanagan model: a trajectory
transcription method that forms the basis for MALTO, the software in use in JPL
for preliminary interplanetary trajectory design.
"""
from PyKEP.planet._planet import *
from PyKEP.planet._planet import _base


def _keplerian_ctor(self, *args):
    """
PyKEP.planet.keplerian(when,orbital_elements, mu_central_body, mu_self,radius, safe_radius [, name = 'unknown'])

PyKEP.planet.keplerian(when,r,v, mu_central_body, mu_self,radius, safe_radius [, name = 'unknown'])

- when: a :py:class:`PyKEP.epoch` indicating the orbital elements reference epoch
- orbital_elements: a sequence of six containing a,e,i,W,w,M (SI units, i.e. meters and radiants)
- r,v: position and velocity of an object at when (SI units)
- mu_central_body: gravity parameter of the central body (SI units, i.e. m^2/s^3)
- mu_self: gravity parameter of the planet (SI units, i.e. m^2/s^3)
- radius: body radius (SI units, i.e. meters)
- safe_radius: mimimual radius that is safe during a fly-by of the planet (SI units, i.e. m)
- name: body name

.. note::

   All classes having Keplerian ephemerides as :py:class:`PyKEP.planet.mpcorb` inherit from this (c++) class

Example::

  earth = planet(epoch(54000,"mjd"),(9.99e-01 * AU, 1.67e-02, 8.85e-04 * DEG2RAD, 1.75e+02 * DEG2RAD, 2.87e+02 * DEG2RAD, 2.57e+02 * DEG2RAD), MU_SUN, 398600e9, 6378000, 6900000,  'Earth')"
    """
    self._orig_init(*args)
keplerian._orig_init = keplerian.__init__
keplerian.__init__ = _keplerian_ctor

def _j2_ctor(self, *args):
    """
PyKEP.planet.j2(when,orbital_elements, mu_central_body, mu_self,radius, safe_radius, J2RG2 [, name = 'unknown'])

PyKEP.planet.j2(when,r,v, mu_central_body, mu_self,radius, safe_radius, J2RG2 [, name = 'unknown'])

- when: a :py:class:`PyKEP.epoch` indicating the orbital elements reference epoch
- orbital_elements: a sequence of six containing a,e,i,W,w,M (SI units, i.e. meters and radiants)
- r,v: position and velocity of an object at when (SI units)
- mu_central_body: gravity parameter of the central body (SI units, i.e. m^2/s^3)
- mu_self: gravity parameter of the planet (SI units, i.e. m^2/s^3)
- radius: body radius (SI units, i.e. meters)
- safe_radius: mimimual radius that is safe during a fly-by of the planet (SI units, i.e. m)
- J2RG2: the product of J2 and the mean radius of the oblate primary
- name: body name

Example::

  deb1 = planet.j2(epoch(54000,"mjd"),(7000000, 1.67e-02, 78.23 * DEG2RAD, 175. * DEG2RAD, 287. * DEG2RAD, 257 * DEG2RAD), MU_EARTH, 1, 1, 1, J2*EARTH_RADIUS**2 'deb1')"
    """
    self._orig_init(*args)
j2._orig_init = j2.__init__
j2.__init__ = _j2_ctor
