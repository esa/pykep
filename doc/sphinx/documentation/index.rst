===================
pykep Documentation
===================

pykep is made of main modules called :ref:`core`, :ref:`planet`, :ref:`simsflanagan`, :ref:`pontryagin`, :ref:`phasing`, :ref:`trajopt`, :ref:`orbitplots`, :ref:`util`. When we import pykep via::

  import pykep as pkp

the core module gets imported in the pykep namespace, while all other modules will appear in dedicated namespaces, e.g.
pykep.trajopt, pykep.phasing etc.

In addition, two important dictionaries will be imported::

  pkp.__version__
  pkp.__extensions__

The first one, pkp.__version__, is the pykep version dictionary that can be queried for major, minor or bugfix number. pykep 1.2.0
means that the major version is 1, the minor is 2 and the bugfix version is 0.

The second one, pkp.__extensions__, is the pykep extension dictionary containing the extension installed in the system as detected by pykep.
For example::

  {'matplotlib': True,
 'mplot3d': True,
 'pygmo': True,
 'scikit-learn': True,
 'scipy': True}

Documentation of the main modules
=================================
.. toctree::
  :maxdepth: 1

  core
  simsflanagan
  pontryagin
  phasing
  trajopt
  orbitplots
  planets
  util

Constants defined within pykep namespace
========================================

=================================       =============================================================================
Name                                    Description
=================================       =============================================================================
pykep.AU                                One astronomical unit in meters
pykep.DAY2SEC                           Conversion factor from days to seconds
pykep.DAY2YEAR                          Conversion factor from days to years
pykep.DEG2RAD                           Conversion factor from degrees to radians
pykep.RAD2DEG                           Conversion factor from radians to degree
pykep.SEC2DAY                           Conversion factor from seconds to days
pykep.MU_SUN                            Sun gravitational constant in m^3/s^2
pykep.EARTH_VELOCITY                    Square root of MU_SUN/AU. Average earth velocity in meters per seconds.
pykep.G0                                The standard gravity acceleration at ground level in m/s^2
=================================       =============================================================================


.. note::

   pykep follows two conventions that are to be kept in mind.

   * :emphasis:`Units`: The S.I. system is used to perform all computations in pykep and it is the default choice for everything. In some cases, to help the user instantiate objects and define problems, classes and functions may accept inputs in different units or return results in different units, but this is only exceptional and documented in the class or function help.

   * :emphasis:`Orbital parameters`: In pykep, the osculating Keplerian elements are always in the following order: a,e,i,W,w,M (E or H), where a is the semi-major axis, e the eccentriciy, i the inclination, W the Right Axcension of the Acending Node, w the argument of perigee and M the mean anomaly (E and H being the eccentric and hyperbolic anomalies).
