======================
Code Documentation
======================

pykep is made of main modules called :ref:`core`, :ref:`planet`, :ref:`simsflanagan`, :ref:`pontryagin`, :ref:`phasing`, :ref:`trajopt`, :ref:`orbitplots`, :ref:`util`. 

When we import pykep via::

  import pykep as pk

the core module gets imported in the pykep namespace, while all other modules will appear in dedicated namespaces, e.g.
pykep.trajopt, pykep.phasing etc.

In addition, two dictionaries will be imported::

  pk.__version__
  pk.__extensions__

The first one, pk.__version__, is the pykep version dictionary that can be queried for major and minor version.

The second one, pk.__extensions__, is the pykep extension dictionary containing relevant third-party dependencies
as detected by pykep.

For example::

  {'matplotlib': True,
 'mplot3d': True,
 'pygmo': True,
 'scikit-learn': True,
 'scipy': True}

----------------------------------------------------------------------

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
pykep.MU_EARTH                          Earth gravitational constant in m^3/s^2
pykep.EARTH_VELOCITY                    Square root of MU_SUN/AU. Average earth velocity in meters per seconds.
pykep.G0                                The standard gravity acceleration at ground level in m/s^2
=================================       =============================================================================


.. note::

   pykep follows two conventions that are to be kept in mind.

   * :emphasis:`Units`: The S.I. system is used to perform all computations in pykep and it is the default choice for everything. In some cases, to help the user instantiate objects and define problems, classes and functions may accept inputs in different units or return results in different units, but this is only exceptional and documented in the class or function help.

   * :emphasis:`Orbital parameters`: In pykep, the osculating Keplerian elements are always in the following order: a,e,i,W,w,M (E or H), where a is the semi-major axis, e the eccentriciy, i the inclination, W the Right Axcension of the Acending Node, w the argument of perigee and M the mean anomaly (E and H being the eccentric and hyperbolic anomalies).

----------------------------------------------------------------------

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
