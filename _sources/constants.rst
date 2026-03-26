.. _constants:

Global constants
##################

In `pykep` the access a number of common constants are provided for convenience. The user can overwrite their values if needed.
These constants are not used in the `pykep` internals, they are only provided for convenience for the user to instantiate / use
the various `pykep` objects and functionalities.

.. list-table:: Pykep global constants
   :widths: 50 25 25 50
   :header-rows: 1

   * - Constant's Name
     - Symbol in pykep
     - Units
     - Value
   * - Astronomical Unit 
     - pykep.AU
     - :math:`m` 
     - 1.495978707e+11
   * - Cavendish constant
     - pykep.CAVENDISH
     - :math:`\frac{m^3}{sec^2 kg}` 
     - 6.67430e-11
   * - Sun's gravitational parameter
     - pykep.MU_SUN
     - :math:`\frac{m^3}{sec^2}` 
     - 1.32712440041279419e20
   * - Earth's gravitational parameter (JPL DE440)
     - pykep.MU_EARTH
     - :math:`\frac{m^3}{sec^2}` 
     - 3.98600435507e+14
   * - Moon's gravitational parameter (JPL DE440)
     - pykep.MU_MOON
     - :math:`\frac{m^3}{sec^2}` 
     - 4.902800118e+12
   * - Earth's velocity
     - pykep.EARTH_VELOCITY
     - :math:`\frac{m}{sec}` 
     - 29784.691831696804
   * - Earth's radius
     - pykep.EARTH_RADIUS
     - :math:`m` 
     - 6378137.0
   * - Earth's :math:`J_2`
     - pykep.EARTH_J2
     - --
     - 0.00108262668
   * - CR3BP: Moon-Earth parameter
     - pykep.CR3BP_MU_EARTH_MOON
     - --
     - 0.01215058439470971
   * - BCP: Moon-Earth parameter (same as for the CR3BP)
     - pykep.BCP_MU_EARTH_MOON
     - --
     - 0.01215058439470971
   * - BCP :cite:p:`simo1995bicircular`: Scaled mass of the Sun
     - pykep.BCP_MU_S
     - --
     - 328900.55970856483
   * - BCP :cite:p:`simo1995bicircular`: Scaled Sun–(Earth + Moon) distance
     - pykep.BCP_RHO_S
     - --
     - 3.88811143E2
   * - BCP :cite:p:`simo1995bicircular`: Scaled angular velocity of the Sun
     - pykep.BCP_OMEGA_S
     - --
     - -9.25195985E-01
   * - Seconds in one day
     - pykep.DAY2SEC
     - --
     - 86400.0
   * - Days in one second
     - pykep.SEC2DAY
     - --
     - 1.1574074074074073e-05
   * - Degrees in one radians
     - pykep.RAD2DEG
     - --
     - 57.29577951308232
   * - Radians in one degree
     - pykep.DEG2RAD
     - --
     - 0.017453292519943295
   * - Days in one year
     - pykep.YEAR2DAY
     - --
     - 365.25
   * - Years in one day
     - pykep.DAY2YEAR
     - --
     - 0.0027378507871321013