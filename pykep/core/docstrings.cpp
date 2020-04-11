#include <string>

#include "docstrings.hpp"

namespace pykep
{

std::string epoch_doc()
{
    return R"(
A precise point in time. Julian dates are supported.
)";
}

std::string epoch_from_string_doc()
{
    return R"(
pykep.epoch_from_string(s)

- s: a string containing the date in the format 'YYYY-MM-DD HH:MM:SS'

Returns a :py:class:`pykep.epoch`

Note::

    Excess digits in fractional seconds will be dropped. 
    Ex: '1:02:03.123456999' => '1:02:03.123456'. This behavior depends on the precision
    defined in astro_constant.hpp and used to compile the keplerian toolbox library.
    
The function is based on the corresponding `boost date_time library function <http://www.boost.org/doc/libs/1_71_0/doc/html/date_time/posix_time.html#ptime_from_string>`_

Example::

    e = pykep.epoch_from_string('2002-01-20 23:59:54.003')
)";
}

std::string epoch_from_iso_string_doc()
{
    return R"(
pykep.epoch_from_iso_string(s)

- s: string containing a date in the ISO format 'YYYYMMDDTHHMMSS'

Returns a :py:class:`pykep.epoch` object constructed from a from a non delimited iso form string containing a date.

The function is based on the corresponding `boost date_time library function <http://www.boost.org/doc/libs/1_71_0/doc/html/date_time/posix_time.html#ptime_from_string>`_

Example::

    e = pykep.epoch_from_iso_string('20020120T235954')
)";
}

std::string lambert_problem_doc()
{
    return R"(
lambert_problem(r1 = [1,0,0], r2 = [0,1,0], tof = pi/2, mu = 1., cw = False, max_revs = 0)

- r1: starting position (x1,y1,z1)
- r2: final position    (x2,y2,z2)
- tof: time of flight
- mu: gravitational parameter (default is 1)
- cw: True for retrograde motion (clockwise), False if counter-clock wise
- max_revs: Maximum number of multirevs to be computed

.. note::

   Units need to be consistent. The multirev Lambert's problem will be solved upon construction and its solution stored in data members.

Example:: 

    l = lambert_problem(r1 = [1,0,0], r2 = [0,1,0], tof = 5 * pi / 2.)
)";
}

std::string propagate_lagrangian_doc()
{
    return R"(

propagate_lagrangian(r0 = [1,0,0], v0 = [0,1,0], tof = pi/2, mu = 1)

- r: start position, x,y,z
- v: start velocity, vx,vy,vz
- tof: propagation time
- mu: central body gravity constant

Returns a tuple (rf, vf) containing the final position and velocity after the propagation.

Example::

  rf,vf = propagate_lagrangian(r0 = [1,0,0], v0 = [0,1,0], tof = pi/2, mu = 1)
)";
}

std::string propagate_taylor_doc()
{
    return R"(

pykep.propagate_taylor(r0 = [1,0,0], v0 = [0,1,0], m0 = 100, thrust = [0,0,0], tof = pi/2, mu = 1, veff = 1, log10tol =-15, log10rtol = -15)

- r: start position, x,y,z.
- v: start velocity, vx,vy,vz.
- m0: starting mass.
- thrust: fixed inertial thrust, ux,uy,uz.
- tof: propagation time.
- mu: central body gravity constant.
- veff: the product (Isp g0) defining the engine efficiency.
- log10tol: the logarithm of the absolute tolerance passed to taylor propagator.
- log10rtol: the logarithm of the relative tolerance passed to taylor propagator.

Returns a tuple (rf, vf, mf) containing the final position, velocity and mass after the propagation.

Example::

  r,v,m = propagate_taylor(r0 = [1,0,0], v0 = [0,1,0], m0 = 100, thrust = [0,0,0], tof = pi/2, mu = 1, veff = 1, log10tol =-15, log10rtol = -15)
)";
}





} // namespace pykep
