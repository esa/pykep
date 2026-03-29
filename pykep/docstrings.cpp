// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <string>

#include "docstrings.hpp"
#include "kep3/core_astro/basic_transfers.hpp"

namespace pykep
{

std::string core_module_doc()
{
    return R"(core is the Pykep module that contains most of its core routines efficiently coded in c++
)";
}

std::string m2e_doc()
{
    return R"(m2e(M, ecc)
    
    Converts from Mean to Eccentric anomaly. Requires ecc < 1.

    Args:
      *M* (:class:`float`): the Mean anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> M = 1.2
      >>> ecc = 0.1
      >>> pk.m2e(M, ecc)
      1.296254963787226
)";
}

std::string e2m_doc()
{
    return R"(e2m(E, ecc)
    
    Converts from Eccentric to Mean anomaly. Requires ecc < 1.

    Args:
      *E* (:class:`float`): the Eccentric anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> E = 0.5
      >>> ecc = 0.1
      >>> pk.e2m(E, ecc)
      0.4520574461395797
)";
}

std::string e2f_doc()
{
    return R"(e2f(E, ecc)
    
    Converts from eccentric to true anomaly. Requires ecc < 1.

    Args:
      *E* (:class:`float`): the Eccentric anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> E = 0.5
      >>> ecc = 0.1
      >>> pk.e2f(E, ecc)
      0.5502639747136633
)";
}

std::string f2e_doc()
{
    return R"(f2e(f, ecc)
    
    Converts from True to Eccentric anomaly. Requires ecc < 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> f = 1.2
      >>> ecc = 0.1
      >>> pk.f2e(f, ecc)
      1.1082931139529482
)";
}

std::string f2m_doc()
{
    return R"(f2m(f, ecc)
    
    Converts from True to Mean anomaly. Requires ecc < 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Mean anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> f = -0.34
      >>> ecc = 0.67
      >>> pk.f2m(f, ecc)
      -0.05065883735669101
)";
}

std::string m2f_doc()
{
    return R"(m2f(M, ecc)
    
    Converts from Mean to True anomaly. Requires ecc < 1.

    Args:
      *M* (:class:`float`): the Mean anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> M = 0.32
      >>> ecc = 0.65
      >>> pk.m2f(M, ecc)
      1.4497431281728277
)";
}

std::string h2n_doc()
{
    return R"(h2n(H, ecc)
    
    Converts from Hyperbolic to Hyperbolic Mean anomaly. Requires ecc > 1.

    Args:
      *H* (:class:`float`): the Hyperbolic anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Hyperbolic Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> H = 1.2
      >>> ecc = 10.32
      >>> pk.h2n(H, ecc)
      14.377641187853621
)";
}

std::string n2h_doc()
{
    return R"(n2h(N, ecc)
    
    Converts from Hyperbolic Mean to Hyperbolic anomaly. Requires ecc > 1.

    Args:
      *N* (:class:`float`): the Hyperbolic Mean anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Hyperbolic anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> N = 1.2
      >>> ecc = 10.32
      >>> pk.n2h(N, ecc)
      0.12836469743916526
)";
}

std::string h2f_doc()
{
    return R"(h2f(H, ecc)
    
    Converts from Hyperbolic to True anomaly. Requires ecc > 1.

    Args:
      *H* (:class:`float`): the Hyperbolic anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> H = 10.32
      >>> ecc = 4.5
      >>> pk.h2f(H, ecc)
      1.7948251330114304
)";
}

std::string f2h_doc()
{
    return R"(f2h(f, ecc)
    
    Converts from True to Hyperbolic anomaly. Requires ecc > 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Hyperbolic anomaly

    Examples:
      >>> import pykep as pk
      >>> f = 1.2
      >>> ecc = 1.1
      >>> pk.f2h(f, ecc)
      0.30083016696826936
)";
}

std::string f2n_doc()
{
    return R"(f2n(f, ecc)
    
    Converts from True to Hyperbolic Mean anomaly. Requires ecc > 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Hyperbolic Mean anomaly

    Examples:
      >>> import pykep as pk
      >>> f = 1.2
      >>> ecc = 5.7
      >>> pk.f2n(f, ecc)
      8.421335633880908
)";
}

std::string n2f_doc()
{
    return R"(n2f(N, ecc)
    
    Converts from Hyperbolic Mean to True anomaly. Requires ecc > 1.

    Args:
      *N* (:class:`float`): the Hyperbolic Mean anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> N = 10.32
      >>> ecc = 13.45
      >>> pk.n2f(N, ecc)
      0.7373697968359353
)";
}

std::string zeta2f_doc()
{
    return R"(zeta2f(zeta, ecc)
    
    Converts from Gudermannian to True anomaly. Requires ecc > 1.

    See Battin: "An Introduction to the Mathematics and Methods of Astrodynamics" for a 
    definition of zeta and the treatment of the resulting equations.

    Args:
      *zeta* (:class:`float`): the Gudermannian (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> zeta = 8.2
      >>> ecc = 2.2
      >>> pk.zeta2f(zeta, ecc)
      2.3290929552114266
)";
}

std::string f2zeta_doc()
{
    return R"(f2zeta(f, ecc)
    
    Converts from True anomaly to Gudermannian. Requires ecc > 1.

    Args:
      *f* (:class:`float`): the True anomaly (rad.)

      *ecc* (:class:`float`): the eccentricity

    Returns:
      :class:`float`: the Gudermannian 

    Examples:
      >>> import pykep as pk
      >>> f = 0.5
      >>> ecc = 3.3
      >>> pk.f2zeta(f, ecc)
      0.36923933496389816
)";
}

std::string m2e_v_doc()
{
    return R"(m2e_v(Ms, eccs)
    
    Converts from Mean to Eccentric anomaly (vectorized version). Requires ecc < 1.

    Args:
      *Ms* (:class:`numpy.ndarray` or :class:`float`): the Mean anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Ms = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.375
      >>> Es = pk.m2e_v(Ms, ecc)
      >>> np.shape(Es)
      (100,)
)";
}

std::string e2m_v_doc()
{
    return R"(e2m_v(Es, eccs)
    
    Converts from Eccentric to Mean anomaly (vectorized version). Requires ecc < 1.

    Args:
      *Es* (:class:`numpy.ndarray` or :class:`float`): the Eccentric anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Es = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.86345
      >>> Ms = pk.e2m_v(Es, ecc)
      >>> np.shape(Ms)
      (100,)
)";
}

std::string e2f_v_doc()
{
    return R"(e2f_v(Es, eccs)
    
    Converts from eccentric to true anomaly (vectorized version). Requires ecc < 1.

    Args:
      *Es* (:class:`numpy.ndarray` or :class:`float`): the Eccentric anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Es = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.0256
      >>> fs = pk.e2f_v(Es, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string f2e_v_doc()
{
    return R"(f2e_v(fs, eccs)
    
    Converts from True to Eccentric anomaly (vectorized version). Requires ecc < 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Eccentric anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.23
      >>> Es = pk.f2e_v(fs, ecc)
      >>> np.shape(Es)
      (100,)
)";
}

std::string f2m_v_doc()
{
    return R"(f2m_v(fs, eccs)
    
    Converts from True to Mean anomaly (vectorized version). Requires ecc < 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Mean anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.4
      >>> Ms = pk.f2m_v(fs, ecc)
      >>> np.shape(Ms)
      (100,)
)";
}

std::string m2f_v_doc()
{
    return R"(m2f_v(Ms, eccs)
    
    Converts from Mean to True anomaly (vectorized version). Requires ecc < 1.

    Args:
      *Ms* (:class:`numpy.ndarray` or :class:`float`): the Mean anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Ms = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 0.4
      >>> fs = pk.m2f_v(Ms, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string h2n_v_doc()
{
    return R"(h2n_v(Hs, eccs)
    
    Converts from Hyperbolic to Hyperbolic Mean anomaly (vectorized version). Requires ecc > 1.

    Args:
      *Hs* (:class:`numpy.ndarray` or :class:`float`): the Hyperbolic anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Hyperbolic Mean anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Hs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 4.5
      >>> Ns = pk.h2n_v(Hs, ecc)
      >>> np.shape(Ns)
      (100,)
)";
}

std::string n2h_v_doc()
{
    return R"(n2h_v(Ns, eccs)
    
    Converts from Hyperbolic Mean to Hyperbolic anomaly (vectorized version). Requires ecc > 1.

    Args:
      *Ns* (:class:`numpy.ndarray` or :class:`float`): the Hyperbolic Mean anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Hyperbolic anomaly (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Ns = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 4.5
      >>> Hs = pk.n2h_v(Ns, ecc)
      >>> np.shape(Hs)
      (100,)
)";
}

std::string h2f_v_doc()
{
    return R"(h2f_v(Hs, eccs)
    
    Converts from Hyperbolic to True anomaly (vectorized version). Requires ecc > 1.

    Args:
      *Hs* (:class:`numpy.ndarray` or :class:`float`): the Hyperbolic anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly in [-pi, pi] (rad.)

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Hs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 4.5
      >>> fs = pk.h2f_v(Hs, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string f2h_v_doc()
{
    return R"(f2h_v(fs, eccs)
    
    Converts from True to Hyperbolic anomaly (vectorized version). Requires ecc > 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Hyperbolic anomaly

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 5.7
      >>> Hs = pk.n2h_v(fs, ecc)
      >>> np.shape(Hs)
      (100,)
)";
}

std::string f2n_v_doc()
{
    return R"(f2n_v(fs, eccs)
    
    Converts from True to Hyperbolic Mean anomaly (vectorized version). Requires ecc > 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Hyperbolic Mean anomaly

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 5.7
      >>> Ns = pk.n2f_v(fs, ecc)
      >>> np.shape(Ns)
      (100,)
)";
}

std::string n2f_v_doc()
{
    return R"(n2f_v(Ns, eccs)
    
    Converts from Hyperbolic Mean to True anomaly (vectorized version). Requires ecc > 1.

    Args:
      *Ns* (:class:`numpy.ndarray` or :class:`float`): the Hyperbolic Mean anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> Ns = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 13.45
      >>> fs = pk.n2f_v(Ns, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string zeta2f_v_doc()
{
    return R"(zeta2f_v(zetas, eccs)
    
    Converts from Gudermannian to True anomaly (vectorized version). Requires ecc > 1.

    See Battin: "An Introduction to the Mathematics and Methods of Astrodynamics" for a 
    definition of zeta and the treatment of the resulting equations.

    Args:
      *zetas* (:class:`numpy.ndarray` or :class:`float`): the Gudermannian (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the True anomaly

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> zetas = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 2.2
      >>> fs = pk.zeta2f_v(zetas, ecc)
      >>> np.shape(fs)
      (100,)
)";
}

std::string f2zeta_v_doc()
{
    return R"(f2zeta_v(fs, eccs)
    
    Converts from True anomaly to Gudermannian (vectorized version). Requires ecc > 1.

    Args:
      *fs* (:class:`numpy.ndarray` or :class:`float`): the True anomaly (rad.)

      *eccs* (:class:`numpy.ndarray` or :class:`float`): the eccentricity

    Returns:
      :class:`numpy.ndarray` or :class:`float`: the Gudermannian 

    Examples:
      >>> import pykep as pk
      >>> import numpy as np
      >>> fs = np.linspace(-np.pi/2, np.pi/2, 100)
      >>> ecc = 10.2
      >>> zetas = pk.f2zeta_v(fs, ecc)
      >>> np.shape(zetas)
      (100,)
)";
}

std::string hohmann_doc()
{
    return R"(hohmann(r1, r2, mu)

    Computes the delta v and transfer time required for a Hohmann transfer between two circular orbits.

    Args:
        *r1* (:class:`float`): radius of the first orbit

        *r2* (:class:`float`): radius of the second orbit

        *mu* (:class:`float`): gravitational parameter of the central body

    Returns:
        [:class:`float`, :class:`float`, [:class:`float`, :class:`float`]]:
            [total delta v, transfer time, [first delta v, second delta v]]

    Examples:
      >>> import pykep as pk
      >>> r1 = 7000000
      >>> r2 = 9000000
      >>> mu = pk.MU_EARTH
      >>> dv_total, t_transfer, [dv1, dv2] = pk.hohmann(r1, r2, mu)
      >>> print("Total delta v:", dv_total, "m/s")
      >>> print("Transfer time:", t_transfer, "s")
      >>> print("First delta v:", dv1, "m/s")
      >>> print("Second delta v:", dv2, "m/s")
)";
}

std::string ic2par_doc()
{
    return R"(ic2par(posvel, mu)

    Converts Cartesian state vectors (position and velocity) to Keplerian osculating orbital elements.

    Args:
        *posvel* (:class:`list` [:class:`list`, :class:`list`]): A list containing two 3D vectors: the position vector [x, y, z] in units L and the velocity vector [vx, vy, vz] in units L/T.

        *mu* (:class:`float`): Gravitational parameter of the central body (in units L^3/T^2).

    Returns:
        [:class:`float`, :class:`float`, :class:`float`, :class:`float`, :class:`float`, :class:`float`]:
            A list of six Keplerian orbital elements:

            - *a*: semi-major axis (in units L, positive for ellipses, negative for hyperbolae)
            - *e*: eccentricity (unitless)
            - *i*: inclination (radians, in [0, π])
            - *Ω*: RAAN (radians, in [0, 2π])
            - *ω*: argument of periapsis (radians, in [0, 2π])
            - *f*: true anomaly (radians, in [0, 2π])

    Examples:
      >>> import pykep as pk
      >>> r = [7000e3, 0, 0]
      >>> v = [0, 7.5e3, 0]
      >>> mu = pk.MU_EARTH
      >>> elements = pk.ic2par([r, v], mu)
      >>> a, e, i, Omega, omega, f = elements
      >>> print(f"Semi-major axis: {a} m, Eccentricity: {e}")
)";
}

std::string par2ic_doc()
{
    return R"(par2ic(par, mu)

    Converts Keplerian osculating orbital elements to Cartesian state vectors (position and velocity).

    Args:
        *par* (:class:`list` of :class:`float`): A list of six Keplerian orbital elements:

            - *a*: semi-major axis (in units L, positive for ellipses, negative for hyperbolae)
            - *e*: eccentricity (unitless)
            - *i*: inclination (radians, in [0, π])
            - *Ω*: longitude of ascending node (radians, in [0, 2π])
            - *ω*: argument of periapsis (radians, in [0, 2π])
            - *f*: true anomaly (radians, in [0, 2π])

        *mu* (:class:`float`): 
            Gravitational parameter of the central body (in units L^3/T^2)

    Returns:
        [:class:`list` of :class:`float`, :class:`list` of :class:`float`]:
            A list containing two 3D vectors:

            - *position*: Cartesian position vector [x, y, z] in L  
            - *velocity*: Cartesian velocity vector [vx, vy, vz] in L/T

    Raises:
        :class:`ValueError`: If the semi-major axis and eccentricity are incompatible  
        :class:`ValueError`: If the true anomaly is beyond the asymptotes for a hyperbolic trajectory

    Examples:
      >>> import pykep as pk
      >>> a = 10000e3
      >>> e = 0.1
      >>> i = 0.1
      >>> Omega = 0.5
      >>> omega = 1.0
      >>> f = 2.0
      >>> mu = pk.MU_EARTH
      >>> r, v = pk.par2ic([a, e, i, Omega, omega, f], mu)
      >>> print("Position vector:", r)
      >>> print("Velocity vector:", v)
)";
}

std::string ic2mee_doc()
{
    return R"(ic2mee(posvel, mu, retrogade)

    Converts Cartesian state vectors (position and velocity) to equinoctial orbital elements.

    Equinoctial elements provide a non-singular representation of orbital motion, especially useful for
    near-circular or near-equatorial orbits. The retrograde flag allows switching between the standard and
    retrograde equinoctial elements. These last are not singular for inclinations of π.

    Args:
        *posvel* (:class:`list` [:class:`list`, :class:`list`]): A list containing two 3D vectors: the position vector [x, y, z] in units L and the velocity vector [vx, vy, vz] in units L/T.

        *mu* (:class:`float`): Gravitational parameter of the central body (in units L^3/T^2).

        *retrogade* (:class:`bool`): Whether to use the retrograde equinoctial frame.

    Returns:
        :class:`list` of :class:`float`:
            A list of six equinoctial orbital elements:

            - *p*: semi-latus rectum (in units L)  
            - *f*: eccentricity vector times cos(Ω+ω)
            - *g*: eccentricity vector times sin(Ω+ω)
            - *h*: tan(i/2) cos Ω
            - *k*: tan(i/2) sin Ω
            - *L*: true longitude (radians, in [0, 2π])

    Examples:
      >>> import pykep as pk
      >>> r = [7000e3, 0.0, 0.0]
      >>> v = [0.0, 7.5e3, 1.0e3]
      >>> mu = pk.MU_EARTH
      >>> retro = False
      >>> eq = pk.ic2mee([r, v], mu, retro)
      >>> print("Equinoctial elements:", eq)
)";
}

std::string mee2ic_doc()
{
    return R"(mee2ic(eq_elem, mu, retrogade)

    Converts equinoctial orbital elements to Cartesian state vectors (position and velocity).

    Equinoctial elements provide a non-singular representation of orbital motion, especially useful for
    near-circular or near-equatorial orbits. The retrograde flag allows switching between the standard and
    retrograde equinoctial frames. These last are not singular for inclinations of \pi.

    Args:
        *eq_elem* (:class:`list` [:class:`float`]): A list of six equinoctial elements:
            - *p*: semi-latus rectum (in units L)  
            - *f*: eccentricity vector times cos(Ω+ω)
            - *g*: eccentricity vector times sin(Ω+ω)
            - *h*: tan(i/2) cos Ω
            - *k*: tan(i/2) sin Ω
            - *L*: true longitude (radians, in [0, 2π])

        *mu* (:class:`float`): Gravitational parameter of the central body (in units L^3/T^2).

        *retrogade* (:class:`bool`): Whether to use the retrograde equinoctial frame.

    Returns:
        :class:`list` of :class:`list`:
            A list containing two 3D vectors:

            - *position* vector [x, y, z] in units L
            - *velocity* vector [vx, vy, vz] in units L/T

    Examples:
      >>> import pykep as pk
      >>> eq = [7000e3, 0.01, 0.01, 0.01, 0.01, 0.0]
      >>> mu = pk.MU_EARTH
      >>> retro = False
      >>> r, v = pk.mee2ic(eq, mu, retro)
      >>> print("Position:", r)
      >>> print("Velocity:", v)
)";
}

std::string mee2par_doc()
{
    return R"(mee2par(eq_elem, retrogade)

    Converts equinoctial orbital elements to classical Keplerian elements.

    This function transforms the non-singular equinoctial elements into classical orbital elements, which are 
    more intuitive but can be singular for certain inclinations or eccentricities. The retrograde flag selects
    the appropriate transformation for orbits with inclination near \pi.

    Args:
        *eq_elem* (:class:`list` [:class:`float`]): A list of six equinoctial elements:
            - *p*: semi-latus rectum (in units L)  
            - *f*: eccentricity vector times cos(Ω+ω)
            - *g*: eccentricity vector times sin(Ω+ω)
            - *h*: tan(i/2) cos Ω
            - *k*: tan(i/2) sin Ω
            - *L*: true longitude (radians, in [0, 2π])

        *retrogade* (:class:`bool`): Whether to use the retrograde equinoctial frame.

    Returns:
        :class:`list` of :class:`float`:
            A list of six Keplerian orbital elements:

            - *a*: semi-major axis (in units L, positive for ellipses, negative for hyperbolae)
            - *e*: eccentricity (unitless)
            - *i*: inclination (radians, in [0, π])
            - *Ω*: longitude of ascending node (radians, in [0, 2π])
            - *ω*: argument of periapsis (radians, in [0, 2π])
            - *f*: true anomaly (radians, in [0, 2π])

    Examples:
      >>> eq = [7000e3, 0.01, 0.02, 0.001, 0.002, 0.5]
      >>> retro = False
      >>> par = mee2par(eq, retro)
      >>> print("Keplerian elements:", par)
)";
}

std::string par2mee_doc()
{
    return R"(par2mee(par, retrogade)

    Converts classical Keplerian orbital elements to equinoctial elements.

    This function provides a non-singular representation of orbits by transforming Keplerian elements into 
    equinoctial elements. The retrograde flag allows conversion to a frame that remains non-singular for 
    inclinations near \pi.

    Args:
        *par* (:class:`list` [:class:`float`]): A list of six Keplerian elements:
            - *a*: semi-major axis (in units L, positive for ellipses, negative for hyperbolae)
            - *e*: eccentricity (unitless)
            - *i*: inclination (radians, in [0, π])
            - *Ω*: longitude of ascending node (radians, in [0, 2π])
            - *ω*: argument of periapsis (radians, in [0, 2π])
            - *f*: true anomaly (radians, in [0, 2π])

        *retrogade* (:class:`bool`): Whether to use the retrograde equinoctial frame.

    Returns:
        :class:`list` of :class:`float`:
            A list of six equinoctial orbital elements:

            - *p*: semi-latus rectum (in units L)  
            - *f*: eccentricity vector times cos(Ω+ω)
            - *g*: eccentricity vector times sin(Ω+ω)
            - *h*: tan(i/2) cos Ω
            - *k*: tan(i/2) sin Ω
            - *L*: true longitude (radians, in [0, 2π])

    Examples:
      >>> par = [7000e3, 0.01, 0.1, 1.0, 0.5, 0.3]
      >>> retro = False
      >>> eq = par2mee(par, retro)
      >>> print("Equinoctial elements:", eq)
)";
}


std::string bielliptic_doc()
{
    return R"(bielliptic(r1, r2, rb, mu)

    Computes the delta v and transfer time required for a bielliptic transfer between two circular orbits.

    Args:
        *r1* (:class:`float`): radius of the first orbit

        *r2* (:class:`float`): radius of the second orbit

        *rb* (:class:`float`): radius of the intermediate orbit

        *mu* (:class:`float`): gravitational parameter of the central body

    Returns:
        [:class:`float`, :class:`float`, [:class:`float`, :class:`float`, :class:`float`]]:
            [total delta v, transfer time, [first delta v, second delta v, third delta v]]

    Examples:
      >>> import pykep as pk
      >>> r1 = 7000000
      >>> r2 = 9000000
      >>> rb = 11000000
      >>> mu = pk.MU_EARTH
      >>> dv_total, t_transfer, [dv1, dv2, dv3] = pk.bielliptic(r1, r2, rb, mu)
      >>> print("Total delta v:", dv_total, "m/s")
      >>> print("Transfer time:", t_transfer, "s")
      >>> print("First delta v:", dv1, "m/s")
      >>> print("Second delta v:", dv2, "m/s")
      >>> print("Third delta v:", dv3, "m/s")
)";
}

std::string mima_doc()
{
    return R"(mima(dv1, dv2, tof, Tmax, veff)
    
    The Maximum Initial Mass Approximation.

    Izzo, D., Hennes, D., Simões, L. F., & Märtens, M. (2016). Designing complex interplanetary trajectories
    for the global trajectory optimization competitions. Space Engineering: Modeling and Optimization with Case Studies, 151-176.
    https://www.esa.int/gsp/ACT/doc/MAD/pub/ACT-RPR-MAD-2016-SSCI-Fast_approximators.pdf

    Having computed a two-impulse transfer, this approximation allows to compute
    the maximum initial mass that a spacecraft can have as to be able to perform
    that transfer in low-thrust.

    Args:
        *dv1* (:class:`list`): First  vectorial delta v (m/s, or Any velocity units)

        *dv2* (:class:`list`): Second vectorial delta v (m/s or Any velocity units)

        *tof* (:class:`float`): Time of flight (Any time units)

        *Tmax* (:class:`float`, optional): Maximum spacecraft thrust.

        *veff* (:class:`float`, optional): Isp*G0.

    Returns:
        :class:`float`, :class:`float`: mima and magnitude of the acceleration required 
        (units induced by the inputs)

    Examples:
      >>> import numpy as np
      >>> import pykep as pk
      >>> dv1 = np.array([320,-345,43]) #m/s
      >>> dv2 = np.array([-510,175,87]) #m/s
      >>> tof = 150*24*60*60 #seconds
      >>> mima, a_required  = pk.mima(dv1, dv2, tof, Tmax = 0.6, veff=3000*pk.G0)
      >>> print("Maximum initial mass:", mima, "kg")
      >>> print("Required acceleration:", a_required*1000, "mm/s)
)";
}

std::string mima_from_hop_doc()
{
    return R"(mima_from_hop(pl_s, pl_t, when_s, when_t, Tmax, veff, mu)

    The Maximum Initial Mass Approximation from a hop transfer. A hop transfer is fixed time rendezvous between two planets.

    Izzo, D., Hennes, D., Simões, L. F., & Märtens, M. (2016). Designing complex interplanetary trajectories
    for the global trajectory optimization competitions. Space Engineering: Modeling and Optimization with Case Studies, 151-176.
    https://www.esa.int/gsp/ACT/doc/MAD/pub/ACT-RPR-MAD-2016-SSCI-Fast_approximators.pdf

    Args:
        *pl_s* (:class:`~pykep.planet`): the source planet

        *pl_t* (:class:`~pykep.planet`): the target planet

        *when_s* (:class:`~pykep.epoch`): the departure epoch

        *when_t* (:class:`~pykep.epoch`): the arrival epoch

        *Tmax* (:class:`float`, optional): Maximum spacecraft thrust.

        *veff* (:class:`float`, optional): Isp*G0.

        *mu* (:class:`float`, optional): gravitational parameter of the central body.

    Returns:
        :class:`float`, :class:`float`: mima and magnitude of the acceleration required

    Examples:
      >>> import pykep as pk
      >>> ... # assuming to have two planets pl_s and pl_t
      >>> when_s = pk.epoch(0.0)
      >>> when_t = pk.epoch(100.0)
      >>> Tmax = 0.6
      >>> veff = 3000*pk.G0
      >>> mu = pk.MU_SUN
      >>> mima, a_required = pk.mima_from_hop(pl_s, pl_t, when_s, when_t, Tmax, veff, mu)
      >>> print("Maximum initial mass:", mima, "kg")
      >>> print("Required acceleration:", a_required*1000, "mm/s)
)";
}

std::string mima2_doc()
{
    return R"(mima2(posvel1, dv1, dv2, tof, Tmax, veff, mu)
    
    More accurate approximation of mima.

    Having computed a two-impulse transfer, this approximation allows to compute
    the maximum initial mass that a spacecraft can have as to be able to perform
    that transfer in low-thrust.

    Izzo, D., ... & Yam, C. H. (2025). Asteroid mining: ACT&Friends’ results for the GTOC12 problem. Astrodynamics, 9(1),
    19-40.  (https://arxiv.org/pdf/2410.20839)

    Args:
        *posvel1* (:class:`list` [:class:`list`, :class:`list`]): initial position and velocty ALONG THE LAMBERT TRANSFER.

        *dv1* (:class:`numpy.ndarray`): First delta v (m/s, or Any velocity units)

        *dv2* (:class:`numpy.ndarray`): Second delta v (m/s or Any velocity units)

        *tof* (:class:`float`): Time of flight (Any time units)

        *Tmax* (:class:`float`): Maximum spacecraft thrust.

        *veff* (:class:`float`): Isp*G0.

        *mu* (:class:`float`): gravitational parameter of the central body.

    Returns:
        :class:`float`, :class:`float`: mima and magnitude of the acceleration required 
        (units induced by the inputs)
)";
}

std::string mima2_from_hop_doc()
{
    return R"(mima2_from_hop(pl_s, pl_t, when_s, when_t, Tmax, veff, mu)

    More accurate approximation of mima from a hop transfer. A hop transfer is fixed time transfer between two planets.

    Izzo, D., ... & Yam, C. H. (2025). Asteroid mining: ACT&Friends’ results for the GTOC12 problem. Astrodynamics, 9(1),
    19-40.  (https://arxiv.org/pdf/2410.20839)

    Args:
        *pl_s* (:class:`~pykep.planet`): the source planet

        *pl_t* (:class:`~pykep.planet`): the target planet

        *when_s* (:class:`~pykep.epoch`): the departure epoch

        *when_t* (:class:`~pykep.epoch`): the arrival epoch

        *Tmax* (:class:`float`, optional): Maximum spacecraft thrust.

        *veff* (:class:`float`, optional): Isp*G0.

        *mu* (:class:`float`, optional): gravitational parameter of the central body.

    Returns:
        :class:`float`, :class:`float`: mima and magnitude of the acceleration required

    Examples:
      >>> import pykep as pk
      >>> ... # assuming to have two planets pl_s and pl_t
      >>> when_s = pk.epoch(0.0)
      >>> when_t = pk.epoch(100.0)
      >>> Tmax = 0.6
      >>> veff = 3000*pk.G0
      >>> mu = pk.MU_SUN
      >>> mima2, a_required = pk.mima2_from_hop(pl_s, pl_t, when_s, when_t, Tmax, veff, mu)
      >>> print("Maximum initial mass:", mima, "kg")
      >>> print("Required acceleration:", a_required*1000, "mm/s)
)";
}

std::string alpha2direct_doc()
{
    return R"(alpha2direct(alphas)

    Converts from alpha encoded to transfer times.

    Args:
        *alphas* (:class:`list`): a sequence of transfer times encoded using the alpha encoding.

        *T* (:class:`float`): the total transfer time.

    Returns:
        :class:`list`:: The encoded transfer times
)";
}

std::string direct2alpha_doc()
{
    return R"(direct2alpha(tofs)

    Converts from transfer times to alpha encoded.

    Args:
        *tofs* (:class:`list`): a sequence of transfer times.

    Returns:
        :class:`list`:, :class:`float`: The alpha-encoded transfer times, the total transfer time (for cenvenience)
)";
}

std::string eta2direct_doc()
{
    return R"(eta2direct(etas, max_tof)

    Converts from eta encoded to transfer times.

    Args:
        *etas* (:class:`list`): a sequence of transfer times encoded using the eta encoding.

        *max_tof* (:class:`float`): maximum time of flight (might actually be less according to the value of the last eta.)

    Returns:
        :class:`list`: The encoded transfer times
)";
}

std::string direct2eta_doc()
{
    return R"(direct2eta(tofs, max_tof)

      Converts from transfer times to eta encoded.

    Args:
        *tofs* (:class:`list`):  a sequence of transfer times

        *max_tof* (:class:`float`): maximum time of flight (might actually be less according to the value of the last eta.)

    Returns:
        :class:`list`: The eta-encoded transfer times
)";
}

std::string epoch_from_float_doc()
{
    return R"(__init__(when: float, julian_type = MJD2000)
    
    Constructs an epoch from a Julian Date.

    Args:
      *when* (:class:`float`): the Julian Date (days since reference)

      *julian_type* (:class:`~pk.epoch.julian_type`): one of MJD2000, JD or MJD

    Examples:
      >>> import pykep as pk
      >>> pk.epoch(12.3, pk.epoch.julian_type.MJD2000)
      2000-01-13T07:12:00.000000
)";
}

std::string epoch_from_datetime_doc()
{
    return R"(**Alternative Constructor:**
    **__init__(** *when: datetime.datetime* **)**
    
    Constructs an epoch from a datetime object.

    Args:
      *when* (:class:`datetime.datetime`): a date

    Examples:
      >>> import pykep as pk
      >>> from datetime import datetime
      >>> pk.epoch(datetime(year=2000, month=1, day=13))
      2000-01-13T00:00:00.000000
)";
}

std::string epoch_from_string_doc()
{
    return R"(**Alternative Constructor:**
    **__init__(** *when: str*, *string_format = pk.epoch.string_format.ISO* **)**
    
    Constructs an epoch from a string.

    Args:
      *when* (:class:`str`): a date

      *string_format* (:class:`~pykep.epoch.string_format`): string format.

    Examples:
      >>> import pykep as pk
      >>> pk.epoch("2000-01-14T00:00:00.000001")
      2000-01-14T00:00:00.000001
)";
}

std::string planet_docstring()
{
    return R"(__init__(udpla)

Planet class.

This type-erasing class represents a generic object moving in space. 
Basically anything which can be defined by its position and velocity in some reference frame.

In order to define a planet in pykep, the user must first define a class
whose methods describe the properties of the planet and allow to compute
its ephemerides (position and velocity), possibly its osculating elements, etc.. 
In pykep, we refer to such a class as a **user-defined planet**, or UDPLA for short. 
Once defined and instantiated, a UDPLA can then be used to construct an instance
of this class, the :class:`~pykep.planet`.

Every UDPLA must implement at least the following method:

.. code-block::

   def eph(self, epoch):
     ...

The ``eph()`` method is expected to return the Cartesian position and velocity at epoch
in some chosen reference frame.

The ``eph()`` method of the UDPLA will then be accessible from the corresponding
:func:`pykep.planet.eph()` method (see its documentation for information on how the method should be implemented
in the UDPLA and other details).

The mandatory method above allow to define a simple planet, which, in a minimal case,
could actually be also just a fixed point in space, for example if its ``eph()`` method
returns a constant position and zero velocity. 

In order to consider more complex cases, the UDPLA may implement one or more of the following methods:

.. code-block::

   def eph_v(self, mjd2000s):
     ...
   def get_mu_central_body(self):
     ...
   def get_mu_self(self):
     ...
   def get_radius(self):
     ...
   def get_safe_radius(self):
     ...
   def period(self, mjd2000):
     ...
   def elements(self, mjd2000, elements_type):
     ...
   def get_name(self):
     ...
   def get_extra_info(self):
     ...

See the documentation of the corresponding methods in this class for details on how the optional
methods in the UDPLA should be implemented and on how they are used by :class:`~pykep.planet`.

Args:
    *udpla*: a user-defined planet, either C++ or Python

Raises:
    *NotImplementedError*: if *udpla* does not implement the mandatory methods detailed above
    unspecified: any exception thrown by methods of the UDP invoked during construction,
    the deep copy of the UDP, the constructor of the underlying C++ class,
    failures at the intersection between C++ and Python (e.g., type conversion errors, mismatched function
    signatures, etc.)

)";
}

std::string planet_get_extra_info_docstring()
{
    return R"(get_extra_info()

Planet's extra info.

If the UDPLA provides a ``get_extra_info()`` method, then this method will return the output of its ``get_extra_info()``
method. Otherwise, an empty string will be returned. 

The string representation of a :class:`~pykep.planet` contains the output of a call to this method.

Returns:
  :class:`str`: extra info about the UDPLA

Raises:
  unspecified: any exception thrown by the ``get_extra_info()`` method of the UDPLA

)";
}

std::string planet_get_name_docstring()
{
    return R"(get_name()

Planet's name.

If the UDPLA provides a ``get_name()`` method, then this method will return the output of its ``get_name()`` method.
Otherwise, an implementation-defined name based on the type of the UDPLA will be returned.

The string representation of a :class:`~pykep.planet` contains the output of a call to this method.

Returns:
    :class:`str`: the problem's name

)";
}

std::string planet_get_mu_central_body_docstring()
{
    return R"(get_mu_central_body()

The gravitational parameter in SI units (m^3/sec^2) of a main body of attraction.

If the UDPLA provides a ``get_mu_central_body()`` method, then this method will return the output of its ``get_mu_central_body()`` method.
Otherwise, -1 will be returned.


Returns:
    :class:`float`: the central body gravitational parameter.

)";
}

std::string planet_get_mu_self_docstring()
{
    return R"(get_mu_self()

The gravitational parameter in SI units (m^3/sec^2) of the planet.

If the UDPLA provides a ``get_mu_self()`` method, then this method will return the output of its ``get_mu_self()`` method.
Otherwise, -1 will be returned.


Returns:
    :class:`float`: the planet's gravitational parameter.

)";
}

std::string planet_get_radius_docstring()
{
    return R"(get_radius()

An average radius in SI units (m) of the planet.

If the UDPLA provides a ``get_radius()`` method, then this method will return the output of its ``get_radius()`` method.
Otherwise, -1 will be returned.


Returns:
    :class:`float`: the planet's average radius.

)";
}

std::string planet_get_safe_radius_docstring()
{
    return R"(get_safe_radius()

The safe radius in SI units (m) of the planet. This is mainly for use in planetary fly-by manoeuvres as to avoid
the planet atmosphere or circumvent its radiation environment.

If the UDPLA provides a ``get_safe_radius()`` method, then this method will return the output of its ``get_safe_radius()`` method.
Otherwise, -1 will be returned.


Returns:
    :class:`float`: the planet's safe radius.

)";
}

std::string planet_eph_docstring()
{
    return R"(eph(when = 0.)

The planet ephemerides, i.e. its position and velocity.

In order to be able to construct a :class:`~pykep.planet` object, the user must provide his own UDPLA (User-Defined-Planet).
This is a class that must implement the method: 

.. code-block::

   def eph(self, mjd2000: float):
      ...
      return [[float, float, float], [float, float, float]]

.. note::
   In the udpla, the signature for eph demands a float as epoch (mjd2000). The planet, instead, constructed from the same udpla, will also allow :class:`~pykep.epoch`.
   

Args:
    *when* (:class:`float` or :class:`~pykep.epoch`): the epoch at which compute the period. When a :class:`float` is passed mjd2000 is assumed.

Returns:
    :class:`list` [:class:`list`, :class:`list`]: r and v, that is the final position and velocity after the propagation.

)";
}

std::string planet_eph_v_docstring()
{
    return R"(eph_v(mjd2000s)

The planet ephemerides, i.e. position and velocity (vectorized version over many epochs).

This method is the vectorized version of its companion :func:`~pykep.planet.eph` and, in its default implementation, it just
calls it in a loop. This behaviour can be changed by the user (for efficiency purposes) who can provide a more efficient version
in his UDPLA by coding a method having the signature: 


.. code-block::

   def eph_v(self, mjd2000s):
      ...
      return np.array((len(mjd2000s), 6))

see, for example, the python implementation of the UDPLAS :class:`~pykep.udpla.tle` and :class:`~pykep.udpla.spice`.

Args:
    *mjd2000s* (:class:`numpy.ndarray` or :class:`list`): the Modified Julian Dates at which to compute the ephemerides.

Returns:
    :class:`list` [:class:`list`, :class:`list`]: r and v, that is the final position and velocity after the propagation.

)";
}

std::string planet_acc_docstring()
{
    return R"(acc(when = 0.)


The planet acceleration, i.e. its inertial acceleration vector at a given epoch.


In order to be able to construct a :class:`~pykep.planet` object, the user must provide his own UDPLA (User-Defined-Planet).
This is a class that must implement the method:

.. code-block::

   def acc(self, mjd2000: float):
      ...
      return [float, float, float]

.. note::
   In the udpla, the signature for acc demands a float as epoch (mjd2000). The planet, instead, constructed from the same udpla, will also allow :class:`~pykep.epoch`.


Args:
    *when* (:class:`float` or :class:`~pykep.epoch`): the epoch at which compute the acceleration. When a :class:`float` is passed mjd2000 is assumed.

Returns:
    :class:`list` [:class:`float`, :class:`float`, :class:`float`]: the acceleration vector at the requested epoch.

)";
}

std::string planet_acc_v_docstring()
{
    return R"(acc_v(mjd2000s)

The planet acceleration (vectorized version over many epochs).

This method is the vectorized version of its companion :func:`~pykep.planet.acc` and, in its default implementation, it just
calls it in a loop. This behaviour can be changed by the user (for efficiency purposes) who can provide a more efficient version
in his UDPLA by coding a method having the signature:

.. code-block::

   def acc_v(self, mjd2000s):
      ...
      return np.array((len(mjd2000s), 3))

Args:
    *mjd2000s* (:class:`numpy.ndarray` or :class:`list`): the Modified Julian Dates at which to compute the accelerations.

Returns:
    :class:`list` [:class:`list`]: a, that is the acceleration vectors at the requested epochs.
)";
}

std::string planet_period_docstring()
{
    return R"(period(when = 0.)

The period of the planet in seconds.

If the UDPLA provides a ``period(float)`` method, ``planet.period`` will call it.
Otherwise, if the UDPLA provides a ``get_mu_self()`` method ``planet.period`` will return the period as computed by the
equation:

.. math::
   T = 2 \pi \sqrt{\frac{a^3}{\mu}}

Else, -1 will be returned.

Args:
    *when* (:class:`float` or :class:`~pykep.epoch`): the epoch at which compute the period. When a :class:`float` is passed mjd2000 is assumed.

Returns:
    :class:`float`: the planet's period.

)";
}

std::string planet_elements_docstring()
{
    return R"(elements(when = 0., el_type = KEP_F)

The elements of the planet at epoch.

If the UDPLA provides a ``elements(float, pk.el_type)`` method, then ``planet.elements`` will call it.
Otherwise, if the UDPLA provides a ``get_mu_self()`` method ``planet.elements`` will return the elements as computed by the
:func:`pykep.ic2par`. Otherwise, -1 will be returned.

Args:
    *when* (:class:`float` or :class:`~pykep.epoch`): the epoch at which compute the period. When a :class:`float` is passed mjd2000 is assumed.

    *el_type* (:class:`~pykep.el_type`): the elements type.

Returns:
    :class:`list`: the planet's elements at epoch.

)";
}

std::string udpla_keplerian_from_posvel_docstring()
{
    return R"(**Alternative Constructor:**
    __init__(ep, posvel, mu_central_body, name = "unknown", added_params = [-1,-1,-1])

Constructs a Keplerian udpla from its position and velocity at epoch.

Args:
    *ep* (:class:`~pykep.epoch`): the epoch at which the orbital elements are provided.

    *posvel* (:class:`list` [:class:`list`, :class:`list`]): the body position and velocty.

    *mu_central_body* (:class:`float`): the gravitational parameter of the main attracting body.

    *name* (:class:`str`): the name of the orbiting body.

    *added_params* (:class:`list`): the body gravitational parameter, its radius and its safe radius. (if -1 they are assumed unknown)

Examples:
    >>> import pykep as pk
    >>> r = [1, 0, 0]
    >>> v = [0, 1, 0]
    >>> ep = pk.epoch("2025-03-22")
    >>> udpla = pk.udpla.keplerian(ep = ep, posvel = [r, v], mu_central_body =1, name = "my_pla")
    >>> pla = pk.planet(udpla)
)";
}

std::string udpla_keplerian_from_elem_docstring()
{
    return R"(__init__(when, elem, mu_central_body, name = "unknown", added_params = [-1,-1,-1], elem_type = KEP_F)

Constructs a Keplerian udpla from its orbital elements at epoch.

Args:
    *when* (:class:`~pykep.epoch`): the epoch at which the orbital elements are provided.

    *elem* (:class:`list` or :class:`numpy.ndarray`): the orbital elements. by default.

    *mu_central_body* (:class:`float`): the gravitational parameter of the main attracting body.

    *name* (:class:`str`): the name of the orbiting body.

    *added_params* (:class:`list`): the body gravitational parameter, its radius and its safe radius. (if -1 they are assumed unknown)

    *el_type* (:class:`~pykep.el_type`): the elements type. Defaults to osculating Keplerian (a ,e ,i, W, w, f) with true anomaly.

Examples:
    >>> import pykep as pk
    >>> elem = [1, 0, 0, 0, 0, 0]
    >>> when = pk.epoch("2025-03-22")
    >>> udpla = pk.udpla.keplerian(when = when, elem = elem, mu_central_body =1, name = "my_pla")
    >>> pla = pk.planet(udpla)
)";
}

std::string udpla_jpl_lp_docstring()
{
    return R"(__init__(name = "earth")

Constructs a solar system planet with ephemerides computed using a low-precision (non Keplerian)
model from JPL (https://ssd.jpl.nasa.gov/planets/approx_pos.html).

Args:
    *name* (:class:`str`): the name of the solar system planet.

Examples:
    >>> import pykep as pk
    >>> udpla = pk.udpla.jpl_lp(name="mercury")
    >>> pla = pk.planet(udpla)
)";
}

std::string udpla_vsop2013_docstring()
{
    return R"(__init__(body = "mercury", thresh = 1e-5)

Constructs a solar system planet with ephemerides computed using the VSOP2013 analytical
theory (https://en.wikipedia.org/wiki/VSOP_model).

Args:
    *body* (:class:`str`): the name of the solar system planet.

    *thresh* (:class:`float`): the truncation threshold for the theory's coefficients.

Examples:
    >>> import pykep as pk
    >>> udpla = pk.udpla.vsop2013(body="venus")
    >>> pla = pk.planet(udpla)
)";
}

std::string lambert_problem_docstring()
{
    return R"(__init__(r0 = [1,0,0], r1 = [0,1,0], tof = pi/2, mu = 1., cw = False, multi_revs = 0)

      Args:
          *r0* (1D array-like): Cartesian components of the first position vector [xs, ys, zs]. Defaults to [1,0,0].

          *r1* (1D array-like): Cartesian components of the second position vector [xf, yf, zf]. Defaults tot [0,1,0].

          *tof* (:class:`float`): time of flight. Defaults to :math:`\frac{\pi}{2}`.

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *cw* (:class:`bool`): True for retrograde motion (clockwise). Defaults to False.

          *multi_revs* (:class:`float`): Maximum number of multiple revolutions to be computed. Defaults to 0.

      .. note::

        Units need to be consistent. The multirev Lambert's problem will be solved upon construction
        and its solution stored in data members.

      Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> r0 = [1,0,0]
        >>> r1 = [0,1,0]
        >>> tof = np.pi/2
        >>> mu = 1.
        >>> lp = pk.lambert_problem(r0, r1, tof, mu)
        >>> lp.v0[0]
        [-4.1028493158958256e-16, 1.0000000000000002, 0.0]
)";
}

std::string get_kep_docstring()
{
    return R"(ta.get_kep(tol)

Returns a Taylor adaptive propagator (Heyoka) for the simple Keplerian dynamics retrieving one from a global cache and making a copy. 

If the requested propagator was never created, a call to this function will trigger its creation, else it will
return the one from a global cache, thus avoiding jitting.

In `pykep`, the simple keplerian dynamics is defined in Cartesian coordinates (thus it is not symplectic as not in a Hamiltonian form). 

The specific dynamics used is that returned by :func:`~pykep.ta.kep_dyn`.

Args:
    *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

Returns:
    :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> ta = pk.ta.get_kep(tol = 1e-16)
  >>> ta.time = 0.
  >>> ta.state[:] = [1.01238082345234, -0.0423523523454,  0.22634376321, -0.1232623614,    0.123462698209365, 0.123667064622]
  >>> mu = 0.01215058560962404
  >>> tof = 5.7856656782589234
  >>> ta.pars[0] = mu
  >>> ta.propagate_until(tof)
)";
}

std::string get_kep_var_docstring()
{
    return R"(ta.get_kep_var(tol)

Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the simple Keplerian dynamics retrieving one from a global cache and making a copy. 

If the requested propagator was never created, a call to this function will trigger its creation, else it will
return the one from a global cache, thus avoiding jitting.

.. note::
   Variations are only considered with respect to initial conditions.

In `pykep`, the simple Keplerian dynamics is defined in Cartesian coordinates (thus it is not symplectic as not in a Hamiltonian form). 

The specific dynamics used is that returned by :func:`~pykep.ta.kep_dyn`.

Args:
    *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

Returns:
    :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> ta = pk.ta.get_kep_var(tol = 1e-16)
  >>> ta.time = 0.
  >>> ta.state[:] = [1.01238082345234, -0.0423523523454,  0.22634376321, -0.1232623614,    0.123462698209365, 0.123667064622]
  >>> mu = 0.01215058560962404
  >>> tof = 5.7856656782589234
  >>> ta.pars[0] = mu
  >>> ta.propagate_until(tof)
)";
}

std::string kep_dyn_docstring()
{
    return R"(kep_dyn()

The dynamics of the simple Keplerian problem (kep).

In `pykep`, the kep is defined in Cartesian coordinates (thus it is not symplectic as not in a Hamiltonian form). 

The parameter :math:`\mu` is the central body gravitational parameter. 


.. math::
   \left\{
   \begin{array}{l}
       \dot{\mathbf r} = \mathbf v \\
       \dot{\mathbf v} = -\frac{mu}{r^3} \mathbf r \\
   \end{array}\right.

where :math:`\mu` is the only parameter.

Returns:
    :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(x, dx), ...]
)";
}

std::string get_zoh_kep_docstring()
{
  return R"(ta.get_zoh_kep(tol)

Returns a Taylor adaptive propagator (Heyoka) for the zero-order-hold Keplerian low-thrust dynamics,
retrieving one from a global cache and making a copy.

If the requested propagator was never created, this creates it; otherwise returns a copy of the cached
version (avoiding re-jitting while keeping independent runtime state).

The specific dynamics used is that returned by :func:`~pykep.ta.zoh_kep_dyn`.

Args:
  *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator.

Returns:
  :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> veff = 1.32
  >>> controls = [0.022, 0.023, -0.21, 0.1]
  >>> ta = pk.ta.get_zoh_kep(tol = 1e-16)
  >>> ta.time = 0.
  >>> ta.state[:] = [1.2, 0.1, 0.1, 0.1, 1., 0.123, 1.]
  >>> ta.pars[:] = controls + [1. / veff]
  >>> ta.propagate_until(1.23)
)";
}

std::string get_zoh_kep_var_docstring()
{
  return R"(ta.get_zoh_kep_var(tol)

Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the zero-order-hold
Keplerian low-thrust dynamics, retrieving one from a global cache and making a copy.

If the requested propagator was never created, this creates it; otherwise returns a copy of the cached
version (avoiding re-jitting while keeping independent runtime state).

.. note::
   Variations are only considered with respect to initial conditions and the thrust parameters
   :math:`[T, i_x, i_y, i_z]`.

The specific dynamics used is that returned by :func:`~pykep.ta.zoh_kep_dyn`.

Args:
  *tol* (:class:`float`): the tolerance of the variational Taylor adaptive propagator.

Returns:
  :class:`hy::taylor_adaptive`: The variational Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> veff = 1.32
  >>> controls = [0.022, 0.023, -0.21, 0.1]
  >>> ta = pk.ta.get_zoh_kep_var(tol = 1e-16)
  >>> ta.time = 0.
  >>> ta.state[:7] = [1.2, 0.1, 0.1, 0.1, 1., 0.123, 1.]
  >>> ta.pars[:] = controls + [1. / veff]
  >>> ta.propagate_until(1.23)
)";
}

std::string zoh_kep_dyn_docstring()
{
  return R"(zoh_kep_dyn()

The dynamics in Cartesian coordinates of a constant-thrust mass-varying spacecraft
orbiting a central body with unitary gravitational parameter :math:`\mu = 1`.

The electric propulsion is characterized by its effective exhaust velocity
:math:`v_{eff} = I_{sp} g_0`, with :math:`c = 1 / v_{eff}`.

The state is: :math:`[x, y, z, v_x, v_y, v_z, m]`

The system parameters are (in this order): :math:`[T, i_x, i_y, i_z] + [c]`.

Returns:
  :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]:
  The dynamics in the form [(x, dx), ...]
)";
}

std::string get_zoh_eq_docstring()
{
  return R"(ta.get_zoh_eq(tol)

Returns a Taylor adaptive propagator (Heyoka) for the zero-order-hold equinoctial low-thrust dynamics,
retrieving one from a global cache and making a copy.

If the requested propagator was never created, this creates it; otherwise returns a copy of the cached
version (avoiding re-jitting while keeping independent runtime state).

The specific dynamics used is that returned by :func:`~pykep.ta.zoh_eq_dyn`.

Args:
  *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator.

Returns:
  :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.
)";
}

std::string get_zoh_eq_var_docstring()
{
  return R"(ta.get_zoh_eq_var(tol)

Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the zero-order-hold
Equinoctial low-thrust dynamics, retrieving one from a global cache and making a copy.

If the requested propagator was never created, this creates it; otherwise returns a copy of the cached
version (avoiding re-jitting while keeping independent runtime state).

.. note::
   Variations are only considered with respect to initial conditions and the thrust parameters
   :math:`[T, i_r, i_t, i_n]`.

The specific dynamics used is that returned by :func:`~pykep.ta.zoh_eq_dyn`.

Args:
  *tol* (:class:`float`): the tolerance of the variational Taylor adaptive propagator.

Returns:
  :class:`hy::taylor_adaptive`: The variational Taylor adaptive propagator.
)";
}

std::string zoh_eq_dyn_docstring()
{
  return R"(zoh_eq_dyn()

The dynamics in equinoctial elements of a constant-thrust mass-varying spacecraft
orbiting a central body with unitary gravitational parameter :math:`\mu = 1`.

The electric propulsion is characterized by its effective exhaust velocity
:math:`v_{eff} = I_{sp} g_0`, with :math:`c = 1 / v_{eff}`.

The state is: :math:`[p, f, g, h, k, L, m]`

The system parameters are (in this order): :math:`[T, i_r, i_t, i_n] + [c]`.

Returns:
  :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]:
  The dynamics in the form [(x, dx), ...]
)";
}

std::string get_zoh_cr3bp_docstring()
{
  return R"(ta.get_zoh_cr3bp(tol)

Returns a Taylor adaptive propagator (Heyoka) for the zero-order-hold CR3BP low-thrust dynamics,
retrieving one from a global cache and making a copy.

If the requested propagator was never created, this creates it; otherwise returns a copy of the cached
version (avoiding re-jitting while keeping independent runtime state).

The specific dynamics used is that returned by :func:`~pykep.ta.zoh_cr3bp_dyn`.

Args:
  *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator.

Returns:
  :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.
)";
}

std::string get_zoh_cr3bp_var_docstring()
{
  return R"(ta.get_zoh_cr3bp_var(tol)

Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the zero-order-hold
CR3BP low-thrust dynamics, retrieving one from a global cache and making a copy.

If the requested propagator was never created, this creates it; otherwise returns a copy of the cached
version (avoiding re-jitting while keeping independent runtime state).

.. note::
   Variations are only considered with respect to initial conditions and the thrust parameters
   :math:`[T, i_x, i_y, i_z]`.

The specific dynamics used is that returned by :func:`~pykep.ta.zoh_cr3bp_dyn`.

Args:
  *tol* (:class:`float`): the tolerance of the variational Taylor adaptive propagator.

Returns:
  :class:`hy::taylor_adaptive`: The variational Taylor adaptive propagator.
)";
}

std::string zoh_cr3bp_dyn_docstring()
{
  return R"(zoh_cr3bp_dyn()

The dynamics of a fixed-thrust mass-varying spacecraft in the circular restricted
three-body problem (CR3BP), with thrust fixed in the rotating frame.

The state is: :math:`[x, y, z, v_x, v_y, v_z, m]`

The system parameters are (in this order): :math:`[T, i_x, i_y, i_z] + [c, \mu]`.

Returns:
  :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]:
  The dynamics in the form [(x, dx), ...]
)";
}

std::string get_zoh_ss_docstring()
{
  return R"(ta.get_zoh_ss(tol)

Returns a Taylor adaptive propagator (Heyoka) for the zero-order-hold solar-sail dynamics,
retrieving one from a global cache and making a copy.

If the requested propagator was never created, this creates it; otherwise returns a copy of the cached
version (avoiding re-jitting while keeping independent runtime state).

The specific dynamics used is that returned by :func:`~pykep.ta.zoh_ss_dyn`.

Args:
  *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator.

Returns:
  :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.
)";
}

std::string get_zoh_ss_var_docstring()
{
  return R"(ta.get_zoh_ss_var(tol)

Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the zero-order-hold
solar-sail dynamics, retrieving one from a global cache and making a copy.

If the requested propagator was never created, this creates it; otherwise returns a copy of the cached
version (avoiding re-jitting while keeping independent runtime state).

.. note::
   Variations are only considered with respect to initial conditions and the sail-angle parameters
   :math:`[\alpha, \beta]`.

The specific dynamics used is that returned by :func:`~pykep.ta.zoh_ss_dyn`.

Args:
  *tol* (:class:`float`): the tolerance of the variational Taylor adaptive propagator.

Returns:
  :class:`hy::taylor_adaptive`: The variational Taylor adaptive propagator.
)";
}

std::string zoh_ss_dyn_docstring()
{
  return R"(zoh_ss_dyn()

The dynamics in Cartesian coordinates of a solar-sail spacecraft around a central body
with unitary gravitational parameter :math:`\mu = 1`.

The state is: :math:`[x, y, z, v_x, v_y, v_z]`

The system parameters are (in this order): :math:`[\alpha, \beta] + [c]`.

Returns:
  :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]:
  The dynamics in the form [(x, dx), ...]
)";
}

std::string get_bcp_docstring()
{
    return R"(ta.get_bcp(tol)

Returns a Taylor adaptive propagator (Heyoka) for the Bicircular Problem (BCP) dynamics
retrieving one from a global cache (making a copy).

If the requested propagator was never created, this creates it; otherwise returns cached version
(avoiding re-jitting).

The BCP is defined in Cartesian coordinates (non-symplectic, non-Hamiltonian, time-dependent)
with Sun on the :math:`x` axis at :math:`t=0`.

Args:
    *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator.

Returns:
    :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
    Import and setup::

        import pykep as pk
        ta = pk.ta.get_bcp(tol=1e-16)
        ta.time = 0.
        ta.state[:] = [1.01238082345234, -0.0423523523454, 0.22634376321, 
                      -0.1232623614, 0.123462698209365, 0.123667064622]
        mu = pk.CR3BP_MU_EARTH_MOON
        mu_s = pk.BCP_MU_S
        rho_s = pk.BCP_RHO_S
        rho_p = pk.BCP_RHO_S
        tof = 5.7856656782589234

    Create propagator and propagate::

        ta.pars[:] = [mu, mu_s, rho_s, rho_p]
        ta.propagate_until(tof)
)";
}


std::string get_bcp_var_docstring()
{
    return R"(ta.get_bcp_var(tol:float)

Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the 
Bicircular Problem (BCP) dynamics retrieving one from a global cache (making a copy).

If the requested propagator was never created, this creates it; otherwise returns cached version
(avoiding re-jitting).

.. note::
   Variations are only considered with respect to initial conditions.

The BCP is defined in Cartesian coordinates (non-symplectic, non-Hamiltonian, time-dependent)
with Sun on the :math:`x` axis at :math:`t=0`.

Args:
    *tol* (:class:`float`): the tolerance of the variational Taylor adaptive propagator.

Returns:
    :class:`hy::taylor_adaptive`: The variational Taylor adaptive propagator.

Examples:
    Import and setup::

        import pykep as pk
        ta = pk.ta.get_bcp_var(tol=1e-16)
        ta.time = 0.
        ta.state[:6] = [1.01238082345234, -0.0423523523454, 0.22634376321, 
                       -0.1232623614, 0.123462698209365, 0.123667064622]
        mu = pk.CR3BP_MU_EARTH_MOON
        mu_s = pk.BCP_MU_S
        rho_s = pk.BCP_RHO_S
        rho_p = pk.BCP_RHO_S
        tof = 5.7856656782589234

    Create propagator and propagate::

        ta.pars[:] = [mu, mu_s, rho_s, rho_p]
        ta.propagate_until(tof)
)";
}


std::string bcp_dyn_docstring()
{
    return R"(bcp_dyn()

The dynamics of the Bicircular Problem (BCP).

In `pykep`, the BCP is defined in Cartesian coordinates (it is not symplectic nor in a Hamiltonian form). It is time-dependent and assumes the Sun
is on the :math:`x` axis at time zero.

The equations are non-dimensional with units :math:`L = r_{12}` (distance between the primaries), :math:`M = m_1 + m_2` (total system mass) and
:math:`T = \sqrt{\frac{r_{12}^3}{m_1+m_2}}` (so that the angular velocity of the primaries is :math:`\omega = 1`
and the period of rotation between the primaries is :math:`2\pi`).
The parameter :math:`\mu` is defined as :math:`\frac{m_2}{m_1+m_2}` where :math:`m_2` is the mass of the
secondary body (i.e. placed on the positive x axis).

The equations of motion are:

.. math::
   \left\{
   \begin{array}{l}
       \dot{\mathbf r} = \mathbf v \\
       \dot v_x = 2v_y + x - (1 - \mu) \frac{x + \mu}{r_1^3} - \mu \frac{x + \mu - 1}{r_2^3} - \frac{\mu_s}{r_s^3} (x- \rho \cos(\omega_s t)) - \frac{\mu_s}{\rho_s^2} \cos(\omega_s t)\\
       \dot v_y = -2 v_x + y - (1 - \mu) \frac{y}{r_1^3} - \mu \frac{y}{r_2^3} - \frac{\mu_s}{r_s^3} (y - \rho \sin(\omega_s t)) - \frac{\mu_s}{\rho_s^2} \sin(\omega_s t)\\
       \dot v_z = -(1 - \mu) \frac{z}{r_1^3} - \mu \frac{z}{r_2^3} - \mu_s \frac{z}{r_s^3} 
   \end{array}\right.

where :math:`\mu, \mu_s, \rho_s` and $\omega_s$ are the system parameters. The spacecraft distance from the primaries
is :math:`r_1, r_2` and its distance from the Sun is :math:`r_s`. The Sun orbits at a distance :math:`\rho` with
angular velocity :math:`\omega_s`:.

Returns:
    :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(x, dx), ...]
)";
}

std::string get_cr3bp_docstring()
{
    return R"(ta.get_cr3bp(tol)

Returns a Taylor adaptive propagator (Heyoka) for the CR3BP problem retrieving one from a global cache and making a copy. 

If the requested propagator was never created, a call to this function will trigger its creation, else it will
return the one from a global cache, thus avoiding jitting.

In `pykep`, the CR3BP is defined in Cartesian coordinates (thus it is not symplectic as not in a Hamiltonian form). 

The specific dynamics used is that returned by :func:`~pykep.ta.cr3bp_dyn`.

Args:
    *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

Returns:
    :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> ta = pk.ta.get_cr3bp(tol = 1e-16)
  >>> ta.time = 0.
  >>> ta.state[:] = [1.01238082345234, -0.0423523523454,  0.22634376321, -0.1232623614,    0.123462698209365, 0.123667064622]
  >>> mu = 0.01215058560962404
  >>> tof = 5.7856656782589234
  >>> ta.pars[0] = mu
  >>> ta.propagate_until(tof)
)";
}

std::string get_cr3bp_var_docstring()
{
    return R"(ta.get_cr3bp_var(tol)

Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the CR3BP problem retrieving one from a global cache and making a copy. 

If the requested propagator was never created, a call to this function will trigger its creation, else it will
return the one from a global cache, thus avoiding jitting.

.. note::
   Variations are only considered with respect to initial conditions.

In `pykep`, the CR3BP is defined in Cartesian coordinates (thus it is not symplectic as not in a Hamiltonian form). 

The specific dynamics used is that returned by :func:`~pykep.ta.cr3bp_dyn`.

Args:
    *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

Returns:
    :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> ta = pk.ta.get_cr3bp_var(tol = 1e-16)
  >>> ta.time = 0.
  >>> ta.state[:] = [1.01238082345234, -0.0423523523454,  0.22634376321, -0.1232623614,    0.123462698209365, 0.123667064622]
  >>> mu = 0.01215058560962404
  >>> tof = 5.7856656782589234
  >>> ta.pars[0] = mu
  >>> ta.propagate_until(tof)
)";
}

std::string cr3bp_dyn_docstring()
{
    return R"(cr3bp_dyn()

The dynamics of the Circular Restricted Three Body Problem (CR3BP).

In `pykep`, the CR3BP is defined in Cartesian coordinates (thus it is not symplectic as not in a Hamiltonian form). 

The parameter :math:`\mu` is defined as :math:`\frac{m_2}{m_1+m_2}` where :math:`m_2` is the mass of the
secondary body (i.e. placed on the positive x axis). 

The equations are non-dimensional with units :math:`L = r_{12}` (distance between the primaries), :math:`M = m_1 + m_2` (total system mass) and
:math:`T = \sqrt{\frac{r_{12}^3}{m_1+m_2}}` (so that the angular velocity of the primaries is :math:`\omega = 1` 
and the period of rotation between the primaries is :math:`2\pi`).

.. math::
   \left\{
   \begin{array}{l}
       \dot{\mathbf r} = \mathbf v \\
       \dot v_x = 2v_y + x - (1 - \mu) \frac{x + \mu}{r_1^3} - \mu \frac{x + \mu - 1}{r_2^3} \\
       \dot v_y = -2 v_x + y - (1 - \mu) \frac{y}{r_1^3} - \mu \frac{y}{r_2^3} \\
       \dot v_z = -(1 - \mu) \frac{z}{r_1^3} - \mu \frac{z}{r_2^3}
   \end{array}\right.

where :math:`\mu` is the only parameter.

Returns:
    :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(x, dx), ...]
)";
}

std::string cr3bp_jacobi_C_docstring()
{
    return R"(cr3bp_jacobi_C()
Jacobi constant of the CR3BP as a `heyoka` expression.

The Jacobi constant is a conserved quantity of the CR3BP. It is defined as:

.. math::

   C = x^2 + y^2 + 2(1 - \mu) \frac{1}{r_1} + 2\mu \frac{1}{r_2} - v_x^2 - v_y^2 - v_z^2 = 2 U - v^2

where :math:`\mu` is the mass ratio of the primaries, :math:`r_1` and :math:`r_2` are the distances to the primaries.

Returns:
    :class:`hy::expression`: the Jacobi constant of the CR3BP as a function of :math:`x,y,z,v_x,v_y,v_z,\mu`.
)";
}

std::string cr3bp_effective_potential_U_docstring()
{
    return R"(cr3bp_effective_potential_U()
Effective potential of the CR3BP as a `heyoka` expression.
The effective potential is defined as:

.. math::

   U = \frac{1}{2} \left(x^2 + y^2 + 2(1 - \mu) \frac{1}{r_1} + 2\mu \frac{1}{r_2}\right)

where :math:`\mu` is the mass ratio of the primaries, :math:`r_1` and :math:`r_2` are the distances to the primaries.

Returns:
    :class:`hy::expression`: the effective potential of the CR3BP as a function of :math:`x,y,z,v_x,v_y,v_z,\mu`.
)";
}

std::string get_pc_docstring()
{
    return R"(ta.get_pc(tol, optimality)

Returns a Taylor adaptive propagator (Heyoka) for the TPBVP problem resulting from the application of 
Pontryagin Maximum Principle (PMC) to the low-trhust problem (constant maximal thrust) in 
Cartesian coordinates. If the requested propagator was never created, a call to this function will
trigger its compilation. Otherwise, it will return the one from a global cache, thus avoiding jitting.
The specific dynamics used is that returned by :func:`~pykep.ta.pc_dyn`.
Both time optimal and mass optimal systems can be returned by setting the *optimality* parameter.

Args:
    *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

    *optimality* (:class:`pykep.optimality_type`): the optimality principle to be used.

Returns:
    :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> ta = pk.ta.get_pc(tol = 1e-16, optimality = pk.optimality_type.TIME)
  >>> ta.time = 0.
  >>> # We set the initial conditions with some arbitrary values (all costates to 1.)
  >>> ta.state[:14] = [1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.]
  >>> ta.pars[:] = [1., 0.01, 1.] # in case of TIME parameters are [mu, Tmax, Isp g0]
  >>> tof = 1.2345
  >>> ta.propagate_until(tof)
)";
}

std::string get_pc_var_docstring()
{
    return R"(ta.get_pc_var(tol, optimality)

Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the TPBVP problem resulting
from the application of Pontryagin Maximum Principle (PMP) to the low-thrust problem
(constant maximal thrust) in Cartesian coordinates. If the requested propagator was never created, 
a call to this function will trigger its compilation. Otherwise, it will return the one from
a global cache, thus avoiding jitting.

.. note::
   Variations are considered with respect to the initial conditions on the costates and to the
   parameters :math:`\epsilon` and :math:`\lambda_0`. In the time optimal case :math:`\epsilon`
   is still considered a parameter (for consistency) but it is not used.

The specific dynamics used is that returned by :func:`~pykep.ta.pc_dyn`.

Args:
    *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

    *optimality* (:class:`pykep.optimality_type`): the optimality principle to be used.

Returns:
    :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> ta_var = pk.ta.get_pc_var(tol = 1e-16, optimality = pk.optimality_type.TIME))
  >>> ta_var.time = 0.
  >>> ta_var.state[:14] = [1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.]
  >>> ta_var.pars[:] = [1., 0.01, 1.] # in case of TIME parameters are [mu, Tmax, Isp g0]
  >>> tof = 1.2345
  >>> ta_var.propagate_until(tof)
)";
}

std::string pc_dyn_docstring()
{
    return R"(ta.pc_dyn(optimality)

The augmented dynamics of the TPBVP originating when applying an indirect method
to the low-thrust transfer problem in Cartesian coordinates.

The (non-augmented) dynamics is,

.. math::
   \left\{
   \begin{array}{l}
        \dot{\mathbf r} = \mathbf f_r = \mathbf v \\
        \dot{\mathbf v} = \mathbf f_v = -\frac{mu}{r^3}\mathbf r + c_1 \frac um \hat{\mathbf i}\\
        \dot{m} = f_m = - \frac{c_1}{c_2} u
   \end{array}
   \right.

The state, containing the spacecraft state and the co-states,  is:

.. math::
   \mathbf x = [\mathbf r, \mathbf v, m] = [x,y,z,v_x,v_y,v_z,l_x,l_y,l_z,l_{vx},l_{vy},l_{vz},l_m]

While the parameters:

.. math::
   \mathbf p = [\mu, c_1, c_2, \epsilon, \lambda_0]

describe, respectively, the gravitational parameter, the maximum thrust, the effective velocity
(product of the specific impulse :math:`I_{sp}` by :math:`g_0`),
an homotopy parameter and a factor multiplying the instantaneous cost (in theory this can be any positive number as 
the solution of the problem will not change)

The controls, representing magnitude and direction of the spacecraft thrust are:

.. math::
   \mathbf u = [u, \hat{\mathbf i}]

The final equations of motion (in the case of mass optimality), are derived from the Hamiltonian:

.. math::
   \mathcal H(\mathbf x, \mathbf \lambda, \mathbf u) = \mathbf \lambda_r \cdot \mathbf f_r + \mathbf \lambda_v \cdot \mathbf f_v + \lambda_m  f_m + \lambda_0 \frac{c_1}{c_2} \left(u + \epsilon\log(u(1-u))\right)

or, in the case of time optimality:

.. math::
   \mathcal H(\mathbf x, \mathbf \lambda, \mathbf u) = \mathbf \lambda_r \cdot \mathbf f_r + \mathbf \lambda_v \cdot \mathbf f_v + \lambda_m  f_m + \lambda_0 \frac{c_1}{c_2}

by taking the derivatives:

.. math::
   \left\{
   \begin{array}{l}
   \dot{\mathbf x} = \frac{\partial \mathcal H}{\partial \mathbf \lambda} \\
   \dot{\mathbf \lambda} = - \frac{\partial \mathcal H}{\partial \mathbf x} \\
   \end{array}\right.

The equation are then made independent of the controls applying Pontryagin minimum principle which takes the general form:

.. math::
   \mathbf u^* = \textbf{argmin}_{\mathbf u \in \mathcal U} \mathcal H(\mathbf x, \mathbf \lambda, \mathbf u)

where :math:`\mathcal U` is the set of admissible controls, i.e. the set of thrust directions and magnitudes.

and substituting in the equations the optimal controls found as function of the augmented system state.

.. note::
   In the case of time optimality, the dynamics will not depend on neither :math:`\epsilon` nor :math:`\lambda_0`
   and the size of the parameter vector will this be three.
   

Args:
    *optimality* (:class:`pykep.optimality_type`): the optimality principle to be used.

Returns:
    :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(x, dx), ...]
)";
}

std::string get_pc_H_cfunc_docstring()
{
    return R"(ta.get_pc_H(optimality)

The Hamiltonian function associated with the two-point boundary value problem (TPBVP)
arising from the application of an indirect method to the low-thrust transfer problem
in Cartesian coordinates.

This returns the compiled function of the Hamiltonian as a sole function of the augmented system state, under the 
assumption of mass-optimality or time-optimality, depending on the provided setting.

The state is composed of both the physical state of the spacecraft and the corresponding costates:

.. math::
   \mathbf x = [\mathbf r, \mathbf v, m] = [x,y,z,v_x,v_y,v_z,m]

.. math::
   \mathbf \lambda = [\lambda_r, \lambda_v, \lambda_m] = [l_x,l_y,l_z,l_{vx},l_{vy},l_{vz},l_m]

The system parameters are:

.. math::
   \mathbf p = [\mu, c_1, c_2, \epsilon, \lambda_0]

In the case of **mass optimality**, the Hamiltonian is given by:

.. math::
   \mathcal H(\mathbf x, \mathbf \lambda, \mathbf u) =
   \lambda_r \cdot \mathbf v +
   \lambda_v \cdot \left(-\frac{\mu}{r^3}\mathbf r + \frac{c_1 u}{m} \hat{\mathbf i} \right) +
   \lambda_m \left(-\frac{c_1}{c_2} u\right) +
   \lambda_0 \frac{c_1}{c_2} \left(u + \epsilon \log(u(1-u))\right)

In the case of **time optimality**, the Hamiltonian becomes:

.. math::
   \mathcal H(\mathbf x, \mathbf \lambda, \mathbf u) =
   \lambda_r \cdot \mathbf v +
   \lambda_v \cdot \left(-\frac{\mu}{r^3}\mathbf r + \frac{c_1 u}{m} \hat{\mathbf i} \right) +
   \lambda_m \left(-\frac{c_1}{c_2} u\right) +
   \lambda_0 \frac{c_1}{c_2}

In both cases the control :math:`u` is assumed to be optimal and is thus a function of the states/co-states. 
Along an optimal trajectory this Hamiltonian must be constant.

Args:
    *optimality* (:class:`pykep.optimality_type`): the optimality principle to be used.

Returns:
    :class:`hy::c_func`: The Hamiltonian :math:`\mathcal H(\mathbf x, \mathbf \lambda)` as a numeric function.

Examples:
    >>> import pykep as pk
    >>> mu = ..
    ...
    >>> H_func([x,y,z,vx,vy,vz,m,lx,ly,lz,lvx,lvy,lvz,lm], pars = [mu, c1, c2, eps, l0])
    )";
}


std::string get_peq_docstring()
{
    return R"(ta.get_peq(tol, optimality)

Returns a Taylor adaptive propagator (Heyoka) for the TPBVP problem resulting from the application of 
Pontryagin Maximum Principle (PMC) to the low-trhust problem (constant maximal thrust) in 
Modified Equinoctial elements. If the requested propagator was never created, a call to this function will
trigger its compilation. Otherwise, it will return the one from a global cache, thus avoiding jitting.
The specific dynamics used is that returned by :func:`~pykep.ta.peq_dyn`.
Both time optimal and mass optimal systems can be returned by setting the *optimality* parameter.

Args:
    *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

    *optimality* (:class:`pykep.optimality_type`): the optimality principle to be used.

Returns:
    :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> ta = pk.ta.get_peq(tol = 1e-16, optimality = pk.optimality_type.TIME)
  >>> ta.time = 0.
  >>> # We set the initial conditions with some arbitrary values (all costates to 1.)
  >>> ta.state[:14] = [1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.]
  >>> ta.pars[:] = [1., 0.01, 1.] # in case of TIME parameters are [mu, Tmax, Isp g0]
  >>> tof = 1.2345
  >>> ta.propagate_until(tof)
)";
}

std::string get_peq_var_docstring()
{
    return R"(ta.get_peq_var(tol, optimality)

Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the TPBVP problem resulting
from the application of Pontryagin Maximum Principle (PMP) to the low-thrust problem
(constant maximal thrust) in Modified Equinoctial Elements. If the requested propagator was never created, 
a call to this function will trigger its compilation. Otherwise, it will return the one from
a global cache, thus avoiding jitting.

.. note::
   Variations are considered with respect to the initial conditions on the costates and to the
   parameters :math:`\epsilon` and :math:`\lambda_0`. In the time optimal case :math:`\epsilon`
   is still considered a parameter (for consistency) but it is not used.

The specific dynamics used is that returned by :func:`~pykep.ta.pc_dyn`.

Args:
    *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

    *optimality* (:class:`pykep.optimality_type`): the optimality principle to be used.

Returns:
    :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

Examples:
  >>> import pykep as pk
  >>> ta_var = pk.ta.get_peq_var(tol = 1e-16, optimality = pk.optimality_type.TIME))
  >>> ta_var.time = 0.
  >>> ta_var.state[:14] = [1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.]
  >>> ta_var.pars[:] = [1., 0.01, 1.] # in case of TIME parameters are [mu, Tmax, Isp g0]
  >>> tof = 1.2345
  >>> ta_var.propagate_until(tof)
)";
}

std::string peq_dyn_docstring()
{
    return R"(ta.peq_dyn(optimality)

The augmented dynamics of the TPBVP originating when applying an indirect method
to the low-thrust transfer problem in **Modified Equinoctial Elements (MEE)**.

The (non-augmented) dynamics is:

.. math::
   \left\{
   \begin{array}{l}
   \dot p = \sqrt{\frac{p}{\mu}} \frac{2p}{w} i_t \cdot \frac{u}{m} \\
   \dot f = \frac{1}{m} \sqrt{\frac{p}{\mu}} \left[ i_r \sin L + \left( (1+w)\cos L + f \right)\frac{i_t}{w} - (h\sin L - k\cos L)\frac{g i_n}{w} \right] \cdot \frac{u}{m} \\
   \dot g = \frac{1}{m} \sqrt{\frac{p}{\mu}} \left[ -i_r \cos L + \left( (1+w)\sin L + g \right)\frac{i_t}{w} + (h\sin L - k\cos L)\frac{f i_n}{w} \right] \cdot \frac{u}{m} \\
   \dot h = \sqrt{\frac{p}{\mu}} \frac{s^2 i_n}{2w} \cos L \cdot \frac{u}{m} \\
   \dot k = \sqrt{\frac{p}{\mu}} \frac{s^2 i_n}{2w} \sin L \cdot \frac{u}{m} \\
   \dot L = \sqrt{\frac{p}{\mu}} \left[ \mu\left(\frac{w}{p}\right)^2 + \frac{1}{w}(h\sin L - k\cos L) \frac{i_n}{m} \cdot \frac{u}{m} \right] \\
   \dot m = -\frac{c_1}{c_2} \cdot \frac{u}{m}
   \end{array}
   \right.

The auxiliary functions and variables used above are:

.. math::
   w = 1 + f \cos L + g \sin L,\quad s^2 = 1 + h^2 + k^2

Introducing the state, containing both the physical state and the co-states, as:

.. math::
   \mathbf x = [p, f, g, h, k, L, m, \lambda_p, \lambda_f, \lambda_g, \lambda_h, \lambda_k, \lambda_L, \lambda_m]

the parameter vector, containing the gravitational parameter, maximum thrust, effective exhaust velocity, a homotopy continuation parameter, 
and a normalizing Lagrange multiplier (see later), as:

.. math::
   \mathbf p = [\mu, c_1, c_2, \epsilon, \lambda_0]

and the control vector as:

.. math::
   \mathbf u = [u, \hat{\mathbf i}_\tau]

the dynamics can be rewritten in compact form as:

.. math::
   \left\{
   \begin{array}{l}
   \dot{\mathbf x} = \frac{c_1 u(t)}{m} \mathbf B(\mathbf x) \hat{\mathbf i}_\tau + \mathbf D(\mathbf x) \\
   \dot m = -c_2 u(t)
   \end{array}
   \right.

Where:

.. math::
   \sqrt{\frac{\mu}{p}} \mathbf B(\mathbf x) = 
   \begin{bmatrix}
   0 & \frac{2p}{w} & 0 \\
   \sin L & \frac{(1+w)\cos L + f}{w} & -\frac{g}{w}(h \sin L - k \cos L) \\
   -\cos L & \frac{(1+w)\sin L + g}{w} & \frac{f}{w}(h \sin L - k \cos L) \\
   0 & 0 & \frac{1}{w} \frac{s^2}{2} \cos L \\
   0 & 0 & \frac{1}{w} \frac{s^2}{2} \sin L \\
   0 & 0 & \frac{1}{w} (h \sin L - k \cos L)
   \end{bmatrix}

.. math::
   \mathbf D(\mathbf x) = 
   \begin{bmatrix}
   0 \\
   0 \\
   0 \\
   0 \\
   0 \\
   \sqrt{\frac{\mu}{p^3}} w^2
   \end{bmatrix}

As in the Cartesian case, the equations of motion for the TPBVP
are derived from the Hamiltonian, which depends on the chosen optimality:

- **Mass Optimality**:

.. math::
   \mathcal H = \boldsymbol \lambda^\top \dot{\mathbf x} + \lambda_m \dot{m} + \lambda_0 \frac{c_1}{c_2} \left( u + \epsilon \log(u(1 - u)) \right)

- **Time Optimality**:

.. math::
   \mathcal H = \boldsymbol \lambda^\top \dot{\mathbf x} + \lambda_m \dot{m} + \lambda_0 \frac{c_1}{c_2}

The full augmented dynamics is thus obtained by:

.. math::
   \left\{
   \begin{array}{l}
   \dot{\mathbf x} = \frac{\partial \mathcal H}{\partial \boldsymbol \lambda} \\
   \dot{\boldsymbol \lambda} = - \frac{\partial \mathcal H}{\partial \mathbf x}
   \end{array}
   \right.

And by applying Pontryagin’s Minimum Principle:

.. math::
   \mathbf u^* = \arg\min_{\mathbf u \in \mathcal U} \mathcal H(\mathbf x, \boldsymbol \lambda, \mathbf u)

where :math:`\mathcal U` is the admissible set of control inputs (direction and throttle).

Args:
    *optimality* (:class:`pykep.optimality_type`): the optimality principle to be used.

Returns:
    :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(x, dx), ...]
)";
}


std::string propagate_lagrangian_docstring()
{
    return R"(propagate_lagrangian(rv = [[1,0,0], [0,1,0]], tof = pi/2, mu = 1, stm = False)

    Propagates (Keplerian) the state for an assigned time and computes the State Transition Matrix (if requested) using the Lagrangian coefficients.

    Args:
          *rv* (2D array-like): Cartesian components of the initial position vector and velocity [[x0, y0, z0], [v0, vy0, vz0]]. Defaults to [[1,0,0], [0,1,0]].

          *tof* (:class:`float`): time of flight. Defaults to :math:`\frac{\pi}{2}`.

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *stm* (:class:`bool`): requests the computations of the State Transition Matrix

    Returns:
          :class:`tuple` (:class:`list`, :class:`list`): r and v, that is the final position and velocity after the propagation. (if *stm* is False)
          :class:`tuple` (:class:`tuple` (:class:`list`, :class:`list`), :class:`numpy.ndarray` (6,6)): (r,v) and the STM. (if *stm* is True)

    Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> r0 = [1,0,0]
        >>> v0 = [0,1,0]
        >>> tof = np.pi/2
        >>> mu = 1
        >>> [r1,v1], stm = pk.propagate_lagrangian(rv=[r0,v0], tof = tof, mu = mu, stm = True)
        >>> [r1,v1] = pk.propagate_lagrangian(rv=[r0,v0], tof = tof, mu = mu, stm = False)

)";
}

std::string propagate_lagrangian_grid_docstring()
{
    return R"(propagate_lagrangian_grid(rv = [[1,0,0], [0,1,0]], tofs = [pi/2], mu = 1, stm = False)

    This function calls :func:`pykep.propagate_lagrangian` in a loop to provide the results on a time grid containing several epochs. 
    The main advantage is in the API which is convenient and reminiscent of the :func:`heyoka.propagate_grid` interface.

    Note that this function is not necessarily more efficient than calling
    :func:`pykep.propagate_lagrangian` in a loop, since there is no parallelization nor SIMD magic implemented atm. 
    Nevertheless we offer this interface for convenience as it may allow more compact code. 

    Args:
          *rv* (2D array-like): Cartesian components of the initial position vector and velocity [[x0, y0, z0], [v0, vy0, vz0]]. Defaults to [[1,0,0], [0,1,0]].

                    *tofs* (1D array-like): time of flight. Defaults to [:math:`\frac{\pi}{2}`].

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *stm* (:class:`bool`): requests the computations of the State Transition Matrix

    Returns:
          :class:`list` [:class:`tuple` ( :class:`list` , :class:`list` ) ]: For each time of flight: [r,v], that is the final position
          and velocity after the propagation and the flattened stm (if requested).

    Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> r0 = [1,0,0]
        >>> v0 = [0,1,0]
        >>> tofs = [np.pi/2, np.pi, 3*np.pi/4]
        >>> mu = 1
        >>> res = pk.propagate_lagrangian_grid(rv = [r0, v0], tofs = tofs, mu = mu, stm = True)
        >>> rs = [it[0][0] for it in res]
        >>> vs = [it[0][1] for it in res]
        >>> stms = [it[1] for it in res]
)";
}

std::string leg_sf_docstring()
{
    return R"(__init__(rvs = [[1,0,0], [0,1,0]], ms = 1., throttles = [0,0,0,0,0,0], rvf = [[0,1,0], [-1,0,0]], mf = 1., tof = pi/2, max_thrust = 1., veff = 1., mu=1., cut = 0.5)

      This class represents an interplanetary low-thrust transfer between a starting and a final point in the augmented state-space :math:`[\mathbf r, \mathbf v, m]`.
      The low-thrust transfer is described by a sequence of equally spaced impulses as described in:

      Sims, J., Finlayson, P., Rinderle, E., Vavrina, M. and Kowalkowski, T., 2006, August. Implementation of a low-thrust trajectory optimization algorithm for preliminary design. 
      In AIAA/AAS Astrodynamics specialist conference and exhibit (p. 6746).

      The low-thrust transfer will be feasible is the state mismatch equality constraints and the throttle mismatch inequality constraints are satisfied.

      Args:
          *rvs* (2D array-like): Cartesian components of the initial position vector and velocity [[xs, ys, zs], [vxs, vys, vzs]]. Defaults to [[1,0,0], [0,1,0]].

          *ms* (:class:`float`): initial mass. Defaults to 1.

          *throttles* (1D array-like): the Cartesan components of the throttle history [ux1, uy1, uz1, ux2, uy2, uz2, .....]. Defaults to a ballistic, two segments profile [0,0,0,0,0,0].

          *rvf* (2D array-like): Cartesian components of the final position vector and velocity [[xf, yf, zf], [vxf, vyf, vzf]]. Defaults to [[0,1,0], [-1,0,0]].

          *mf* (:class:`float`): final mass. Defaults to 1.

          *tof* (:class:`float`): time of flight. Defaults to :math:`\frac{\pi}{2}`.

          *max_thrust* (:class:`float`): maximum level for the spacecraft thrust. Defaults to 1.

          *veff* (:class:`float`): effective velocity of the propulsion system. Defaults to 1.

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *cut* (:class:`float`): the leg cut, in [0,1]. It determines the number of forward and backward segments. Defaults to 0.5.

      .. note::

        Units need to be consistent. 

      Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> sf = pk.leg.sims_flanagan()
)";
}

std::string leg_sf_alpha_docstring()
{
    return R"(__init__(rvs = [[1,0,0], [0,1,0]], ms = 1., throttles = [0,0,0,0,0,0], talphas = [0,0], rvf = [[0,1,0], [-1,0,0]], mf = 1., tof = pi/2, max_thrust = 1., isp = 1., mu=1., cut = 0.5)

      This class represents an interplanetary low-thrust transfer between a starting and a final point in the augmented state-space :math:`[\mathbf r, \mathbf v, m]`.
      The low-thrust transfer is described by a sequence impulses (not necessarily equally spaced, but centered in the time-intervals defined by *talphas*):

      The low-thrust transfer will be feasible is the state mismatch equality constraints and the throttle mismatch inequality constraints are satisfied.

      Note: NO GRADIENTS impelemented at the moment.

      Args:
          *rvs* (2D array-like): Cartesian components of the initial position vector and velocity [[xs, ys, zs], [vxs, vys, vzs]]. Defaults to [[1,0,0], [0,1,0]].

          *ms* (:class:`float`): initial mass. Defaults to 1.

          *throttles* (1D array-like): the Cartesan components of the throttle history [ux1, uy1, uz1, ux2, uy2, uz2, .....]. Defaults to a ballistic, two segments profile [0,0,0,0,0,0].
          
          *talphas* (1D array-like): the time-intervals where the impulses are centred [ta1, ta2, ..., tnseg]. Defaults to two equally spaced impulses  [tof/2, tof/2].

          *rvf* (2D array-like): Cartesian components of the final position vector and velocity [[xf, yf, zf], [vxf, vyf, vzf]]. Defaults to [[0,1,0], [-1,0,0]].

          *mf* (:class:`float`): final mass. Defaults to 1.

          *tof* (:class:`float`): time of flight. Defaults to :math:`\frac{\pi}{2}`.

          *max_thrust* (:class:`float`): maximum level for the spacecraft thrust. Defaults to 1.

          *isp* (:class:`float`): specific impulse of the propulasion system. Defaults to 1.

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *cut* (:class:`float`): the leg cut, in [0,1]. It determines the number of forward and backward segments. Defaults to 0.5.

      .. note::

        Units need to be consistent. 

      Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> sf = pk.leg.sims_flanagan()
)";
}
std::string leg_sf_rvs_docstring()
{
    return "The initial position vector and velocity: [[xs, ys, zs], [vxs, vys, vzs]].";
};
std::string leg_sf_ms_docstring()
{
    return "Initial mass.";
};
std::string leg_sf_throttles_docstring()
{
    return "The Cartesan components of the throttle history [ux1, uy1, uz1, ux2, uy2, uz2, .....].";
};
std::string leg_sf_talphas_docstring()
{
    return "The time interval of each segment [ta1, ta2, ..., tanseg]. Sums to time-of-flight.";
};
std::string leg_sf_rvf_docstring()
{
    return "The final position vector and velocity: [[xs, ys, zs], [vxs, vys, vzs]].";
};
std::string leg_sf_mf_docstring()
{
    return "Final mass.";
};
std::string leg_sf_tof_docstring()
{
    return "Time of flight.";
};
std::string leg_sf_max_thrust_docstring()
{
    return "Maximum spacecraft thruet.";
};
std::string leg_sf_veff_docstring()
{
    return "Effective velocity of the propulsion system (Isp*G0 in the V units of the dynamics)";
};
std::string leg_sf_mu_docstring()
{
    return "Central body gravitational parameter.";
};
std::string leg_sf_cut_docstring()
{
    return "The leg cut: it determines the number of forward and backward segments.";
};
std::string leg_sf_nseg_docstring()
{
    return "The total number of segments";
};
std::string leg_sf_nseg_bck_docstring()
{
    return "The total number of backward segments";
};
std::string leg_sf_nseg_fwd_docstring()
{
    return "The total number of forward segments";
};
std::string leg_sf_mc_docstring()
{
    return R"(compute_mismatch_constraints()

      In the Sims-Flanagan trajectory leg model, a forward propagation is performed from the starting state as well as a backward from the final state.
      The state values thus computed need to match in some middle control point. This is typically imposed as 7 independent constraints called mismatch-constraints
      computed by this method. 

      Returns:
          :class:`list` [:class:`float`]: The seven mismatch constraints in the same units used to construct the leg.
      
      Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> sf = pk.leg.sims_flanagan()
        >>> sf.compute_mismatch_constraints()
)";
};
std::string leg_sf_tc_docstring()
{
    return R"(compute_throttle_constraints()

      In the Sims-Flanagan trajectory leg model implemented in pykep, we introduce the concept of throttles. Each throttle is defined by three numbers
      :math:`[u_x, u_y, u_z] \in [0,1]` indicating that a certain component of the thrust vector has reached a fraction of its maximum allowed value. 
      As a consequence, along the segment along which the throttle is applied, the constraint  :math:`u_x ^2 + u_y ^2 + u_z^2 = 1`, called a throttle constraint,
      has to be met. 

      Returns:
          :class:`list` [:class:`float`]: The throttle constraints.
      
      Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> sf = pk.leg.sims_flanagan()
        >>  sf.throttles = [0.8]*3
        >>> sf.compute_throttle_constraints()
)";
};
std::string leg_sf_mc_grad_docstring()
{
    return R"(compute_mc_grad()

Computes the gradients of the mismatch constraints. Indicating the initial augmented state with :math:`\mathbf x_s = [\mathbf r_s, \mathbf v_s, m_s]`, the
final augmented state with :math:`\mathbf x_f = [\mathbf r_f, \mathbf v_f, m_f]`, the total time of flight with :math:`T` and introducing the throttle vector
:math:`\mathbf u = [u_{x0}, u_{y0}, u_{z0}, u_{x1}, u_{y1}, u_{z1} ]` and :math:`\mathbf {\tilde u} = [\mathbf u, T]` (note the time of flight at the end), this method computes the following gradients:

.. math::
  \frac{\partial \mathbf {mc}}{\partial \mathbf x_s}  \rightarrow (7\times7)

.. math::
  \frac{\partial \mathbf {mc}}{\partial \mathbf x_f} \rightarrow (7\times7)

.. math::
  \frac{\partial \mathbf {mc}}{\partial \mathbf {\tilde u}} \rightarrow (7\times(3\mathbf{nseg} + 1))

Returns:
    :class:`tuple` [:class:`numpy.ndarray`, :class:`numpy.ndarray`, :class:`numpy.ndarray`]: The three gradients. sizes will be (7,7), (7,7) and (7, 3nseg + 1)

Examples:
  >>> import pykep as pk
  >>> import numpy as np
  >>> sf = pk.leg.sims_flanagan()
  >>  sf.throttles = [0.8]*3
  >>> sf.compute_mc_grad()
)";
};

std::string leg_sf_tc_grad_docstring()
{
    return R"(compute_tc_grad()

Computes the gradients of the throttles constraints. Indicating the total time of flight with :math:`T`  and introducing the throttle vector
:math:`\mathbf u = [u_{x0}, u_{y0}, u_{z0}, u_{x1}, u_{y1}, u_{z1} ]`, this method computes the following gradient:

.. math::
  \frac{\partial \mathbf {tc}}{\partial \mathbf u} \rightarrow (\mathbf{nseg} \times3\mathbf{nseg})

Returns:
    :class:`tuple` [:class:`numpy.ndarray`]: The gradient. Size will be (nseg,nseg*3).

Examples:
  >>> import pykep as pk
  >>> import numpy as np
  >>> sf = pk.leg.sims_flanagan()
  >>  sf.throttles = [0.8]*3
  >>> sf.compute_tc_grad()
)";

};

// ------------------- ZOH LEG DOCSTRINGS -------------------
std::string leg_zoh_docstring() {
  return R"(__init__(state0, controls, state1, tgrid, cut, tas, max_steps=None)

This class implements an interplanetary low-thrust transfer between a starting and final state in the augmented state-space :math:`[\mathbf{r}, \mathbf{v}, m]`. The transfer is modelled as a sequence of non-uniform segments along which a continuous and constant (zero-order hold) control acts. The time intervals defining these segments are also provided in `tgrid`.

The formulation generalises :class:`pykep.leg.sims_flanagan` to arbitrary dynamics and non-uniform time grids. The dynamics are assumed to be zero-order hold and must be provided as compatible Taylor-adaptive integrators (`tas`).

.. note::
   The requirements on the `tas` passed are: a) the first four *heyoka* parameters must be :math:`T, i_x, i_y, i_z`, b) the system dimension must be 7 c) for the variational integrator, variations on the state and the four parameters only are considered. These requirements are all fulfilled by :class:`pykep.ta.zoh_kep`, :class:`pykep.ta.zoh_eq`, :class:`pykep.ta.zoh_cr3bp` and their variational versions.

A transfer is feasible when the state mismatch equality constraints are satisfied. In the intended usage, throttle equality constraints are also enforced to ensure a proper thrust representation as :math:`T \hat{\mathbf{i}}` with :math:`|\hat{\mathbf{i}}| = 1`.

.. math::
   i_x^2 + i_y^2 + i_z^2 = 1, \quad \forall \text{segments}

Examples:
  >>> import numpy as np
  >>> import heyoka as hy
  >>> state0 = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0])
  >>> state1 = np.array([1.2, 0.1, 0.0, 0.0, 0.9, 0.1, 0.95])
  >>> controls = np.array([0.022, 0.7, 0.7, 0.1, 0.025, -0.3, 0.8, 0.4, 0.015, -0.2, 0.8, 0.4])
  >>> tgrid = np.array([0.0, 0.5, 1.0, 1.23])
  >>> ta = pk.ta.get_zoh_eq(tol=1e-16)
  >>> ta_var = pk.ta.get_zoh_eq_var(tol=1e-16)
  >>> leg = pk.leg.zoh(state0, controls, state1, tgrid, cut=0.5, tas=(ta, ta_var))

)";
}

std::string leg_zoh_state0_docstring() {
  return "Initial state vector [r0, v0, m0] (length 7).";
}
std::string leg_zoh_controls_docstring() {
  return "Control parameters [T, i_x, i_y, i_z] for each segment (length 4*nseg).";
}
std::string leg_zoh_state1_docstring() {
  return "Final state vector [r1, v1, m1] (length 7).";
}
std::string leg_zoh_tgrid_docstring() {
  return "Non-uniform time grid (length nseg+1).";
}
std::string leg_zoh_cut_docstring() {
  return "Forward/backward segment split ratio (0 ≤ cut ≤ 1).";
}
std::string leg_zoh_max_steps_docstring() {
  return "Maximum number of steps for the integrator (optional).";
}
std::string leg_zoh_nseg_docstring() {
  return "The total number of segments.";
}
std::string leg_zoh_nseg_fwd_docstring() {
  return "The total number of forward segments.";
}
std::string leg_zoh_nseg_bck_docstring() {
  return "The total number of backward segments.";
}
std::string leg_zoh_mc_docstring() {
  return R"(compute_mismatch_constraints()

Propagates forward from *state0* and backward from *state1* and returns the 7-component
state mismatch at the midpoint.

Returns:
  :class:`list`: Mismatch vector of length 7. All entries are zero for a feasible transfer.
)";
}
std::string leg_zoh_tc_docstring() {
  return R"(compute_throttle_constraints()

Computes the throttle unit-norm constraints :math:`i_x^2 + i_y^2 + i_z^2 - 1` for every segment.

Returns:
  :class:`list`: Constraint values of length nseg. All entries are zero when the direction vector is unit-norm on every segment.
)";
}
std::string leg_zoh_mc_grad_docstring() {
  return R"(compute_mc_grad()

Computes the gradients of the mismatch constraints. Indicating the initial augmented state with :math:`\mathbf x_s = [\mathbf r_s, \mathbf v_s, m_s]`, the final augmented state with :math:`\mathbf x_f = [\mathbf r_f, \mathbf v_f, m_f]`, the time grid as :math:`T_{grid}` and the introducing the control vector :math:`\mathbf u = [T_0, i_{x0}, i_{y0}, i_{z0}, T_1, i_{x1}, i_{y1}, i_{z1}, \ldots]`, this method computes the following gradients:

.. math::
   
   \frac{\partial \mathbf {mc}}{\partial \mathbf x_s}  \rightarrow (7\times7)

.. math::
   
   \frac{\partial \mathbf {mc}}{\partial \mathbf x_f}  \rightarrow (7\times7)

.. math::
   
   \frac{\partial \mathbf {mc}}{\partial \mathbf u}  \rightarrow (7\times(4\mathbf{nseg}))

.. math::
   
   \frac{\partial \mathbf {mc}}{\partial \mathbf T_{grid}}  \rightarrow (7\times(\mathbf{nseg} + 1))

Returns:
  :class:`tuple` [:class:`numpy.ndarray`, :class:`numpy.ndarray`, :class:`numpy.ndarray`, :class:`numpy.ndarray`]: The four gradients. Sizes will be (7,7), (7,7), (7,4nseg), and (7,nseg+1).
)";
}
std::string leg_zoh_tc_grad_docstring() {
  return R"(compute_tc_grad()

Computes the gradients of the throttle constraints. Introducing the control vector as :math:`\mathbf u = [T_0, i_{x0}, i_{y0}, i_{z0}, T_1, i_{x1}, i_{y1}, i_{z1}, ...]`, this method computes the following gradient:

.. math::
   
   \frac{\partial \mathbf {tc}}{\partial \mathbf u} \rightarrow (\mathbf{nseg} \times 4\mathbf{nseg})

Returns:
  :class:`tuple` [:class:`numpy.ndarray`]: The gradient. Size will be (nseg,4nseg).
)";
}

std::string fb_con_docstring()
{
    return R"(fb_con(v_rel_in, v_rel_out, mu, safe_radius)
Alternative signature: fb_con(v_rel_in, v_rel_out, planet)

Computes the constraint violation during a fly-by modelled as an instantaneous rotation of the incoming and outgoing relative velocities (Mivovitch).
The two must be identical in magnitude (equality constraint) and the angle :math:`\alpha` between them must be less than the 
:math:`\alpha_{max}`: the maximum value allowed for that particular *planet* (inequality constraint), as computed from its gravitational
parameter and safe radius using the formula:

.. math::
  \alpha_{max} = - 2 \arcsin\left(\frac{1}{e_{min}}\right)

where:

.. math::
  e_{min} = 1 + V_{\infty}^2 \frac{R_{safe}}{\mu}

.. note::
  This function is often used in the multiple gravity assist low-thrust (MGA-LT) encoding of an interplanetary trajectory where multiple
  low-thrust legs are patched at the fly-by planets by forcing satisfaction for these constraints.

Args:
  *v_rel_in* (:class:`list` (3,)): Cartesian components of the incoming relative velocity.

  *v_rel_out* (:class:`list` (3,)): Cartesian components of the outgoing relative velocity.

  *mu* (:class:`float`): planet gravitational parameter

  *safe_radius* (:class:`float`): planet safe radius

  *planet* (:class:`~pykep.planet`): planet (in which case *mu* and *safe_radius* will be extracted from this object). Note: this signature is slower and to be avoided.

Returns:
  :class:`tuple` [:class:`float`, :class:`float`]: The equality constraint violation (defined as the difference between the squared velocities) and the inequality constraint violation 
  (negative if satisfied).

Examples:
  >>> import pykep as pk
  >>> eq, ineq = pk.fb_con(v_rel_in = [10.,1.,-4.], v_rel_out = [10.,1.,-4.], mu = 1., safe_radius = 1.)

)";
};

std::string fb_dv_docstring()
{
    return R"(fb_dv(v_rel_in, v_rel_out, mu, safe_radius)
Alternative signature: fb_dv(v_rel_in, v_rel_out, planet)

Computes the :math:`\Delta V` necessary to perform a fly-by modelled as an instantaneous rotation of the incoming and outgoing
relative velocities (Mivovitch). If planetary gravity is not enough to patch the incoming and outcoming conditions
(i.e. the :func:`~pykep.fb_con()` returns some constraint violation) a :math:`\Delta V` is assumed
at the end of the planetocentric hyperbola.

.. note::
  This function is often used in the multiple gravity assist (MGA) encoding of an interplanetary trajectory where multiple
  Lambert arcs are patched at the fly-by planets by applying a hopefully vanishing :math:`\Delta V`.

Args:
  *v_rel_in* (:class:`list` (3,)): Cartesian components of the incoming relative velocity.

  *v_rel_out* (:class:`list` (3,)): Cartesian components of the outgoing relative velocity.

  *mu* (:class:`float`): planet gravitational parameter

  *safe_radius* (:class:`float`): planet safe radius

  *planet* (:class:`~pykep.planet`): planet (in which case *mu* and *safe_radius* will be extracted from this object). Note: this signature is slower and to be avoided.

Returns:
  :class:`float`: The magnitude of the required :math:`\Delta V`

Examples:
  >>> import pykep as pk
  >>> DV = pk.fb_dv(v_rel_in = [10.,1.,-4.], v_rel_out = [10.,1.,-4.], mu = 1., safe_radius = 1.)
)";
};

std::string fb_vout_docstring()
{
    return R"(fb_vout(v_in, v_pla, rp, beta, mu)

Propagates incoming conditions through a planetary encounter (fly-by) assuming an instantaneous rotation of magnitude :math:`\delta` of 
the incoming and outgoing relative velocities (Mivovitch). The planetocentric hyperbola (hence the angle :math:`\delta`) is fully
determined by the incoming condition :math:`v_{\infty} = \mathbf v_{in} - \mathbf v_{pla}` as well as by its pericenter :math:`r_p`
and an angle :math:`\beta` defining the orientation of the orbital plane. Eventually the outgoing conditions are computed as:

.. math::
  \mathbf v_{out} = \mathbf v_{in} + \mathbf v_{\infty}^{out}

.. math::
  \mathbf v_{\infty}^{out} = |\mathbf v_{\infty}^{in}|
  \left(
  \cos\delta\hat{\mathbf b}_1
  +\sin\delta\cos\beta\hat{\mathbf b}_2
  +\sin\delta\sin\beta\hat{\mathbf b}_3
  \right)

where the :math:`[\hat{\mathbf b}_1, \hat{\mathbf b}_2, \hat{\mathbf b}_3]` frame is defined by the incoming 
relative velocity (normalized), :math:`\hat{\mathbf b}_1 \times \mathbf v_{pla}` (normalized) and completing the right handed frame.

.. note::
  This function is often used in the multiple gravity assist with DSM (MGA-DSM) encoding of an interplanetary trajectory where multiple
  impulsive manoeuvres are allowed between planetary fly-bys.

Args:
  *v_in* (:class:`list` (3,)): Cartesian components of the incoming velocity 
  (note: this is NOT the relative velocity, rather the absolute and in the same frame as *v_pla*).

  *v_pla* (:class:`list` (3,)): Cartesian components of the planet velocity.

  *rp* (:class:`float`): planetocentric hyperbola pericenter radius.

  *beta* (:class:`float`): planetocentric hyperbola plane angle.

  *mu* (:class:`float`): planet gravitational parameter.

Returns:
  :class:`list` (3,): The outgoing velocity in the frame of *v_pla* (inertial).

Examples:
  >>> import pykep as pk
  >>> DV = pk.fb_vout(v_in = [1.,1.,1], v_pla = [10.,1.,-4.], rp = 1., mu = 1., beta=1.2)
)";
};

} // namespace pykep
