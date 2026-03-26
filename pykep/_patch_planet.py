## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)
## This file is part of the kep3 library.
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

from .core import planet

def _planet_extract(self, t):
    """Extract the user-defined planet.

    This method allows to extract a reference to the user-defined planet (UDPLA) stored within the
    :class:`~pykep.planet` instance. The behaviour of this function depends on the value
    of *t* (which must be a :class:`type`) and on the type of the internal UDPLA:

    * if the type of the UDPLA is *t*, then a reference to the UDP will be returned,
    * if *t* is :class:`object` and the UDP is a Python object (as opposed to an
      exposed C++ planet), then a reference to the
      UDPLA will be returned (this allows to extract a Python UDPLA without knowing its type),
    * otherwise, :data:`None` will be returned.

    Args:
        *t* (:class:`type`): the type of the user-defined planet to extract

    Returns:
        a reference to the internal user-defined planet, or :data:`None` if the extraction fails

    Raises:
        *TypeError*: if *t* is not a :class:`type`

    Examples:
        >>> import pykep as pk
        >>> udpla = pk.udpla.keplerian(pk.epoch(0), [1,0,0,0,0,0], 1)
        >>> pla = pk.planet(udpla)
        >>> type(pla.extract(pk.udpla.keplerian)) 
        <class 'udpla._keplerian'>
        >>> pla.extract(pk.udpla.jpl_lp) is None
        True
        >>> class my_udpla:
        ...     def eph(self, ep):
        ...         return [[1,0,0],[0,1,0]]
        ...     def get_name(self):
        ...         return "my_udpla"
        ...     def get_mu_central_body(self):
        ...         return 1.
        >>> pla2 = pk.planet(my_udpla())
        >>> p2.extract(object) # doctest: +SKIP
        <__main__.my_udpla at 0x7ff68b63d210>
        >>> pla2.extract(my_udpla) # doctest: +SKIP
        <__main__.my_udpla at 0x7f8f7241c350>
        >>> pla2.extract(pk.udpla.keplerian) is None
        True
    """
    if not isinstance(t, type):
        raise TypeError("the 't' parameter must be a type")
    # This happens if the udpla is cpp
    if hasattr(t, "_pykep_cpp_udpla"):
        return self._cpp_extract(t())
    # Else we extract the UDPLA from the python class, check its type 
    # and see if we want to return None or not as to have a consistent behaviour
    udpla = self._py_extract()
    if (type(udpla) == t):
        return udpla
    else:
        return None


def _planet_is(self, t):
    """Check the type of the user-defined planet.

    This method returns :data:`False` if :func:`~pykep.planet.extract()` returns
    :data:`None`, and :data:`True` otherwise.

    Args:
        *t*  (:class:`type`): the type that will be compared to the type of the UDPLA

    Returns:
        bool: whether the UDPLA is of type *t* or not

    Raises:
        *unspecified*: any exception thrown by :func:`~pykep.planet.extract()`

    """
    return not self.extract(t) is None


# Do the actual patching.
setattr(planet, "extract", _planet_extract)
setattr(planet, "is_", _planet_is)