from .. import core as _core
from ._tle import tle
from ._spice import spice, de440s

# expose original names for pickle and imports
_keplerian = _core._keplerian
_null_udpla = _core._null_udpla
_jpl_lp = _core._jpl_lp
_vsop2013 = _core._vsop2013

# alias with proper module and name for docs & usage
keplerian = _keplerian
keplerian.__name__ = "keplerian"
keplerian.__module__ = "pykep.udpla"

null_udpla = _null_udpla
null_udpla.__name__ = "null_udpla"
null_udpla.__module__ = "pykep.udpla"

jpl_lp = _jpl_lp
jpl_lp.__name__ = "jpl_lp"
jpl_lp.__module__ = "pykep.udpla"

vsop2013 = _vsop2013
vsop2013.__name__ = "vsop2013"
vsop2013.__module__ = "pykep.udpla"

del _core
