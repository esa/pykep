# Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
#                          Advanced Concepts Team, European Space Agency (ESA)
#
# This file is part of the pykep library.
#
# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Pykep is a toolbox for interplanetary trajectory design developed by ESA's Advanced Concepts Team. Its main
 purpose is fast prototyping of research ideas, and is not intended for operational usage.

 Some important conventions followed:
 1 - All units expected are S.I. (m,sec,kg,N) unless explicitly stated.
 2 - The default set of osculating orbital parameters is, in this order: [sma, ecc, incl, W, w, f], where f is the true anomaly
 3 - The default option to represent epochs as floats is the modified julian date 2000 (MJD2000). By default, time durations are in days."""

# Version setup.
from ._version import __version__

del _version

import os as _os

if _os.name == "posix":
    # NOTE: on some platforms Python by default opens extensions
    # with the RTLD_LOCAL flag, which creates problems because
    # public symbols used by heyoka (e.g., sleef functions, quad
    # precision math) are then not found by the LLVM jit machinery.
    # Thus, before importing core, we temporarily flip on the
    # RTLD_GLOBAL flag, which makes the symbols visible and
    # solves these issues. Another possible approach suggested
    # in the llvm discord is to manually and explicitly add
    # libheyoka.so to the DL search path:
    # DynamicLibrarySearchGenerator::Load(“/path/to/libheyoka.so”)
    # See:
    # https://docs.python.org/3/library/ctypes.html
    import ctypes as _ctypes
    import sys as _sys

    _orig_dlopen_flags = _sys.getdlopenflags()
    _sys.setdlopenflags(_orig_dlopen_flags | _ctypes.RTLD_GLOBAL)

    try:
        # Importing cpp functionalities
        from .core import *
    finally:
        # Restore the original dlopen flags whatever
        # happens.
        _sys.setdlopenflags(_orig_dlopen_flags)

        del _ctypes
        del _sys
        del _orig_dlopen_flags
else:
    # Importing cpp functionalities
    from .core import *

del _os

import numpy as _np
_np.set_printoptions(legacy='1.13')
del _np

# Importing user defined planets
from . import udpla

# Importing trajectory legs
from . import leg

# Importing Taylor adaptive integrators
from . import ta

# Importing the python utils
from .utils import *

# Patch the planet class.
from . import _patch_planet

# Patch the epoch class.
from . import _patch_epoch

# Import the plot module
from . import plot

# Import the trajopt module
from . import trajopt
mim_from_hop = trajopt.mim_from_hop # we want mim also in the same namespace as mima(s)

# We import the unit test submodule
from . import test
