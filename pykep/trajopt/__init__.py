## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)## 
## This file is part of the kep3 library.## 
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
User defined problems (compatible to pagmo) that represent interplanetary optimization problems 
"""

# Direct methods for low-thrust problems
from ._sf_point2point import sf_point2point
from ._sf_pl2pl import sf_pl2pl
from ._sf_pl2pl_alpha import sf_pl2pl_alpha
from ._zoh_point2point import zoh_point2point
from ._zoh_ss_point2point import zoh_ss_point2point
from ._zoh_pl2pl import zoh_pl2pl

# Evolutionary encodings for high energy transfers (chemical propulsion)
from ._mga import mga
from ._mga_1dsm import mga_1dsm
from ._pl2pl_N_impulses import pl2pl_N_impulses

# Indirect methods for low-thrust problems
from ._pontryagin_cartesian import pontryagin_cartesian_mass, pontryagin_cartesian_time
from ._pontryagin_equinoctial import pontryagin_equinoctial_mass, pontryagin_equinoctial_time
from ._mim import mim_from_hop

# MIT (multiple Impulse Trajectories)
from ._primer_vector import primer_vector, primer_vector_surrogate
from ._min_Bu_bu import minBu_bu_p, minBu_bu

# The launchers models
from ._launchers import _launchers
launchers = _launchers()

# The interplanetary trajectory gym
from . import gym