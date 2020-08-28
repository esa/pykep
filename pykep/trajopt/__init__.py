"""
This module contains interplanetary trajectories problems in the form of pygmo.problem classes. The use
of pygmo.archipelago and pygmo.island with these may result in slow execution as the parallel execution will
be handled by python multiprocessing module, which is far slower than boost::threads (used instead when the pygmo problems
are implemented in c++)
"""
from pykep.trajopt._lt_margo import lt_margo
from pykep.trajopt._mga_1dsm import mga_1dsm
from pykep.trajopt._mga import mga
from pykep.trajopt._mga_lt_nep import mga_lt_nep
from pykep.trajopt._mr_lt_nep import mr_lt_nep
from pykep.trajopt._pl2pl_N_impulses import pl2pl_N_impulses
# The launchers models
from pykep.trajopt._launchers import _launchers
launchers = _launchers()
# The open traj gym
from . import gym 
from pykep.trajopt._indirect import indirect_or2or, indirect_pt2or, indirect_pt2pt, indirect_pt2pl
from pykep.trajopt._direct import direct_pl2pl
from pykep.trajopt._lambert import lambert_problem_multirev, lambert_problem_stochastic
