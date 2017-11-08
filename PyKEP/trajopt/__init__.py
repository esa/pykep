"""
This module contains interplanetary trajectories problems in the form of PyGMO.problem classes. The use
of PyGMO.archipelago and PyGMO.island with these may result in slow execution as the parallel execution will
be handled by python multiprocessing module, which is far slower than boost::threads (used instead when the PyGMO problems
are implemented in c++)
"""
from PyKEP import __extensions__

if (__extensions__['pygmo']):
    from PyKEP.trajopt._lt_margo import lt_margo
    #from PyKEP.trajopt._mga_1dsm import mga_1dsm
    #from PyKEP.trajopt._mga_lt_nep import mga_lt_nep
    #from PyKEP.trajopt._mr_lt_nep import mr_lt_nep
    #from PyKEP.trajopt._pl2pl_N_impulses import pl2pl_N_impulses

if (__extensions__['pygmo'] and __extensions__['mplot3d'] and __extensions__['scipy']):
    from PyKEP.trajopt._indirect import indirect_or2or, indirect_pt2or, indirect_pt2pt

if (__extensions__['pygmo'] and __extensions__['mplot3d']):
    from PyKEP.trajopt._direct import direct_pl2pl
