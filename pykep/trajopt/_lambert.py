import numpy as np
import math
import random
from pykep.core import fb_vel

class lambert_problem_multirev:
    
    r"""
    This class converts a lambert_problem instance - a number of solutions to a multi revolution Lambert problem
    to an instance representing only "the best" solution. Criteria is the delta velocity of the incoming
    velocity v_in to the outgoing solution velocity get_v1(). Can be used as a replacement for lambert_problem
    in _mga.py and _mga_1dsm.py. lambert_problem_multirev delivers an optimal unique solution comparing all
    Lambert solutions of a multi revolution lambert problem which can improve optimization results for trajectories
    including inner planets.
    """

    def __init__(self, v_in, lambert_problem):
        best_i = 0        
        n = len(lambert_problem.get_v1())
        if n > 0: 
            best_dv = math.inf
            for i in range(n):
                dv = np.linalg.norm([a - b for a, b in zip(lambert_problem.get_v1()[i], v_in)])
                if dv < best_dv:
                    best_dv = dv
                    best_i = i 
        self.best_i = best_i
        self.lambert_problem = lambert_problem
        
    def get_v1(self):
        return [self.lambert_problem.get_v1()[self.best_i]]
   
    def get_v2(self):
        return [self.lambert_problem.get_v2()[self.best_i]]

    def get_r1(self):
        return self.lambert_problem.get_r1()  

    def get_r2(self):
        return self.lambert_problem.get_r2()  

    def get_mu(self):
        return self.lambert_problem.get_mu()

    def get_x(self):
        return [self.lambert_problem.get_x()[self.best_i]]

    def get_iters(self):
        return [self.lambert_problem.get_iters()[self.best_i]]

    def get_tof(self):
        return self.lambert_problem.get_tof()

    def get_Nmax(self):
        return self.lambert_problem.get_Nmax()

class lambert_problem_multirev_ga:    
    r"""
    This class converts a lambert_problem instance - a number of solutions to a multi revolution Lambert problem
    to an instance representing only "the best" solution. Criteria is the delta velocity of the incoming
    velocity v_in to the outgoing solution velocity get_v1() considering a GA maneuver at a planet.
    Can be used as a replacement for lambert_problem in _mga.py and _mga_1dsm.py. lambert_problem_multirev delivers 
    an optimal unique solution comparing all Lambert solutions of a multi revolution lambert problem which can 
    improve optimization results for trajectories including inner planets.
    """
 
    def __init__(self, v_in, lp, planet, v_planet):
        best_i = 0        
        n = len(lp.get_v1())
        if n > 0: 
            best_dv = math.inf           
            for i in range(n):
                vin = [a - b for a, b in zip(v_in, v_planet)]
                vout = [a - b for a, b in zip(lp.get_v1()[i], v_planet)]
                dv = fb_vel(vin, vout, planet)
                if dv < best_dv:
                    best_dv = dv
                    best_i = i 
        self.best_i = best_i
        self.lambert_problem = lp
        
    def get_v1(self):
        return [self.lambert_problem.get_v1()[self.best_i]]
   
    def get_v2(self):
        return [self.lambert_problem.get_v2()[self.best_i]]

    def get_r1(self):
        return self.lambert_problem.get_r1()  

    def get_r2(self):
        return self.lambert_problem.get_r2()  

    def get_mu(self):
        return self.lambert_problem.get_mu()

    def get_x(self):
        return [self.lambert_problem.get_x()[self.best_i]]

    def get_iters(self):
        return [self.lambert_problem.get_iters()[self.best_i]]

    def get_tof(self):
        return self.lambert_problem.get_tof()

    def get_Nmax(self):
        return self.lambert_problem.get_Nmax()

class lambert_problem_stochastic:
    r"""
    This class converts a lambert_problem instance - a number of solutions to a multi revolution Lambert problem -
    to an instance representing a single "good" solution. Criteria is the delta velocity of the incoming
    velocity v_in to the outgoing solution velocity get_v1(). Can be used as a replacement for lambert_problem
    in _mga.py and _mga_1dsm.py. The used criteria is not the only one relevant for optimization problems aiming at 
    reducing the total delta velocity. I doesn't consider the incoming solution velocity get_v2() because 
    the corresponding outgoing target velocity is not available yet. By randomly choosing a "good" solution 
    the value of the objective function calling lambert_problem_stochastic is no longer deterministic. 
    Combined with a non derivative optimization algorithm this way search can be replaced by optimization. 

   .. note::

    The formula used to randomly skip over superior solutions was tested for _tandem and _messenger and GTOC1. 
    Other problems may require different parameters.    

    """
    
    def __init__(self, v_in, lambert_problem):
        best_i = 0        
        n = len(lambert_problem.get_v1())
        if n > 0: 
            best_dv = math.inf
            for i in range(n):
                # convert to km/s
                dv = 0.001 * np.linalg.norm([a - b for a, b in zip(lambert_problem.get_v1()[i], v_in)])
                if dv < best_dv:
                    # skip improvement randomly if dv is only slightly better
                    if i > 0 and random.random() < 1.0 / (1.0 + 2*(0.001*best_dv - dv)**2):
                        continue
                    best_dv = dv
                    best_i = i

        self.best_i = best_i
        self.lambert_problem = lambert_problem
        
    def get_v1(self):
        return [self.lambert_problem.get_v1()[self.best_i]]
   
    def get_v2(self):
        return [self.lambert_problem.get_v2()[self.best_i]]

    def get_r1(self):
        return self.lambert_problem.get_r1()

    def get_r2(self):
        return self.lambert_problem.get_r2()

    def get_mu(self):
        return self.lambert_problem.get_mu()

    def get_x(self):
        return [self.lambert_problem.get_x()[self.best_i]]

    def get_iters(self):
        return [self.lambert_problem.get_iters()[self.best_i]]

    def get_tof(self):
        return self.lambert_problem.get_tof()

    def get_Nmax(self):
        return self.lambert_problem.get_Nmax()
