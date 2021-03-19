from math import pi, sqrt
import math
from typing import Any, List, Tuple

from pykep import AU, DAY2SEC, SEC2DAY, DEG2RAD, RAD2DEG, epoch, ic2par
from pykep.core import fb_vel, lambert_problem
from pykep.trajopt import lambert_problem_multirev_ga
from pykep.planet import jpl_lp

import numpy as np
from pykep.trajopt._resonance import resonance
from pykep.trajopt._rvt import rvt_planet, rotate_vector, rvt

class _solo_mgar_udp:
    """
    See https://www.esa.int/Science_Exploration/Space_Science/Solar_Orbiter . 
    This UDP represents some of the objectives of the original solar orbiter mission:
    - ballistic transfers
    - multiple GA manoevers at Venus and Earth
    - limited mission time
    - inclination > 33 deg at the end of the mission
    - low perhelion near 0.28 AU during the EMP phase. 
    thereby fulfilling its constraints. 
    To keep it as simple as possible, this UDP currently works with the fixed original planet sequence EVVEVVVVVV, 
    more flexible variants may be created later. Each EV / VE transfer is performed as MGA, 
    each VV transfer is performed as resonance. This means, for a free beta parameter the corresponding orbiter
    period is computed. From the given resonance alternatives the one is chosen, which minimizes the timing error 
    - how much earlier / later than venus the orbiter arrives after the resonance transfer.
    Instead of "caching" of position / velocity data 'rvt' objects are used: Orbit representations 
    storing position, velocity, time and mu.  
    The UDP provides a hand crafted 'unconstrain' activated if use_constraints == False which was used to compute the
    reference results at https://github.com/dietmarwo/fast-cma-es/blob/master/examples/data/solo_results.txt
    An example for test/optimize is here: https://gist.github.com/dietmarwo/5b6bbadbd696d36806177be2df3c6bb7
    """

    def __init__(
        self,
        t0_limits=[7000, 8000],
        tof_limits=[[50, 420], [50, 400], [50, 400]],
        max_revs: int=2,
        resonances=
            [[[1, 1], [5, 4], [4, 3]],
            [[1, 1], [5, 4], [4, 3]],
            [[1, 1], [5, 4], [4, 3]],
            [[4, 3], [3, 2], [5, 3]],
            [[4, 3], [3, 2], [5, 3]],
            [[4, 3], [3, 2], [5, 3]]],
        safe_distance=350000,
        max_mission_time=11.0 * 365.25,
        max_dv0=5600,
        multi_objective=False,
        use_constraints=True
    ):
        """
        Args:
            - t0_limits (``list`` of ``float``): start time bounds. 
            - tof_limits (``list`` of ``list`` of ``float``): time of flight bounds. 
            - max_revs (``int``): maximal number of revolutions for Lambert transfer.
            - resonances (``list`` of ``list`` of ``int``): resonance options. 
            - safe_distance: (``float``): safe distance from planet at GA maneuver in m.
            - max_mission_time: (``float``): max mission time in days.
            - max_dv0: (``float``): maximal delta velocity at start.
            - multi_objective (``bool``): when True the problem fitness will return separate objectives
            - use_constraints (``bool``): when False the problem fitness will use a customized "unconstrain" function
        """
        
        self._safe_distance = safe_distance
        self._max_mission_time = max_mission_time
        self._max_dv0 = max_dv0
        self._min_beta = -math.pi
        self._max_beta = math.pi
        
        self._earth = jpl_lp("earth")
        self._venus = jpl_lp("venus")
        self._seq = [self._earth, self._venus, self._venus, self._earth, self._venus,
                     self._venus, self._venus, self._venus, self._venus, self._venus]
        
        assert len(self._seq) - 4 == len(resonances)  # one resonance option selection for each VV sequence

        self._resonances = resonances
        self._t0 = t0_limits
        self._tof = tof_limits
        self._max_revs = max_revs       
        self._multi_objective = multi_objective
        self._use_constraints = use_constraints

        self._n_legs = len(self._seq) - 1
 
        # initialize data to compute heliolatitude
        t_plane_crossing = epoch(7645)
        rotation_axis = self._seq[0].eph(t_plane_crossing)[0]
        self._rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        self._theta = -7.25 * DEG2RAD  # fixed direction of the rotation
         
    def _compute_dvs(self, x: List[float],
                     lps=None,
                     resos=None
    ) -> Tuple[
        List[Any],
        List[Any],
        List[Any],
        List[float],
        List[float],
    ]:
        
        t0 = x[0] * DAY2SEC  # only direct encoding
        tof01 = x[1] * DAY2SEC
        tof23 = x[2] * DAY2SEC
        tof34 = x[3] * DAY2SEC
        betas = x[4:]
       
        rvt_outs = []
        rvt_ins = [None]  # no rvt_in at first planet
        rvt_pls = []
        dvs = []
        used_resos = []
        reso_dts = []
        
        rvt_pls.append(rvt_planet(self._seq[0], t0))
 
        _dv_mga(self._seq[0], self._seq[1], tof01, self._max_revs, rvt_outs, rvt_ins, rvt_pls, dvs, lps)
         
        compute_resonance(self._seq[1], self._resonances[0], betas[0], self._safe_distance,
                           used_resos, reso_dts, rvt_outs, rvt_ins, rvt_pls, dvs, resos)
        
        _dv_mga(self._seq[2], self._seq[3], tof23, self._max_revs, rvt_outs, rvt_ins, rvt_pls, dvs, lps)
  
        _dv_mga(self._seq[3], self._seq[4], tof34, self._max_revs, rvt_outs, rvt_ins, rvt_pls, dvs, lps)

        for i in range(1, 6): 
            compute_resonance(self._seq[i + 3], self._resonances[i], betas[i], self._safe_distance,
                               used_resos, reso_dts, rvt_outs, rvt_ins, rvt_pls, dvs, resos)
        
        n = len(rvt_ins)
        assert len(rvt_outs) == n - 1
        assert len(rvt_pls) == n 
        assert len(rvt_ins) == n  
        assert len(dvs) == n - 1
        return rvt_outs, rvt_ins, rvt_pls, reso_dts, dvs

    # Objective function
    def fitness(self, x):
        if len(x) != len(self.get_bounds()[0]):
            raise ValueError(
                "Expected "
                +str(len(self.get_bounds()[0]))
                +" parameters but got "
                +str(len(x))
            )

        rvt_outs, rvt_ins, _, reso_dts, dvs = self._compute_dvs(x)
        # compute final flyby and resulting trajectory
        rvt_out = rvt_outs[-1].rotate(self._rotation_axis, self._theta)  # rotate
        _, _, incl, _, _, _ = rvt_out.kepler()
        # orbit should be as polar as possible, but we do not care about prograde/retrograde
        corrected_inclination = abs(abs(incl) % pi - pi / 2) * RAD2DEG

        # check perihelion and aphelion bounds during the flight
        emp_perhelion = 2 * AU
        min_sun_distance = 2 * AU
        max_sun_distance = 0

        for i in range(len(rvt_outs)):
            orb = rvt_outs[i]
            tof = orb.tof(rvt_ins[i + 1])
            transfer_a, transfer_e, _, _, _, _ = orb.kepler()
            transfer_period = 2 * pi * sqrt(transfer_a ** 3 / orb._mu)
            perhelion = transfer_a * (1 - transfer_e)
            # update min and max sun distance
            if i >= len(rvt_outs) - 3:
                emp_perhelion = min(emp_perhelion, perhelion)
            min_sun_distance = min(min_sun_distance, perhelion)
            if tof > transfer_period:
                max_sun_distance = max(max_sun_distance, transfer_a * (1 + transfer_e))

        time_all = SEC2DAY * (rvt_ins[-1]._t - rvt_outs[0]._t)
        time_limit = self._max_mission_time  
        # wrong reso timing in seconds
        reso_penalty = np.sum(reso_dts)   
        # permit start dv _max_dv0  
        dvs[0] = max(0, dvs[0] - self._max_dv0) 
        
        if not self._use_constraints:
            # use customized 'unconstrain' to obtain a single objective
            dv_val = np.sum(dvs)                   
            time_val = time_all
            if time_val > time_limit:
                time_val += 10 * (time_val - time_limit) 
            distance_val = min_sun_distance / AU
            # penalty for perhelion < 0.28
            if distance_val < 0.28:
                distance_val += 10 * (0.28 - distance_val);           
            # wrong minimal / maximal distance
            distance_penalty = max(0, 0.28 - min_sun_distance / AU)  
            distance_penalty += max(0, max_sun_distance / AU - 1.2)                
            value = (100 * dv_val + 
                    100 * corrected_inclination + 
                    5000 * (max(0, distance_val - 0.28)) + 
                    5000 * (max(0, emp_perhelion / AU - 0.28)) + 
                    0.5 * time_val + 
                    reso_penalty + 
                    50000 * distance_penalty
                    )
        else:
            value = corrected_inclination
   
        return (
            # objectives
            [value]
            +[emp_perhelion / AU] * self._multi_objective
            +[time_all] * self._multi_objective
            # constraints
            +[time_all - time_limit] * self._use_constraints
            +[np.sum(dvs) - 0.1] * self._use_constraints
            +[reso_penalty - 0.1] * self._use_constraints
            +[0.28 - min_sun_distance / AU] * self._use_constraints
            +[max_sun_distance / AU - 1.2] * self._use_constraints
            )
 
    def get_nobj(self):
        return self._multi_objective * 2 + 1
    
    def get_nic(self):
        return self._use_constraints * 5
    
    def get_bounds(self):
        t0 = self._t0
        tof = self._tof
        n_legs = self._n_legs
        
        lb = [t0[0]]
        ub = [t0[1]]
        lb += [t[0] for t in tof]
        ub += [t[1] for t in tof]
        
        nbetas = n_legs - len(tof)
        lb += [self._min_beta] * nbetas
        ub += [self._max_beta] * nbetas
        return (lb, ub)

    def pretty(self, x):
        lambert_legs = []
        resos = []
        rvt_outs, rvt_ins, rvt_pls, _, dvs = self._compute_dvs(x, lambert_legs, resos)  
        rvt_outs = [rvt.rotate(self._rotation_axis, self._theta) for rvt in rvt_outs]
        rvt_ins[1:] = [rvt.rotate(self._rotation_axis, self._theta) for rvt in rvt_ins[1:]]
        rvt_pls = [rvt.rotate(self._rotation_axis, self._theta) for rvt in rvt_pls]
            
        ep = [epoch(rvt_pl._t * SEC2DAY) for rvt_pl in rvt_pls]
        b_legs = [[rvt_out._r, rvt_out._v] for rvt_out in rvt_outs]
        Vinfx, Vinfy, Vinfz = [
            a - b for a, b in zip(b_legs[0][1], self._seq[0].eph(ep[0])[1])
        ]
        common_mu = rvt_outs[0]._mu
      
        lambert_indices = [lam.best_i for lam in lambert_legs]
        
        transfer_ang = _angle(rvt_outs[0]._r, rvt_outs[1]._r)
  
        print("Multiple Gravity Assist (MGA) + Resonance problem: ")
        print("Planet sequence: ", [pl.name for pl in self._seq])

        print("Departure: ", self._seq[0].name)
        print("\tEpoch: ", ep[0], " [mjd2000]")
        print("\tSpacecraft velocity: ", b_legs[0][1], "[m/s]")
        print("\tLaunch velocity: ", [Vinfx, Vinfy, Vinfz], "[m/s]")
        _, _, transfer_i, _, _, _ = ic2par(*(b_legs[0]), common_mu)
        print("\tTransfer Angle: ", np.degrees(transfer_ang), "deg")
        print("\tOutgoing Inclination:", transfer_i * RAD2DEG, "[deg]")
        print("\tNumber of Revolutions:", int((lambert_indices[0] + 1) / 2))
        print("\tLambert Index:", int(lambert_indices[0]))
        
        lambert_i = 0
        reso_i = 0
        for i in range(1, len(self._seq) - 1):
            pl = self._seq[i]
            e = ep[i]
            dv = dvs[i]
            leg = b_legs[i]
            rtv_in = rvt_ins[i]
            rtv_out = rvt_outs[i]
            rtv_pl = rvt_pls[i]
            vr_in = [a - b for a, b in zip(rtv_in._v, rtv_pl._v)]
            vr_out = [a - b for a, b in zip(rtv_out._v, rtv_pl._v)]
            v_inf = np.linalg.norm(vr_out)
            deflection = _angle(vr_in, vr_out)
            transfer_ang = _angle(rtv_out._r, rvt_outs[i + 1]._r) if i < len(self._seq) - 2 else 0 
            print("Fly-by: ", pl.name)
            print("\tEpoch: ", e, " [mjd2000]")
            print("\tDV: ", dv, "[m/s]")
            print("\tV_inf: ", v_inf, "[m/s]")
            print("\tTransfer Angle: ", np.degrees(transfer_ang), "deg")
            print("\tGA deflection: ", np.degrees(deflection), "deg")
            eph = [rotate_vector(v, self._rotation_axis, self._theta) for v in pl.eph(e)]
            if i < len(self._seq) - 1: 
                assert np.linalg.norm([a - b for a, b in zip(leg[0], eph[0])]) < 0.01
            _, _, transfer_i, _, _, _ = ic2par(eph[0], leg[1], common_mu)
            print("\tOutgoing Inclination:", transfer_i * RAD2DEG, "[deg]")
            if pl != self._seq[i - 1]:  # no lamberts for resonances
                print("\tLambert Index:", str(lambert_indices[lambert_i]))
                lambert_i += 1
            else:  # resonance at Venus
                print("\tResonance:", str(resos[reso_i]._resonance))
                print("\tResonance time error:", str(resos[reso_i]._timing_error) + " sec")
                reso_i += 1
               
        print("Final Fly-by: ", self._seq[-1].name)
        print("\tEpoch: ", ep[-1], " [mjd2000]")
        print("\tSpacecraft velocity: ", rvt_outs[-1]._v, "[m/s]")
        print("\tBeta: ", x[-1])
        print("\tr_p: ", self._seq[-1].radius + self._safe_distance)

        print("Resulting Solar orbit:")
        a, e, i, _, _, _ = rvt_outs[-1].kepler()
        print("Perihelion: ", (a * (1 - e)) / AU, " AU")
        print("Aphelion: ", (a * (1 + e)) / AU, " AU")
        print("Inclination: ", i * RAD2DEG, " degrees")
        T = [SEC2DAY * (rvt_outs[i + 1]._t - rvt_outs[i]._t) for i in range(len(rvt_outs) - 1)]
        print("Time of flights: ", T, "[days]")
        
    def plot(self, x, axes=None, units=AU, N=60):
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from pykep.orbit_plots import plot_planet

        rvt_outs, rvt_ins, rvt_pls, _, _ = self._compute_dvs(x)    
        rvt_outs = [rvt.rotate(self._rotation_axis, self._theta) for rvt in rvt_outs]
        rvt_ins[1:] = [rvt.rotate(self._rotation_axis, self._theta) for rvt in rvt_ins[1:]]
        rvt_pls = [rvt.rotate(self._rotation_axis, self._theta) for rvt in rvt_pls]
          
        ep = [epoch(rvt_pl._t * SEC2DAY) for rvt_pl in rvt_pls]
 
        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams["legend.fontsize"] = 10
            fig = plt.figure()
            axes = fig.gca(projection="3d")

        plt.xlim([-1, 1])
        # planets
        for pl, e in zip(self._seq, ep):
            plot_planet(
                pl, e, units=units, legend=True, color=(0.7, 0.7, 1), axes=axes
            )
            
        # lamberts and resonances    
        for i in range(0, len(self._seq) - 1):
            pl = self._seq[i]
            # stay at planet: it is a resonance colored black
            is_reso = pl == self._seq[i + 1]
            rvt_out = rvt_outs[i]
            tof = rvt_ins[i + 1]._t - rvt_out._t
            rvt_out.plot(tof,
                units=units,
                N=4 * N,
                color="k" if is_reso else "r",
                axes=axes)
        return axes

    def eph(self, rvts, t):
        for i in range(0, len(rvts)):
            orb = rvts[i]
            if i == len(rvts) - 1 or rvts[i + 1]._t > t:
                tof = t - orb._t
                orb = orb.propagate_lagrangian(tof)
                return orb._r, orb._v
                
    def plot_distance_and_flybys(self, x, axes=None, N=1200, extension=0):
        import matplotlib.pyplot as plt
        rvt_outs, rvt_ins, rvt_pls, _, _ = self._compute_dvs(x)              
        ep = [rvt_pl._t * SEC2DAY for rvt_pl in rvt_pls]
        T = [SEC2DAY * (rvt_ins[i + 1]._t - rvt_outs[i]._t) for i in range(len(rvt_outs))]
        timeframe = np.linspace(0, sum(T) + extension, N)            
        earth = self._seq[0]
        venus = self._seq[-1]

        distances = []
        edistances = []
        vdistances = []

        for day in timeframe:
            t = x[0] + day
            pos, _ = self.eph(rvt_outs, t * DAY2SEC)
            epos, _ = earth.eph(t)
            vpos, _ = venus.eph(t)
            distances.append(np.linalg.norm(pos) / AU)
            edistances.append(np.linalg.norm(epos) / AU)
            vdistances.append(np.linalg.norm(vpos) / AU)

        fl_times = list()
        fl_distances = list()
        for pl, t in zip(self._seq, ep):
            fl_times.append(t - x[0])
            pos, _ = pl.eph(t)
            fl_distances.append(np.linalg.norm(pos) / AU)

        if axes is None:
            _, axes = plt.subplots()
        plt.ylim([0, 1.3])

        axes.plot(list(timeframe), distances, label="Solar Orbiter")
        axes.plot(list(timeframe), edistances, label="Earth")
        axes.plot(list(timeframe), vdistances, label="Venus")
        plt.scatter(fl_times, fl_distances, marker="o", color="r")
        axes.set_xlabel("Days")
        axes.set_ylabel("AU")
        axes.set_title("Distance to Sun")
        axes.legend()
        return axes


# propagate rvt_outs, rvt_ins, rvt_pls, dvs using MGA / Lambert
def _dv_mga(pl1, pl2, tof, max_revs, rvt_outs, rvt_ins, rvt_pls, dvs, lps=None):
    rvt_pl = rvt_pls[-1]  # current planet
    v_in = rvt_pl._v if rvt_ins[-1] is None else rvt_ins[-1]._v
    rvt_pl2 = rvt_planet(pl2, rvt_pl._t + tof)        
    rvt_pls.append(rvt_pl2)
    r = rvt_pl._r
    vpl = rvt_pl._v
    r2 = rvt_pl2._r
    lp = lambert_problem(r, r2, tof, rvt_pl._mu, False, max_revs)
    lp = lambert_problem_multirev_ga(v_in, lp, pl1, vpl)
    if not lps is None:
        lps.append(lp)
    v_out = lp.get_v1()[0]
    rvt_out = rvt(r, v_out, rvt_pl._t, rvt_pl._mu)
    rvt_outs.append(rvt_out)
    rvt_in = rvt(r2, lp.get_v2()[0], rvt_pl._t + tof, rvt_pl._mu)
    rvt_ins.append(rvt_in)
    vr_in = [a - b for a, b in zip(v_in, vpl)]
    vr_out = [a - b for a, b in zip(v_out, vpl)]
    dv = fb_vel(vr_in, vr_out, pl1)
    dvs.append(dv)

    
# propagate rvt_outs, rvt_ins, rvt_pls, dvs using resonances        
def compute_resonance(pl, resonances, beta, safe_distance, used_resos, reso_dts, rvt_outs,
                       rvt_ins, rvt_pls, dvs, resos=None):
    rvt_in = rvt_ins[-1]  # current spaceship
    rvt_pl = rvt_pls[-1]  # current planet
    reso = resonance(pl, rvt_in, rvt_pl, resonances)
    reso_dt, used_reso = reso.select_resonance(beta, safe_distance)
    if not resos is None:
        resos.append(reso)
    used_resos.append(used_reso)
    reso_dts.append(reso_dt)
    rvt_outs.append(reso._rvt_out)
    tof = reso.tof()
    time2 = reso._time + tof
    rvt_pl2 = rvt_planet(pl, time2)
    rvt_pls.append(rvt_pl2) 
    rvt_in2 = rvt(rvt_pl2._r, reso._rvt_out._v, time2, rvt_pl2._mu)  # its a resonance, we arrive with same v as we started
    rvt_ins.append(rvt_in2)
    dvs.append(0)  # # its a resonance, we don't need an impulse


def _angle(v1, v2):
    ca = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)  # -> cosine of the angle
    return np.arccos(np.clip(ca, -1, 1)) 


solo_mgar = _solo_mgar_udp(max_revs=2)
solo_mgar_unconstrained = _solo_mgar_udp(max_revs=2, use_constraints=False)
solo_mgar_multi_objective = _solo_mgar_udp(max_revs=2, multi_objective=True)

