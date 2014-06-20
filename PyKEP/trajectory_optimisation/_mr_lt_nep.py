from PyGMO.problem._base import base_problem
import PyKEP
import gtoc7

from copy import deepcopy


class mr_lt_nep(base_problem):

    """
    This class represents, as a global optimization problem (linearly constrained, 
    high diemensional), a Multiple Randezvous trajectory of a low-thrust spacecraft equipped
    with a nuclear electric propulsion engine.

    The decision vector (chromosome) is:
    [t1, tof1, rest1, m_f1] + [throttles1] + 
    [t2, tof2, rest2, m_f2] + [throttles2] + ....

    """

    def __init__(
            self,
            mass=[800, 2000],
            Tmax=0.3,
            Isp=3000,
            nseg=5,
            seq=[PyKEP.planet_gtoc7(3413),PyKEP.planet_gtoc7(234), PyKEP.planet_gtoc7(11432)]
            tof=[1, 365.25 * 3], 
            rest=[30., 365.25], 
            max_total_time=365.25 * 6, 
            t0=[13000, 13200],
            objective='mass',
            **kwargs
    ):

        n = len(seq) - 1
        dim = (4 + nseg * 3) * n + 1
        dim_eq = 7 * n
        dim_ineq = n * nseg + n
        
        super(mr_lt_nep, self).__init__(dim, 0, 1, 
                                            dim_eq + dim_ineq, # constraint dimension
                                            dim_ineq, # inequality constraints
                                            1e-5) # constrain tolerance
        
        self.__seq_ids = seq
        self.__seq = seq
        self.__num_legs = n
        self.__nseg = nseg
        self.__dim_leg = 4 + nseg * 3
        self.__start_mass = mass[1]
        self.__max_total_time = max_total_time
        
        if not objective in ['mass', 'time']:
            raise ValueError("Error in defining the objective. Was it one of mass or time?")
        self.__objective = objective
        
        # set leg bounds
        lb_leg = [t0[0],                  tof[0], rest[0], mass[0]] + [-1] * nseg * 3
        ub_leg = [t0[1]+max_total_time*n, tof[1], rest[1], mass[1]] + [1] * nseg * 3

        # set n leg bounds
        lb = lb_leg * n
        ub = [t0[1], tof[1], rest[1], mass[1]] + [1] * nseg * 3 + ub_leg * (n-1)

        # set total time bounds
        lb += [1.]
        ub += [self.__max_total_time]

        self.set_bounds(lb, ub)

    def _objfun_impl(self, x):
        final_mass = x[-1-self.__dim_leg + 3]
        tof = x[-1-self.__dim_leg - x[0]
        if self.__objective == 'mass':
            return ( -final_mass, )
        elif self.__objective == 'time':
            return (  tof, )


    def _compute_constraints_impl(self, x_full):
        sc_mass = self.__start_mass
        eqs = []
        ineqs = []

        for i in range(self.__num_legs):
            x = x_full[i*self.__dim_leg:(i+1)*self.__dim_leg]
            
            start = PyKEP.epoch(x[0])
            end = PyKEP.epoch(x[0] + x[1])
            
            # Computing starting spaceraft state
            r, v = self.__seq[i].eph(start)
            x0 = PyKEP.sims_flanagan.sc_state(r, v, sc_mass)
            
            # Computing ending spaceraft state
            r, v = self.__seq[i+1].eph(end)
            xe = PyKEP.sims_flanagan.sc_state(r, v, x[3])
            
            # Building the SF leg
            leg = PyKEP.sims_flanagan.leg()
            leg.high_fidelity = True
            leg.set_mu(PyKEP.MU_SUN)
            leg.set_spacecraft(PyKEP.sims_flanagan.spacecraft(sc_mass, .3, 3000.))
            leg.set(start, x0, x[-3 * self.__nseg:], end, xe)
                
            # Setting all constraints
            eqs.extend(leg.mismatch_constraints())
            ineqs.extend(leg.throttles_constraints())

            eqs[-7] /= PyKEP.AU
            eqs[-6] /= PyKEP.AU
            eqs[-5] /= PyKEP.AU
            eqs[-4] /= PyKEP.EARTH_VELOCITY
            eqs[-3] /= PyKEP.EARTH_VELOCITY
            eqs[-2] /= PyKEP.EARTH_VELOCITY
            eqs[-1] /= self.__start_mass

            sc_mass = x[3] # update mass to final mass of leg
            
            if i < self.__num_legs - 1:
                x_next = x_full[(i+1)*self.__dim_leg:(i+2)*self.__dim_leg]
                time_ineq = x[0] + x[1] + x[2] - x_next[0]
                ineqs.append(time_ineq/365.25)
            else:
                final_time_ineq = x[0] + x[1] + x[2] - x_full[0] - x_full[-1] # <- total time
                ineqs.append(final_time_ineq/365.25)

        retval = eqs + ineqs
        return retval


    def initial_guess(self, x):
        """Creates an initial guess based on the compact trajectory format."""

        # set last resting time to 0
        x = deepcopy(x)
        x[-1] = list(x[-1])
        x[-1][0] = sum(x[-2][:2])
        x[-1] = tuple(x[-1])

        n = len(x)-1
        initial_guess = []
        for i in range(n):
            initial_guess.extend([x[i][0], x[i][1], 30., x[i+1][2]] + [0] * self.__nseg * 3)
        initial_guess.append(x[-1][0] - x[0][0])
        return initial_guess


    def to_compact(self, x):
        """Returns the compact trajectory format of a chromosome."""
        epochs = list(x[::self.__dim_leg])[:-1] # remove last because that's total time
        tofs =  list(x[1::self.__dim_leg])
        #resting_times = x[2::self.__dim_leg]
        masses = list(x[3::self.__dim_leg])
        epochs.append(epochs[-1] + tofs[-1] + 30.) # compute last epoch
        tofs.append(0.) # set last tof
        masses.insert(0, self.__start_mass) # add starting mass
        return zip(epochs, tofs, masses, self.__seq_ids)


    def resting_times(self, x):
        return list(x[2::self.__dim_leg])


def optimize_full_daughter(x, start_time_slack=30., without_resting=False, 
                           nseg=5, major_iter=1000, opt_tol=1e-03, feas_tol=1e-09,
                           mass=[1., 2000.], verbose=True, **kwargs):
    """Optimizes a full daughter sequence."""
    import PyGMO
    max_total_time = 365.25 * 6 - 30.
    if without_resting:
        max_total_time += 60.

    seq = [t[-1] for t in x]
    t0 = x[0][0]

    algo = PyGMO.algorithm.snopt(
        screen_output=True,
        major_iter=major_iter,
        opt_tol=opt_tol,
        feas_tol=feas_tol)
    prob = full_daughter(t0=[t0-start_time_slack, t0+start_time_slack],
                         seq=seq, 
                         nseg=nseg, 
                         mass=mass, 
                         #max_total_time=max_total_time, # we don't care if it can score first and last asteroid
                         **kwargs)
    pop = PyGMO.population(prob, 0)
    pop.push_back(prob.initial_guess(x))
    pop = algo.evolve(pop)
    if prob.feasibility_x(pop.champion.x):
        print "FEASIBILE!!!"
    else:
        print "INFEASIBLE :("
    return prob.to_compact(pop.champion.x), prob.feasibility_x(pop.champion.x)


if __name__ == "__main__":
    import gtoc7
    # mountain first probe, not optimized by hippo
    x = [(11966.174890358025, 219.26961918210563, 2000.0, 14253),
 (12216.44450954013, 165.3323384041047, 1806.8158331604805, 12589),
 (12411.776847944235, 141.15577960424727, 1661.1522881358278, 8476),
 (12583.932627548482, 178.60902451557592, 1536.789132925773, 9655),
 (12793.541652064057, 111.33031949556417, 1379.4283653016366, 326),
 (12935.871971559622, 105.42961415434141, 1281.342478637386, 1888),
 (13071.301585713964, 146.75609663509258, 1188.4553187715778, 9649),
 (13248.057682349056, 103.63751082607308, 1059.1580893475061, 1039),
 (13382.695193175128, 104.30789736711017, 980.1205785006173, 699),
 (13517.003090542239, 79.46649008620508, 888.2216901694909, 11002),
 (13626.469580628444, 98.67099816776766, 819.7771999999624, 10549),
 (13755.14057879621, 64.49967763581671, 732.844616072239, 8487),
 (13850.640256432027, 97.739531883387, 676.0181549693853, 8514),
 (13978.379788315415, 0.0, 589.9062252610506, 3734)]

    new_x = gtoc7.optimize_full_daughter(x, opt_tol=1e-03, feas_tol=1e-09)
    print "x=%s" % gtoc7.utils.format_trajectory(new_x)
