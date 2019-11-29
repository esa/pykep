/*****************************************************************************
 *   Copyright (C) 2004-2018 The pagmo development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

// Workaround for http://mail.python.org/pipermail/new-bugs-announce/2011-March/010395.html
#ifdef _WIN32
#include <cmath>
#endif

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/self.hpp>
#include <boost/utility.hpp>

#include <keplerian_toolbox/keplerian_toolbox.hpp>
#include "../boost_python_container_conversions.h"
#include "../utils.h"

using namespace boost::python;

static inline kep_toolbox::array7D get_mismatch_wrapper(const kep_toolbox::sims_flanagan::leg &l)
{
    kep_toolbox::sims_flanagan::sc_state x;
    l.get_mismatch_con(x);
    return x.get_state();
}

static inline std::vector<double> get_throttles_con_wrapper(const kep_toolbox::sims_flanagan::leg &l)
{
    std::vector<double> ceq(l.get_throttles_size());
    l.get_throttles_con(ceq.begin(), ceq.end());
    return ceq;
}

BOOST_PYTHON_MODULE(sims_flanagan)
{
    // Disable docstring c++ signature to allow sphinx autodoc to work properly
    docstring_options doc_options;
    doc_options.disable_signatures();

    // Exposing vectors of throttles into python lists of throttles
    to_tuple_mapping<std::vector<kep_toolbox::sims_flanagan::throttle>>();
    from_python_sequence<std::vector<kep_toolbox::sims_flanagan::throttle>, variable_capacity_policy>();

    // Spacecraft class
    class_<kep_toolbox::sims_flanagan::spacecraft>(
        "spacecraft", "Contains design parameters of a NEP spacecraft",
        init<const double &, const double &, const double &>(
            "pykep.sims_flanagan.spacecraft(mass,thrust,isp)\n\n"
            "- mass: the spacecraft mass\n"
            "- thrust: the maximum thrust of the spacecraft propulsion system\n"
            "- isp: the specific impulse of the spacecarft propulsion system\n\n"
            "Examples::\n\n"
            " sc = sims_flanagan.spacecraft(4500,0.05,2500)"))
        .def("__repr__", &kep_toolbox::sims_flanagan::spacecraft::human_readable)
        .add_property("mass", &kep_toolbox::sims_flanagan::spacecraft::get_mass,
                      &kep_toolbox::sims_flanagan::spacecraft::set_mass, "The spacecraft mass\n\n"
                                                                         "Example::\n\n"
                                                                         "  mass = sc.mass"
                                                                         "  sc.mass = 2500")
        .add_property("thrust", &kep_toolbox::sims_flanagan::spacecraft::get_thrust,
                      &kep_toolbox::sims_flanagan::spacecraft::set_thrust,
                      "The spacecraft propulsion system maximum thrust\n\n"
                      "Example::\n\n"
                      "  T = sc.thrust"
                      "  sc.thrust = 0.05")
        .add_property("isp", &kep_toolbox::sims_flanagan::spacecraft::get_isp,
                      &kep_toolbox::sims_flanagan::spacecraft::set_isp,
                      "The spacecraft propulsion system specific impulse\n\n"
                      "Example::\n\n"
                      "  T = sc.isp"
                      "  sc.isp = 2000")
        .def_pickle(python_class_pickle_suite<kep_toolbox::sims_flanagan::spacecraft>())
        .def(init<>());

    // Spacecraft state class
    class_<kep_toolbox::sims_flanagan::sc_state>(
        "sc_state", "Represents the state of a spacecraft as position, velocity and mass",
        init<const kep_toolbox::array3D &, const kep_toolbox::array3D &, const double &>(
            "pykep.sims_flanagan.sc_state(r,v,m)\n\n"
            "- r: triplet containing the position vector in cartesian coordiantes\n"
            "- v: triplet containing the velocity vector in cartesian coordiantes\n"
            "- m: mass\n\n"
            "Examples::\n\n"
            " x0 = sims_flanagan.sc_state((1,0,0),(0,1,0),1000)"))
        .add_property("r", make_function(&kep_toolbox::sims_flanagan::sc_state::get_position,
                                         return_value_policy<copy_const_reference>()),
                      &kep_toolbox::sims_flanagan::sc_state::set_position,
                      "The spacecraft position in cartesian coordinates\n\n"
                      "Example::\n\n"
                      "  r = x0.r"
                      "  x0.r = (1.0,2.0,0.0)")
        .add_property("v", make_function(&kep_toolbox::sims_flanagan::sc_state::get_velocity,
                                         return_value_policy<copy_const_reference>()),
                      &kep_toolbox::sims_flanagan::sc_state::set_velocity,
                      "The spacecraft velocity in cartesian coordinates\n\n"
                      "Example::\n\n"
                      "  v = x0.v"
                      "  x0.v = (0,1.2,0)")
        .add_property("m", &kep_toolbox::sims_flanagan::sc_state::get_mass,
                      &kep_toolbox::sims_flanagan::sc_state::set_mass, "The spacecraft mass\n\n"
                                                                       "Example::\n\n"
                                                                       "  m = x0.m"
                                                                       "  x0.m = 1200")
        .def("set", &kep_toolbox::sims_flanagan::sc_state::set_state,
             "Sets the whole spacecraft state at once using a 7 dimensional tuple \n\n"
             "Example::\n\n"
             " x0.set((1,0,0,0,1,0,1000))")
        .def("get", &kep_toolbox::sims_flanagan::sc_state::get_state,
             "Gets the whole spacecraft state at once putting it into a 7 dimensional tuple \n\n"
             "Example::\n\n"
             " state = x0.get()")
        .def("__repr__", &kep_toolbox::sims_flanagan::sc_state::human_readable)
        .def_pickle(python_class_pickle_suite<kep_toolbox::sims_flanagan::sc_state>())
        .def(init<>());

    // Spacecraft throttle class
    class_<kep_toolbox::sims_flanagan::throttle>(
        "throttle", "represents one thrusting phase (segment). If the throttle is, say, (1,0,0) th thrust will be "
                    "(T,0,0), where T is the maximum value possible for a particular spacecraft",
        init<kep_toolbox::epoch, kep_toolbox::epoch, const kep_toolbox::array3D &>(
            "pykep.sims_flanagan.throttle(start,end,value)\n\n"
            "- start: starting epoch\n"
            "- end: ending epoch\n"
            "- value: triplet containing the cartesian values of the throttle\n\n"
            "Examples::\n\n"
            " s = epoch(0)"
            " e = epoch(30)"
            " t1 = sims_flanagan.throttle(s,e,(1,0,0))"))
        .add_property("start", make_function(&kep_toolbox::sims_flanagan::throttle::get_start,
                                             return_value_policy<copy_const_reference>()),
                      &kep_toolbox::sims_flanagan::throttle::set_start,
                      "The starting epoch of the throttle\n\n"
                      "Example::\n\n"
                      "  s = t1.start"
                      "  t1.start = epoch_from_string('2002-01-01 00:00:00')")
        .add_property(
            "end",
            make_function(&kep_toolbox::sims_flanagan::throttle::get_end, return_value_policy<copy_const_reference>()),
            &kep_toolbox::sims_flanagan::throttle::set_end, "The final epoch of the throttle\n\n"
                                                            "Example::\n\n"
                                                            "  e = t1.end"
                                                            "  t1.end = epoch_from_string('2002-01-23 00:00:00')")
        .add_property("value", make_function(&kep_toolbox::sims_flanagan::throttle::get_value,
                                             return_value_policy<copy_const_reference>()),
                      &kep_toolbox::sims_flanagan::throttle::set_value, "The cartesian components of the throttle\n\n"
                                                                        "Example::\n\n"
                                                                        "  components = t1.value"
                                                                        "  t1.value = (0.3,0.3,0.2)")
        .def("norm", &kep_toolbox::sims_flanagan::throttle::get_norm, "Calculates the throttle norm. If greater than 1 "
                                                                      "the resulting thrust will be greater than its "
                                                                      "maximum allowed value\n\n"
                                                                      "Example::\n\n"
                                                                      " c = t1.norm()")
        .def("__repr__", &kep_toolbox::sims_flanagan::throttle::human_readable)
        .def_pickle(python_class_pickle_suite<kep_toolbox::sims_flanagan::throttle>())
        .def(init<>());

    // Leg class
    typedef void (kep_toolbox::sims_flanagan::leg::*leg_setter)(
        const kep_toolbox::epoch &, const kep_toolbox::sims_flanagan::sc_state &, const std::vector<double> &,
        const kep_toolbox::epoch &, const kep_toolbox::sims_flanagan::sc_state &);
    typedef const std::vector<kep_toolbox::sims_flanagan::throttle> &(
        kep_toolbox::sims_flanagan::leg::*throttles_getter)();

    class_<kep_toolbox::sims_flanagan::leg>(
        "leg", "represents one low-thrust  trajectory leg in the sims-flanagan model",
        init<const kep_toolbox::epoch &, const kep_toolbox::sims_flanagan::sc_state &, const std::vector<double> &,
             const kep_toolbox::epoch &, const kep_toolbox::sims_flanagan::sc_state &,
             const kep_toolbox::sims_flanagan::spacecraft &, const double>(
            "pykep.sims_flanagan.leg(start,x0,throttles,end,xe,spacecraft,mu)\n\n"
            "- start: starting epoch\n"
            "- x0: starting sc_state\n"
            "- throttles: tuple containing the 3N cartesian components of the throttles\n"
            "- end: final epoch\n"
            "- xe: final sc_state\n"
            "- spacecarft: spacecraft\n"
            "- mu: central body gravity parameter\n\n"
            "Example::\n\n"
            " import pykep as pk"
            " start = pk.epoch(0)\n"
            " end = pk.epoch(340)\n"
            " earth = pk.planet.jpl_lp('earth')\n"
            " mars = pk.planet.jpl_lp('mars')\n"
            " sc = pk.sims_flanagan.spacecraft(4500,0.05,2500)\n"
            " r,v = earth.eph(start)\n"
            " x0 = pk.sims_flanagan.sc_state(r,v,sc.mass)\n"
            " r,v = mars.eph(start)\n"
            " xe = pk.sims_flanagan.sc_state(r, v ,sc.mass)\n"
            " mu = pk.MU_SUN\n"
            " l = pk.sims_flanagan.leg(start,x0,(1,0,0,0,0,1,0,0,0,0,1,0),end,xe,sc,mu)\n"))
        .def(init<>())
        .def("set", leg_setter(&kep_toolbox::sims_flanagan::leg::set_leg),
             "Sets leg's data, leaving the spacecraft and the central body gravitational parameter unchanged\n\n"
             "Example::\n\n"
             " l.set(start,x0,(0,0,0,1,0,0,1,0,0,0,0,0),end,xe)")
        .def("set_mu", &kep_toolbox::sims_flanagan::leg::set_mu, "Sets the leg central body gravitational parameter\n\n"
                                                                 "Example::\n\n"
                                                                 " l.set_mu(pykep.MU_SUN)")
        .def("set_spacecraft", &kep_toolbox::sims_flanagan::leg::set_spacecraft,
             "Sets the leg spacecraft\n\n"
             "Example::\n\n"
             " sc = pykep.sims_flanagan.spacecraft(4500,0.05,2500)\n"
             " l = pykep.sims_flanagan.leg()"
             " l.set_spacecraft(sc)")
        .def("get_mu", &kep_toolbox::sims_flanagan::leg::get_mu, "Gets the leg central body gravitational parameter\n\n"
                                                                 "Example::\n\n"
                                                                 " mu = l.get_mu()")
        .def("get_spacecraft", &kep_toolbox::sims_flanagan::leg::get_spacecraft,
             return_value_policy<copy_const_reference>(), "Gets the leg spacecraft\n\n"
                                                          "Example::\n\n"
                                                          " sc = pykep.sims_flanagan.spacecraft(4500,0.05,2500)\n"
                                                          " l = pykep.sims_flanagan.leg()"
                                                          " l.set_spacecraft(sc)"
                                                          " sc2 = l.get_spacecraft()")
        .def("get_throttles", throttles_getter(&kep_toolbox::sims_flanagan::leg::get_throttles),
             return_value_policy<copy_const_reference>(), "Gets the leg central body gravitational parameter\n\n"
                                                          "Example::\n\n"
                                                          " t = l.get_throttles()")
        .def("get_xi", &kep_toolbox::sims_flanagan::leg::get_x_i, return_value_policy<copy_const_reference>(),
             "Gets the initial spacecraft state\n\n"
             "Example::\n\n"
             " xi = l.get_xi()")
        .def("get_xf", &kep_toolbox::sims_flanagan::leg::get_x_f, return_value_policy<copy_const_reference>(),
             "Gets the final spacecraft state\n\n"
             "Example::\n\n"
             " xf = l.get_xf()")
        .def("get_ti", &kep_toolbox::sims_flanagan::leg::get_t_i, return_value_policy<copy_const_reference>(),
             "Gets the initial leg epoch\n\n"
             "Example::\n\n"
             " ti = l.get_ti()")
        .def("get_tf", &kep_toolbox::sims_flanagan::leg::get_t_f, return_value_policy<copy_const_reference>(),
             "Gets the final leg epoch\n\n"
             "Example::\n\n"
             " tf = l.get_tf()")
        .add_property("high_fidelity", &kep_toolbox::sims_flanagan::leg::get_high_fidelity,
                      &kep_toolbox::sims_flanagan::leg::set_high_fidelity,
                      "If True propagation is not impulsive, but continuous\n\n"
                      "Example::\n\n"
                      " l.high_fidelity = True\n")
        .def("mismatch_constraints", &get_mismatch_wrapper, "Returns a tuple containing the state mismatch of the leg "
                                                            "x,y,z,vx,vy,vz,m (needs to be all zeros for the leg to be "
                                                            "feasible)\n\n"
                                                            "Example::\n\n"
                                                            " ceq = l.mismatch_constraints()\n")
        .def("throttles_constraints", &get_throttles_con_wrapper,
             "Returns a tuple containing the throttle magnitudes minus one\n\n"
             "Example::\n\n"
             " c = l.throttles_constraints()\n")
        .def("__repr__", &kep_toolbox::sims_flanagan::leg::human_readable)
        .def_pickle(python_class_pickle_suite<kep_toolbox::sims_flanagan::leg>());

    typedef void (kep_toolbox::sims_flanagan::leg_s::*leg_s_setter)(
        const kep_toolbox::epoch &, const kep_toolbox::sims_flanagan::sc_state &, const std::vector<double> &,
        const kep_toolbox::epoch &, const kep_toolbox::sims_flanagan::sc_state &, const double &);
    class_<kep_toolbox::sims_flanagan::leg_s>(
        "leg_s", "Represents an interplanetary leg (using the Sundmann variable)",
        init<unsigned, double , double , optional<double >>(
            "pykep.sims_flanagan.leg_s(n_seg, c, alpha, tol = -10)\n\n"
            "- n_seg: number of segments\n"
            "- c: constant in the SUndmann transformation dt = cr^(alpha ds)\n"
            "- alpha: exponent in the SUndmann transformation dt = cr^(alpha ds)\n"
            "- tol: log 10 tolerance set in the Taylor integration of the leg\n\n"
            "Example::\n\n"
            " from pykep import *\n"
            " l = sims_flanagan.leg_s(25,1.0/AU,1.0)\n"))
        .def(init<>())
        .def("set", leg_s_setter(&kep_toolbox::sims_flanagan::leg_s::set_leg),
             "Sets leg's data\n\n"
             "Example::\n\n"
             " import pykep as pk\n"
             " l = pk.sims_flanagan.leg_s(4,1.0/pk.AU,1.0)\n"
             " start = pk.epoch(0)\n"
             " end = pk.epoch(340)\n"
             " earth = pk.planet.jpl_lp('earth')\n"
             " mars = pk.planet.jpl_lp('mars')\n"
             " sc = pk.sims_flanagan.spacecraft(4500,0.05,2500)\n"
             " r,v = earth.eph(start)\n"
             " x0 = pk.sims_flanagan.sc_state(r,v,sc.mass)\n"
             " r,v = mars.eph(end)\n"
             " xe = pk.sims_flanagan.sc_state(r, v ,sc.mass)\n"
             " ds = 1.2 * 365.25 * pk.DAY2SEC\n"
             " l.set_mu(pk.MU_SUN)\n"
             " l.set_spacecraft(sc)\n"
             " l.set(start, x0, (1,0,0,0,0,1,0,0,0,0,1,0), end, xe, ds)\n")
        .def("set_mu", &kep_toolbox::sims_flanagan::leg_s::set_mu,
             "Sets the leg central body gravitational parameter\n\n"
             "Example::\n\n"
             " l.set_mu(pykep.MU_SUN)")
        .def("set_spacecraft", &kep_toolbox::sims_flanagan::leg_s::set_sc,
             "Sets the leg spacecraft\n\n"
             "Example::\n\n"
             " sc = pykep.sims_flanagan.spacecraft(4500,0.05,2500)\n"
             " l = pykep.sims_flanagan.leg()"
             " l.set_spacecraft(sc)")
        .def("get_mu", &kep_toolbox::sims_flanagan::leg_s::get_mu,
             "Gets the leg central body gravitational parameter\n\n"
             "Example::\n\n"
             " mu = l.get_mu()")
        .def("get_spacecraft", &kep_toolbox::sims_flanagan::leg_s::get_spacecraft,
             return_value_policy<copy_const_reference>(), "Gets the leg spacecraft\n\n"
                                                          "Example::\n\n"
                                                          " sc = l.get_spacecraft()")
        .def("get_xi", &kep_toolbox::sims_flanagan::leg_s::get_xi, return_value_policy<copy_const_reference>(),
             "Gets the initial spacecraft state\n\n"
             "Example::\n\n"
             " xi = l.get_xi()")
        .def("get_xf", &kep_toolbox::sims_flanagan::leg_s::get_xf, return_value_policy<copy_const_reference>(),
             "Gets the final spacecraft state\n\n"
             "Example::\n\n"
             " xf = l.get_xf()")
        .def("get_ti", &kep_toolbox::sims_flanagan::leg_s::get_ti, return_value_policy<copy_const_reference>(),
             "Gets the initial leg epoch\n\n"
             "Example::\n\n"
             " ti = l.get_ti()")
        .def("get_tf", &kep_toolbox::sims_flanagan::leg_s::get_tf, return_value_policy<copy_const_reference>(),
             "Gets the final leg epoch\n\n"
             "Example::\n\n"
             " tf = l.get_tf()")
        .def("get_throttles", &kep_toolbox::sims_flanagan::leg_s::get_throttles,
             return_value_policy<copy_const_reference>(), "Reurns a tuple containing the leg's throttles\n\n"
                                                          "Example::\n\n"
                                                          " th = l.get_throttles()")
        .def("mismatch_constraints", &kep_toolbox::sims_flanagan::leg_s::compute_mismatch_con,
             return_value_policy<copy_const_reference>(), "Returns an 8-dim tuple containing the state mismatch of the "
                                                          "leg (x,y,z,vx,vy,vz,m,t) (needs to be all zeros for the leg "
                                                          "to be feasible)\n\n"
                                                          "Example::\n\n"
                                                          " ceq = l.mismatch_constraints()\n")
        .def("throttles_constraints", &kep_toolbox::sims_flanagan::leg_s::compute_throttles_con,
             return_value_policy<copy_const_reference>(), "Returns a nseg-dim tuple containing the throttle magnitudes "
                                                          "minus one (needs to be all negative for the leg to be "
                                                          "feasible)\n\n"
                                                          "Example::\n\n"
                                                          " c = l.throttles_constraints()\n")
        .def(
            "states", &kep_toolbox::sims_flanagan::leg_s::compute_states, return_value_policy<copy_const_reference>(),
            "Returns a nseg+2 tuple containing 11-dim tuples of the spacecraft state: (t,x,y,z,vx,vy,vz,m,ux,uy,uz)\n\n"
            "Example::\n\n"
            " states = l.states()\n")
        .def("__repr__", &kep_toolbox::sims_flanagan::leg_s::human_readable)
        .def_pickle(python_class_pickle_suite<kep_toolbox::sims_flanagan::leg_s>());
}
