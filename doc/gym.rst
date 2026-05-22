.. _gym:

Trajectory Optimization Gym
###########################

A number of interplanetary trajectory optimization benchmarks are provided in
``pykep`` as User Defined Problems (UDPs), compatible with
`pygmo <https://esa.github.io/pygmo2/>`_ :cite:p:`pagmo`.

All gym problems are instantiated when importing ``pykep`` and can therefore be
used directly. The full collection is referred to as the ``pykep`` gym and is
intended as a consistent benchmark set for both evolutionary and gradient-based
trajectory optimization methods.

.. currentmodule:: pykep.trajopt.gym
 
MGA Problems
************

The following problems are based on Multiple Gravity Assist (MGA)
transcriptions.

.. autoattribute:: pykep.trajopt.gym.cassini1

MGA benchmark inspired by the historical Cassini-Huygens interplanetary
transfer to Saturn. The trajectory combines inner-planet gravity assists to
increase heliocentric energy before the final outer-planet leg.
The fly-by sequence is ``E-VVEJ-S`` and reproduces the classic mission-design
logic of chaining resonant planetary encounters to reduce propellant demand.
The objective is the total :math:`\Delta V`, including launch excess
:math:`\Delta V`, intermediate maneuvers, and final Saturn orbit insertion
(:math:`r_p = 108950\ \mathrm{km}`, :math:`e = 0.98`).

.. autoattribute:: pykep.trajopt.gym.cassini1_a

Same physical problem as :class:`~pykep.trajopt.gym.cassini1`, but using
``alpha`` time-of-flight encoding.
This encoding is useful when per-leg time-of-flight bounds are not known
a priori. It introduces an additional decision variable for total
time-of-flight, making it convenient for multi-objective formulations where
total mission duration is an explicit objective.

.. autoattribute:: pykep.trajopt.gym.cassini1_n

Same physical problem as :class:`~pykep.trajopt.gym.cassini1`, but using
``eta`` time-of-flight encoding.
As with the alpha variant, this parametrization is useful when direct,
per-leg time-of-flight bounds are not available a priori.

MGA-1DSM Problems
*****************

The following benchmarks are MGA problems with one Deep Space Maneuver (DSM)
per leg.

.. autoattribute:: pykep.trajopt.gym.cassini2

MGA-1DSM benchmark inspired by the Cassini transfer to Saturn.
Compared to :class:`~pykep.trajopt.gym.cassini1`, one DSM is allowed on each
leg, and launch :math:`\Delta V` is modeled as a bounded quantity rather than
being directly minimized. This shifts more optimization freedom to in-flight
control decisions, making the benchmark representative of preliminary design
studies where launch capability is preallocated and cruise shaping is critical.

.. autoattribute:: pykep.trajopt.gym.rosetta

MGA-1DSM benchmark inspired by Rosetta's transfer to comet
67P/Churyumov-Gerasimenko.
The fly-by sequence is ``E-EMEE-comet``.
The objective is rendezvous (position and velocity matching at arrival) while
minimizing total maneuver :math:`\Delta V`; launch :math:`\Delta V` is bounded.
The instance captures a key challenge of small-body rendezvous missions:
constructing a long, weakly forced transfer that arrives with low relative
velocity after multiple inner-Solar-System gravity assists.

.. autoattribute:: pykep.trajopt.gym.eve_mga1dsm

Compact MGA-1DSM benchmark for algorithm testing.
The planetary sequence is ``E-V-E``.
The objective is Earth rendezvous (position and velocity matching) with minimum
total maneuver :math:`\Delta V`; launch :math:`\Delta V` is bounded.
Despite its compact size, this benchmark preserves the main coupling patterns
of larger MGA-1DSM problems and is therefore useful for rapid solver
prototyping, tuning, and regression testing.

.. autoattribute:: pykep.trajopt.gym.eve_mga1dsm_a

Same physical problem as :class:`~pykep.trajopt.gym.eve_mga1dsm`, with
``alpha`` time-of-flight encoding (see
:class:`~pykep.trajopt.gym.cassini1_a`).

.. autoattribute:: pykep.trajopt.gym.eve_mga1dsm_n

Same physical problem as :class:`~pykep.trajopt.gym.eve_mga1dsm`, with
``eta`` time-of-flight encoding (see
:class:`~pykep.trajopt.gym.cassini1_n`).

.. autoattribute:: pykep.trajopt.gym.juice

MGA-1DSM benchmark inspired by ESA's JUICE transfer to Jupiter.
The sequence is ``E-EVEME-J`` and includes an Ariane 5 launcher model from
Kourou. The instance reflects the mission-design context of delivering a large
science payload to the Jovian system through an energy-constrained,
multi-assist trajectory with strict launch and cruise trade-offs.

.. autoattribute:: pykep.trajopt.gym.juice_mo

Multi-objective variant of :class:`~pykep.trajopt.gym.juice`.
It uses ``alpha`` encoding for time-of-flight variables (see
:class:`~pykep.trajopt.gym.cassini1_a`) and includes total mission
time-of-flight as a second objective. This formulation is aimed at constructing
trade-off fronts between propellant usage and transfer duration in the same
mission context.

.. autoattribute:: pykep.trajopt.gym.messenger

MGA-1DSM benchmark inspired by the MESSENGER transfer to Mercury.
The fly-by sequence is ``E-VVMeMeMe-Me`` (with the first Earth correction
fly-by omitted because no launcher model is used).
This benchmark is intentionally difficult: resonant Mercury fly-bys and
multiple-revolution opportunities create a rugged optimization landscape and
high sensitivity to time-of-flight choices. It is a representative stress test
for global optimization methods in chemically propelled, gravity-assist-heavy
inner-planet missions.

Multiple-Impulse Problems
*************************

These benchmarks represent Earth-Moon transfer problems with a fixed number of
impulses.
They provide compact testbeds for studying how optimization difficulty scales
with control dimensionality as the number of allowed impulses increases.

.. autoattribute:: pykep.trajopt.gym.em3imp

.. autoattribute:: pykep.trajopt.gym.em5imp

.. autoattribute:: pykep.trajopt.gym.em7imp

TOPS Problems
*************

The TOPS (Trajectory Optimisation Problems in Space) benchmarks
:cite:p:`izzo2026zoh` are a collection of low-thrust trajectory optimization problems
specified in JSON format.

In ``pykep`` these datasets are available as dictionaries under:
:data:`~pykep.trajopt.gym.tops_twobody_json`,
:data:`~pykep.trajopt.gym.tops_mee_json`,
:data:`~pykep.trajopt.gym.tops_ss_json`, and
:data:`~pykep.trajopt.gym.tops_cr3bp_json`, corresponding to two-body
Cartesian, two-body MEE, solar-sail, and CR3BP dynamics.

NLP transcriptions of selected TOPS instances are also provided as UDP classes
in ``pykep.trajopt.gym``.

All formulations use zero-order-hold control over a user-defined number of
segments and a forward-backward shooting scheme (:cite:p:`izzo2026zoh`). 
This induces nonlinear mismatch constraints; non-solar-sail variants
also include throttle constraints. The resulting NLP size and difficulty
can be tuned through ``nseg`` and the time-grid encoding
(for example ``uniform`` or ``softmax``).

Fixed Boundaries
----------------

A first set of TOPS instances as NLPs is provided with fixed initial/final
states and, in most cases, fixed time-of-flight.
These benchmarks are suitable for controlled algorithm comparisons and as an
entry point to the ``pykep`` gym.
Non-solar-sail variants are built on :class:`~pykep.trajopt.zoh_point2point`,
while the solar-sail variant is built on
:class:`~pykep.trajopt.zoh_ss_point2point`.

.. autoclass:: pykep.trajopt.gym.tops_twobody

.. autoclass:: pykep.trajopt.gym.tops_mee

.. autoclass:: pykep.trajopt.gym.tops_ss

.. autoclass:: pykep.trajopt.gym.tops_cr3bp

Moving Boundaries
-----------------

Moving-boundary NLP formulations extend the fixed-boundary setting by allowing
departure and arrival epochs to vary, with endpoint states tied to
:class:`~pykep.planet` ephemerides and optional relative-velocity constraints.
This setup is closer to preliminary mission design use cases and introduces a
harder optimization problem because boundary epochs become decision variables and
two constraints are added to control the relative velocity at both endpoints.

In terms of base transcriptions, non-solar-sail variants are formulated on
:class:`~pykep.trajopt.zoh_pl2pl`, while solar-sail variants use
:class:`~pykep.trajopt.zoh_ss_pl2pl`.

.. note::

	At the moment, only the fixed-boundary TOPS classes are exposed under
	:mod:`pykep.trajopt.gym`.