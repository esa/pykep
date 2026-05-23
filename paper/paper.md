---
title: 'pykep: A coolbox for space flight mechanics'
tags:
  - Python
  - astrodynamics
  - interplanetary trajectories
  - preliminary design
  - optimization
authors:
  - name: Dario Izzo
    orcid: 0000-0002-9846-8423
    affiliation: 1
affiliations:
  - name: European Space Agency's Advanced Concepts Team, The Netherlands
    index: 1
date: 23 May 2026
bibliography: paper.bib
---

# Summary

`pykep` is a Python coolbox covering a broad range of space flight mechanics problems, with its
current emphasis on the preliminary analysis of interplanetary spacecraft trajectories. It provides
the mathematical building blocks for orbit propagation, orbital element conversions, Lambert arc
solvers, gravity-assist flyby modelling, and both direct and indirect low-thrust transcriptions.
Given a mission scenario—a spacecraft departing from one body in the solar system and arriving at
another, possibly via intermediate gravity-assist flybys and with low-thrust propulsion phases—`pykep`
allows researchers to compute transfer costs (in terms of propellant mass and time of flight), set
up the corresponding optimization problems, and interface with state-of-the-art numerical optimizers.
The library's C++ core can evaluate thousands of candidate trajectories per second, enabling the
large-scale search over launch dates, flyby sequences, and thrust profiles that is essential during
the earliest stages of mission design. `pykep` is not an operational flight dynamics tool; it is a
research instrument designed to let scientists and engineers prototype new ideas and push the
frontier of what is feasible across the full landscape of space flight mechanics problems.

# Statement of Need

The preliminary design of an interplanetary mission is fundamentally an exploration problem: analysts
must survey a vast combinatorial space of departure dates, planetary flyby sequences, and propulsion
strategies to identify the most promising candidates before committing to costly detailed analysis.
Certified tools for operational flight dynamics prioritize physical fidelity and validated models over
the speed and composability needed to iterate over millions of trajectory options rapidly.

`pykep` was developed at the European Space Agency's Advanced Concepts Team (ACT) specifically to
fill this gap. Its target audience is researchers and scientists working on astrodynamics, orbital
mechanics, space flight mechanics, and related computational fields who need a flexible, Pythonic
toolkit to prototype algorithms, benchmark new methods, and conduct large-scale parametric studies.
While the library's current focus lies on the preliminary analysis of interplanetary trajectories,
its scope is deliberately broader: the same core primitives—propagators, element conversions, and
trajectory leg models—are equally applicable to Earth-orbit problems, lunar transfers, and other
space flight mechanics domains. The library is explicitly *not* intended for operational mission
design or for use in certified flight dynamics pipelines; it is a research tool for pushing the
frontier of what is computable at the earliest stages of mission analysis.

`pykep` deliberately exposes its mathematical primitives—Lambert arc solvers, low-thrust trajectory
legs, orbital element conversions—as first-class objects that can be composed freely and embedded in
optimization loops. By integrating natively with `pagmo` [@pagmo] (a parallel global optimization
library) and `Heyoka` [@biscaniheyoka1] (a Taylor integration suite), `pykep` provides a full stack
from ephemeris evaluation to global trajectory search within a single Python environment—a
coolbox where each tool is a well-defined, composable primitive ready to be assembled into novel
research workflows. The many algorithms available in the library are largely the result of original
research carried out at the ACT over more than a decade, complemented by implementations of
best-in-class methods from the broader astrodynamics literature.

# State of the Field

Several open-source tools address related aspects of spacecraft trajectory design, and it is
important to articulate how `pykep` relates to and differs from each.

**Orekit** [@orekit] is a Java-based orbital mechanics library (with a Python facade) targeted at
operational applications. It provides high-fidelity force models and validated propagators, making
it the tool of choice when a precise, certified answer to a specific trajectory question is required.
Its design is not optimized for the rapid, repeated evaluations over large search spaces that
characterize preliminary mission design.

**GMAT** [@gmat] is NASA's General Mission Analysis Tool, a comprehensive high-fidelity simulator
for a wide range of mission types. It is primarily a GUI- and script-driven application oriented
toward operational analysis rather than algorithmic exploration and programmatic access.

**GPOPS-II** [@gpops2] is a MATLAB-based optimal control software using hp-adaptive Gaussian
quadrature collocation. It excels at locally converging a single, complex continuous trajectory but
is proprietary and not designed for global search over mission scenarios or large parameter sweeps.

**poliastro** [@poliastro] is a Python library for general orbital mechanics. It shares some
features with `pykep` (orbit propagation, Lambert solvers) but is oriented toward education and
general orbital analysis rather than the specialized requirements of preliminary interplanetary
mission design, specifically multi-gravity-assist sequences, low-thrust transcriptions, and tight
integration with global optimization algorithms.

`pykep` occupies a distinct niche by combining, in a single research-oriented Python package: fast
multi-revolution Lambert solvers [@izzolambert]; both direct
(Sims-Flanagan [@sims], zero-order-hold [@izzo2026zoh]) and indirect (Pontryagin-based) low-thrust
transcriptions; native coupling with `pagmo` for global and multi-objective optimization;
Taylor-based high-accuracy propagation via `Heyoka` [@biscaniheyoka1; @biscaniheyoka2]; and a
growing benchmark suite (`gym`) enabling reproducible comparisons of trajectory optimization
algorithms. No existing open-source package provides this combination in a research-grade Python
library explicitly designed to support the development of novel algorithms.

# Software Design

`pykep 3` represents a complete redesign of the library from its earlier Python 2 origins. The
computational core is implemented in modern C++23 and exposed to Python via `pybind11`. This
architecture ensures that the innermost evaluation loops—Lambert arcs, Sims-Flanagan mismatch
constraints [@sims], zero-order-hold propagation [@izzo2026zoh], orbital element conversions
[@equinoctialelems]—run at native speed while remaining fully accessible from Python.

A central abstraction is the *user-defined planet-like object* (`udpla`): any object that provides
position and velocity as a function of epoch can serve as a trajectory endpoint or flyby body,
regardless of whether it is backed by SPICE kernels, the SGP4 model, an analytic Keplerian
approximation, or a user-supplied function. This makes all higher-level trajectory solvers
ephemeris-agnostic and straightforward to test with analytical models before switching to full SPICE
data.

Trajectory optimization problems are exposed as *user-defined problem* (`udp`) objects compatible
with `pagmo` [@pagmo], so that any optimizer in the `pagmo` ecosystem—differential evolution,
particle swarm, self-adaptive evolution strategies, and branch-and-bound—can be applied without
modification. The design separates the mathematical primitives (C++) from the optimization problem
definitions (Python/C++) and auxiliary tools such as visualization and low-thrust approximations
[@approximations] (pure Python). This layering keeps the library maintainable and makes it
straightforward for researchers to contribute new trajectory models or problem formulations. A set
of benchmark problems (`gym`), modeled on challenges from the Global Trajectory Optimization
Competitions [@izzo2016designing], is included to support reproducible algorithm development and
direct comparison of competing approaches.

# Research Impact Statement

`pykep` and its predecessors have been the primary computational tool of the ACT for over a decade
of competition-level and research-grade interplanetary trajectory work. The library was instrumental
in the ACT's participation in multiple editions of the Global Trajectory Optimization Competition
(GTOC) [@izzo2016designing; @izzo2013search; @gtoc12]—a series of open benchmarks widely regarded
as among the most demanding trajectory optimization problems in the field, and a venue where the
library has consistently enabled top-tier results against leading astrodynamics teams worldwide.

Beyond competition, `pykep` has underpinned preliminary mission analysis for several ESA concepts,
including the Hera asteroid deflection mission, the M-ARGO interplanetary CubeSat, and early
feasibility studies for the Titan and Enceladus Mission (TandEM) and the Laplace mission to the
outer solar system. It has supported peer-reviewed research on novel indirect methods via surrogate
primer vectors [@beauregard], comparative studies of machine-learning and astrodynamical approaches
to low-thrust transfers [@acciarini2024computing], approximation methods for asteroid-belt transfer
cost estimation [@approximations], and high-order Taylor methods for astrodynamics
[@biscaniheyoka1; @biscaniheyoka2].

The library is distributed under an open-source license, available on PyPI and conda-forge, and
is actively maintained. Its documentation includes tutorials covering basic propagation, Lambert
problems, multi-gravity-assist trajectory optimization, and low-thrust mission design, lowering the
barrier for adoption by the broader research community.

# AI Usage Disclosure

AI-assisted writing tools (GitHub Copilot) were used in an assistive capacity during the drafting
of this paper. All technical content, bibliographic references, and scientific claims were reviewed
and verified by the authors.

# Acknowledgements

The authors thank the many ACT researchers and interns who have used and improved `pykep` over the
years. Special thanks are due to Francesco Biscani for the co-development of `pagmo` and `Heyoka`, and
to the open-source communities behind `pybind11` and `spiceypy`. Development has been supported by
the European Space Agency through the Advanced Concepts Team's internal research programme.

# References

