---
title: 'PyKEP 3: A coolbox for interplanetary trajectory design'
tags:
  - Python
  - astrodynamics
  - interplanetary trajectories
  - preliminary design
  - optimization
authors:
  - name: Dario Izzo
    orcid: 0000-0002-9846-8423
    equal-contrib: false
    affiliation: "1," # (Multiple affiliations must be quoted)
affiliations:
 - name: European Space Agency's Advanced Concepts Team, The Netherlands
   index: 1
 - name: European Space Agency's Advanced Concepts Team, The Netherlands
   index: 2
date: 28 March 2025
bibliography: paper.bib
---

# Summary

`PyKEP 3` is a Python toolbox developed at the [European Space Agency](https://www.esa.int) by the 
[Advanced Concepts Team](https://www.esa.int/act) to perform
quick analysis of interplanetary trajectory design problems. It is designed to be used by researchers
and engineers to prototype and test new ideas in the field of astrodynamics. The library provides
efficient implementations of algorithms for solving the multiple revolutions Lambert's problem, low-thrust
problems, multiple asteroid rendezvous problems, and more. It also provides support for [JPL SPICE](https://naif.jpl.nasa.gov/naif/toolkit.html),
SGP4 propagation, and the [Heyoka](https://bluescarni.github.io/heyoka.py/index.html) Taylor integration suite.

# Introduction

`PyKEP 3` is the third iteration of the PyKEP library, which has been used by the [Advanced Concepts Team](https://www.esa.int/act)
at the European Space Agency for over a decade. The library was originally developed to support the participation
of non astrodynamics experts in the [Global Trajectory Optimization Competition](http://www.esa.int/gtoc) and has since
grown into a generic tool for the preliminary analysis of interplanetary trajectory design problems. 
Previous versions of the code were used internally to perform some preliminary analysis for the design of the 
[Hera mission](https://www.heramission.space), the M-ARGO interplanetary cubesat concept,
the preliminary analysis of the the initial concepts for the Titan and Enceladus Mission (TandEM) as well as the Laplace mission.
At the core of the library are the necessary building blocks to perform a preliminary design of interplanetary trajectory design problems,
with the aim to provide the general scientific community with a tool to prototype and test new ideas in the field of
astrodynamics focussing on novel ideas rather than on on the implementation of known algorithms. Most
techniques implemented in the library are largely based on original research performed at the ACT, as well as on the state-of-the art
algorithms available from the literature. A Lambert solver as well as low-thrust optimization algorithms (both based on direct and indirect methods)
are available together with a number of approximate methods to estimate transfer costs. 
The library interfaces to a number of open source projects providing relevant functionalities including [spicepy](https://spiceypy.readthedocs.io/en/stable/) for interfacing to JPL SPICE,
[SGP4](https://pypi.org/project/sgp4/) to get standard Low Earth Orbit orbits, 
[Heyoka](https://bluescarni.github.io/heyoka.py/) to perform efficient Taylor integration and [pygmo](https://esa.github.io/pygmo2/) to interface to state-of-the art optimization algorithms.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from .... and support from ...

# References

