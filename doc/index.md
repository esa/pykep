Welcome to pykep's documentation!
=================================

**pykep** is a coolbox developed at the [European Space Agency](https://www.esa.int)
by its [Advanced Concepts Team](https://www.esa.int/gsp/ACT/index.html). 
Its main purpose is to allow for fast prototyping of research ideas on interplanetary trajectory design.
At the library core is the implementation of efficient algorithms allowing to solve the multiple revolutions
Lambert's problem, low-thrust problems, multiple asteroid randezvous problems and more. Support for 
[JPL SPICE](https://naif.jpl.nasa.gov/naif/toolkit.html), SGP4 propagation and the 
[Heyoka](https://bluescarni.github.io/heyoka.py/index.html) Taylor integration suite is provided.

**pykep** has been used by the European Space Agency's [Advanced Concepts Team](https://www.esa.int/gsp/ACT/index.html)
during different Global Trajectory Optimization Competitions [GTOC](https://sophia.estec.esa.int/gtoc_portal), 
several research papers as well as for the optimization of preliminary mission scenarios such as in 
the M-ARGO interplanetary cubesat concept, the phase 0 study for the the Titan and Enceladus Mission (TandEM)
and for preliminary mission analysis of the HERA mission.

```{note}
If you use **pykep** in your research, please cite:

*Izzo, D. (2026). pykep: A research toolbox for space flight mechanics and interplanetary trajectory design.* [paper PDF](https://github.com/esa/pykep/blob/master/paper/paper.pdf)
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: Main

installation
changelog
bibliography
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: Tutorials

tut_basic
tut_trajopt
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: API

api
```

```{image} _static/esa_logo.png
:alt: ESA patch
:width: 55%
:align: left
``` 
```{image} _static/ACT_logo.png
:alt: ACT patch
:width: 35%
:align: right
``` 