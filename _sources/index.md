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

```{toctree}
:maxdepth: 1
:caption: Main

api
bibliography
```

```{toctree}
:maxdepth: 1
:caption: Tutorials

tut_basic
tut_trajopt
```

```{image} _static/esa_logo.png
:alt: ESA patch
:width: 400px
:align: center
``` 