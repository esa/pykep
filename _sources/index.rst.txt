.. pykep documentation master file, created by
   sphinx-quickstart on Thu Nov  4 12:34:23 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==========================================
Welcome to pykep
==========================================

.. image:: images/traj1.gif
.. image:: images/traj2.gif
.. image:: images/traj4.gif
.. image:: images/traj5.gif
.. image:: images/traj6.gif
.. image:: images/github.png
   :target: https://github.com/esa/pykep
   :scale: 50%


**pykep** is a scientific library developed at the `European Space Agency <https://www.esa.int/gsp/ACT/index.html>`_
to provide basic tools for astrodynamics research. Algoritmic efficiency is
a main focus of the library, written in C++ and exposed to Python. At the library core
is the implementation of an efficient solver for the multiple revolutions Lambert's problem, objects representing 
direct (Sims-Flanagan), indirect (Pontryagin) and hybrid methods to represent low-thrust optimization problems, 
efficient keplerian propagators, Taylor-integrators, a SGP4 propagator, TLE and SATCAT support, 
`JPL SPICE <https://naif.jpl.nasa.gov/naif/toolkit.html>`_ code interface and more.

**pykep** has been used by the European Space Agency's `Advanced Concepts Team <https://www.esa.int/gsp/ACT/index.html>`_ 
during different Global Trajectory Optimization Competitions `GTOC <http://sophia.estec.esa.int/gtoc_portal>`_, 
several research papers as well as for the optimization of preliminary mission scenarios such as in 
the M-ARGO interplanetary cubesat concept or the phase 0 study for the The Titan and Enceladus Mission (TandEM).

Please quote **pykep** code using the information in Zenodo.

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2575462.svg
   :target: https://doi.org/10.5281/zenodo.2575462
   :align: left

---------------------------------------------------------------------------------

Content
=================================
.. toctree::
  :maxdepth: 2

  downloading
  system
  installation
  documentation/index 
  examples/index 
  credits
