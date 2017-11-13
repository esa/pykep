.. pykep documentation master file, created by
   sphinx-quickstart on Thu Nov  4 12:34:23 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==========================================
Welcome
==========================================

.. image:: images/traj1.gif
.. image:: images/traj2.gif
.. image:: images/traj3.gif
.. image:: images/traj4.gif
.. image:: images/traj5.gif
.. image:: images/traj6.gif
.. image:: images/logo_pykep.png


pykep is a scientific library providing basic tools for astrodynamics research. Algoritmic efficiency is
a main focus of the library, which is written in C++ and exposed to Python using the boost::python library. At the library core
is the implementation of an efficient solver for the multiple revolutions Lambert's problem, objects representing 
direct (Sims-Flanagan), indirect (Pontryagin) and hybrid methods to represent low-thrust optimization problems
, efficient keplerian propagators, Taylor-integrators, a SGP4 propagator, TLE and SATCAT support and more,  ....

pykep has been compiled and installed successfully on different platforms. pykep is also present in
the  `Python Index <https://pypi.python.org/pypi/pykep>`_ providing some precompiled modules for widely used architectures.

pykep has been used by the European Space Agency's Advanced Concepts Team during
different Global Trajectory Optimization Competitions `GTOC <http://sophia.estec.esa.int/gtoc_portal>`_
and several research papers as well as for the optimization of preliminary mission scenarion for the M-ARGO interplanetary
cubesat concept.
