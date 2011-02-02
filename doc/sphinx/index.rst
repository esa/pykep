.. PyKEP documentation master file, created by
   sphinx-quickstart on Thu Nov  4 12:34:23 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================================
Welcome to PyKEP's documentation!
==================================
+----------------------------+--------------------------------+
|                            |                                |
+----------------------------+                                |
| .. image:: images/logo.png | .. image:: images/logo_act.png |
+----------------------------+                                |
|                            |                                |
+----------------------------+--------------------------------+

PyKEP is a scientific library providing basic tools for astrodynamics research. Algoritmic efficiency is 
a main focus of the library, which is written in C++ and exposed to Python using the boost::python library. At the library core
is the implementation of an efficient solver for the multiple revolutions Lambert's problem, objects representing the Sims-Flanagan low-thrust model,
efficient keplerian propagators, Taylor-integrators and more .... in its current version the library is only supporting unix systems (tested on various linux
distributions and MacOSX).

PyKEP has been used by the European SPace Agency's Advanced Concepts Team during the Global Trajectory Optimization
Competitions (GTOC), in particular `GTOC4 <http://cct.cnes.fr/cct02/gtoc4/index.htm>`_ and `GTOC5 <http://gtoc5.math.msu.su/>`_.

Contents of these documentation pages:

.. toctree::
  :maxdepth: 2

  system
  downloading
  installation
  documentation
  credits
