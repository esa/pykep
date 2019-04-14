pykep
=====

[![Join the chat at https://gitter.im/esa/pykep](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/esa/pykep?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/esa/pagmo.svg?branch=master)](https://travis-ci.org/esa/pykep) [![Code Health](https://landscape.io/github/esa/pykep/master/landscape.svg?style=flat)](https://landscape.io/github/esa/pykep/master)

pykep is a scientific library providing basic tools for astrodynamics research. Algorithmic efficiency is
a main focus of the library, which is written in C++ and exposed to Python using the boost::python library. At the library core
is the implementation of an efficient solver for the multiple revolutions Lambert's problem, objects representing 
direct (Sims-Flanagan), indirect (Pontryagin) and hybrid methods to represent low-thrust optimization problems
, efficient keplerian propagators, Taylor-integrators, a SGP4 propagator, TLE and SATCAT support and more.

Check the official documentation at https://esa.github.io/pykep/
