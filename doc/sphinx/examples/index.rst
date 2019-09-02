Examples
======================

In this page we provide several examples that should just run out of the box if you have all the dependencies
installed correctly. On top of pykep you may need to install pygmo, matplotlib, scipy, pygmo_plugins_nonfree. 
Many example run using the commercial optimizer snopt7 (wrapped by pygmo_plugins_nonfree), 
in case you hacve no access to it the NLOPT optimizer slsqp as well as ipopt may be used. Both of which
are shipped in the pygmo optimization suite.

Some of the examples are the outcome of research we performed for the European Space Agency and contain rather
advanced optimization techniques and methods. Please quote our papers in case you find them of use in your research.

Generic
^^^^^^^
.. toctree::
  :maxdepth: 2

  ex1.ipynb
  ex2.ipynb
  ex12.ipynb
  gravity_spherical_harmonic.ipynb

Advanced Research Topics (requires pygmo)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
  :maxdepth: 2

  ex3.ipynb
  ex4
  ex5
  ex6
  ex7
  ex8
  ex9
  ex10
  ex11