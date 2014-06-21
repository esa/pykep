MGA 1DSM global optimization
==========================================================

.. figure:: ../images/gallery5.png
   :align: left
   
.. figure:: ../images/gallery5b.png
   :alt: "Earth-Venus-Earth transfer via impulsive maneuvers"
   :align: right

These plots are produced by the following code:

.. code-block:: python

   from PyKEP import *
   examples.run_example5()

after the self-adaptive differential evolution algorithm concludes its computations. 

The example demonstrates the use of the mga_1dsm problem constructor of the *interplanetary* module. 
This helper class allows to construct a PyGMO global optimization problem (`PyGMO project <http://pagmo.sourceforge.net/pygmo/index.html>`_) 
representing a Multiple Gravity Assist Interplanetary Trajectory where only one Deep Space Manouvre is allowed in each leg. The problem is quickly constructed as either a single objective (total DV) or
even a multiobjective optimization problem (DV and time of flight)

In the particular instance constructed in the example, we have an Earth-Venus-Earth transfer with a Venus intermediate
fly-by. he selected launch windows is 2016-2017, the launch hyperbolic velocity is 2.5 km/s and the allowed time of flight is in [0.5,3] years.

To solve such a global optimization problem, we use jDE, a self-adaptive version of differential evolution and the generalized migration
operator, thus computing the solution in eight paralel threads.

The code for this example can be studied `here. Feel free to leave comments.
<https://github.com/esa/pykep/blob/master/PyKEP/examples/_ex5.py>`_
