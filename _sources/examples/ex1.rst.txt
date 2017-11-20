Multiple impulses transfer between Earth and Venus
============================================================

.. figure:: ../images/gallery1.png
   :alt: "Earth-Venus 3 impulse transfer"
   :align: left

This plot is produced by the following code:

.. code-block:: python

   import pykep as pk
   pk.examples.run_example1(impulses = 3)

after a call to the Covariance Matrix Adaptation Evolutionary Strategy (CMA-ES) concludes its computations. 
Different images are produced at each run as the algorithm is non-deterministic and problem has several local optima.

This example demonstrates the use of pykep building blocks (used in the class pl2pl_N_impulses available from the trajopt module) to assemble an unconstrained optimization problem representing
an impulsive transfer between the Earth and Venus. Multiple impulses are allowed (the optimal number is in this case 3) allowing
to study the complexification of the objective fuinction landscape as the number of impulses increase.

The problem is solved using the pygmo optimization framework, and in particular the algorithm CMA-ES.

The code for this example can be studied `here. 
<https://github.com/esa/pykep/blob/master/pykep/examples/_ex1.py>`_ Feel free to leave comments.
