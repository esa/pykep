1DSM-MGA global optimization (Nuclear electric propulsion)
==========================================================

.. figure:: ../images/gallery5.png
   :alt: "Earth-Venus-Earth transfer via impulsive maneuvers"
   :align: left

This plot is produced by the following code:

.. code-block:: python

   from PyKEP import *
   kep_examples.run_example5()

after the self-adaptive differential evolution algorithm concludes its computations. Different images are actually
produced each time as the algorithm is non-deterministic.

The example demonstrates the use of deep space maneuvers to perform global optimization for multiple
gravity assist trajectories. In particular, it defines an Earth-Venus-Earth transfer with a Venus intermediate
fly-by as a global optimization problem using a differential evolution algorithm (using the open source
`PyGMO project <http://pagmo.sourceforge.net/pygmo/index.html>`_)
The launch window is selected in the year 2016 for this specific example. 

The code for this example can be studied `here <http://keptoolbox.git.sourceforge.net/git/gitweb.cgi?p=keptoolbox/keptoolbox;a=blob;f=PyKEP/kep_examples/_ex5.py>`_