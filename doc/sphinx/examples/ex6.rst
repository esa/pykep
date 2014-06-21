Asteroid hopping in the main belt (requires PyGMO)
==========================================================

.. figure:: ../images/gallery6.jpg
   :align: left
   

These plots are produced by the following code:

.. code-block:: python

   from PyKEP import *
   examples.run_example6()

if a feasible solution is reached. Repeat several times in case unfeasible solutions are returned.

This example demonstrates the use of the mr_lt_nep class in the trajopt module. The class derives from
PyGMO problem.base and represents the optimization of a multiple randezvous mission of a low-thrust spacecraft.
The initial guess is taken from one of the trajectory submitted by the ACT/ISAS team to the 7th edition
of the global trajectory optimisation competition <http://sophia.estec.esa.int/gtoc_portal/>`_.


The code for this example can be studied `here. Feel free to leave comments.
<https://github.com/esa/pykep/blob/master/PyKEP/examples/_ex6.py>`_
