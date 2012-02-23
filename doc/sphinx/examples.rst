.. _imagegallery:

Examples
======================

Low-thrust global optimization (Nuclear electric propulsion)
------------------------------------------------------------

.. figure:: images/gallery1.png
   :alt: "Eart-Mars low-thrust transfer"
   :align: left

This plot is produced by the following code:

.. code-block:: python

   from PyKEP import *
   kep_examples.run_example1()

after the monotonic basin hopping algorithm concludes its computations. Different images are actually produced each time as
the algorithm is non-deterministic.

The example demonstrates the use of the sims-flanagan module of PyKEP to perform global optimization of a single leg interplanetary trajectory over
large launch windows. In particular, it defines a
transfer between the Earth and Mars as an NLP global optimization problem (using the open source `PyGMO project <http://pagmo.sourceforge.net/pygmo/index.html>`_)
and it then attempts to find one solution using the Monotonic Basin Hopping meta-algorithm connected to an SQP local optimization technique (from SciPy). In case
the user has a license for SNOPT, the use of this typically result in a performance gain

The code for this example can be studied `here <http://keptoolbox.git.sourceforge.net/git/gitweb.cgi?p=keptoolbox/keptoolbox;a=blob;f=PyKEP/kep_examples/_ex1.py>`_

A reference where the methodology is studied more in depth is:

C H Yam, D D Lorenzo, D Izzo: `Low-thrust trajectory design as a constrained global optimization problem <http://pig.sagepub.com/content/early/2011/08/09/0954410011401686.abstract>`_  doi: 10.1177/0954410011401686
Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering August 10, 2011

Multi-revs Lambert Problem
--------------------------

.. figure:: images/gallery2.png
   :alt: multi-rev Lambert Problem
   :align: right

This plot is produced by the following code:

.. code-block:: python

   from PyKEP import *
   kep_examples.run_example2()

The example demonstrates the basic use of the multiple revolution Lambert solver. In particualar, it defines a fixed geometry of the
Earth-Mars constellation and plots all possible resulting solution to the Lambert's problem. The chosen geometry is such that one revolution solutions
exists and are thus visualized.

While not published yet, the multi-rev solver implemented in PyKEP is extremely efficient, if anyone knows about a faster implementation of
a Lambert's solver, please contact us. In our solver we use a transformation that deforms the time of flight curves into lines, thus allowing uniform
consistent convergence properties for all possible geometries.

The code for this example can be studied `here <http://keptoolbox.git.sourceforge.net/git/gitweb.cgi?p=keptoolbox/keptoolbox;a=blob;f=PyKEP/kep_examples/_ex2.py>`_

LT-MGA global optimization (Nuclear electric propulsion)
---------------------------------------------------------

.. figure:: images/gallery3.png
   :alt: "Eart-Venus-Mercury low-thrust transfer"
   :align: left

This plot is produced by the following code:

.. code-block:: python

   from PyKEP import *
   kep_examples.run_example3()

after the monotonic basin hopping algorithm concludes its computations. Different images are actually produced each time as
the algorithm is non-deterministic.

The example demonstrates the use of the sims-flanagan module of PyKEP to perform global optimization of a multiple leg interplanetary trajectory over
large launch windows. In particular, it defines a
transfer between Earth and Mercury making use of a Venus intermediate fy-by as an NLP global optimization problem (using the open source `PyGMO project <http://pagmo.sourceforge.net/pygmo/index.html>`_)
and it then attempts to find one solution using the Monotonic Basin Hopping meta-algorithm connected to an SQP local optimization technique (from SciPy). In case
the user has a license for SNOPT, the use of this typically result in a performance gain

The code for this example can be studied `here <http://keptoolbox.git.sourceforge.net/git/gitweb.cgi?p=keptoolbox/keptoolbox;a=blob;f=PyKEP/kep_examples/_ex3.py>`_

