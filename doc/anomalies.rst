.. _anomalies:

Anomalies Conversions
########################

In `pykep` we adopt the following naming for the various anomalies: 
:math:`M` is the Mean Anomaly,
:math:`E` is the Eccentric Anomaly,
:math:`L` is the True Longitude,
:math:`\lambda` is the Mean Longitude,
:math:`H` is the Hyperbolic Anomaly,
:math:`N` is the Mean Hyperbolic Anomaly,
:math:`\zeta` is the Gudermannian,
in functions or variable names these symbols can be spelled out (like zeta for :math:`\zeta`)
or made lowercase like m for :math:`M`.

Below we list various conversion functions that allow to convert from one anomaly
to another, and their vectorized versions.


.. currentmodule:: pykep

Normal
###############################################################################

.. autofunction:: m2e

.. autofunction:: e2m

.. autofunction:: m2f

.. autofunction:: f2m

.. autofunction:: e2f

.. autofunction:: f2e

.. autofunction:: n2h

.. autofunction:: h2n

.. autofunction:: n2f

.. autofunction:: f2n

.. autofunction:: h2f

.. autofunction:: f2h

.. autofunction:: zeta2f

.. autofunction:: f2zeta

Vectorized
###############################################################################
 
.. autofunction:: m2e_v
    
.. autofunction:: e2m_v

.. autofunction:: m2f_v

.. autofunction:: f2m_v

.. autofunction:: e2f_v

.. autofunction:: f2e_v

.. autofunction:: n2h_v

.. autofunction:: h2n_v

.. autofunction:: n2f_v

.. autofunction:: f2n_v

.. autofunction:: h2f_v

.. autofunction:: f2h_v


.. autofunction:: zeta2f_v

.. autofunction:: f2zeta_v


