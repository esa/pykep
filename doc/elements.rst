.. _elements:

Orbital Elements
###################

In `pykep` the default osculating orbital elements used are the classical set :math:`[a, e, i, \Omega, \omega, f]`
(:math:`f` is the True Anomaly) together with the Cartesian position and velocity :math:`[\mathbf r, \mathbf v]`. Support is given also for
the set :math:`[a, e, i, \Omega, \omega, M]` (:math:`M` is the Mean Anomaly) is supported as well as
the Mean Equinoctial Elements :cite:p:`equinoctialelems` defined as:

.. math::
    \left\{
    \begin{array}{l}
    p = a (1 - e^2) \\
    f = e\cos(\omega + \Omega) \\
    g = e\sin(\omega + \Omega) \\
    h = \tan\left(\frac i2\right)\cos\Omega \\
    k = \tan\left(\frac i2\right)\sin\Omega \\
    L = \Omega + \omega + f 
    \end{array}
    \right.

These are avoid of singularities, except at :math:`i = \pi`, in which case
the retrogade version of the elements is to be used.

.. note::
    In `pykep` the convention :math:`a<0` for hyperbolas is enforced. The user will thus not be able to instantiate
    orbital elements where :math:`a(1-e) < 0`

A number of functions are provided to convert to and from the various orbital parameters. SOme also have the symbolic 
version, which can be used to compute the Jacobian of the transformation. The symbolic version of the functions
return a pair of vectors, the first one is the transformation itself, while the second one is the Jacobian,
optionally returned if the user set the flag ``jacobian`` to ``True``. 

.. currentmodule:: pykep

-----------------------------------------

.. autoclass:: el_type
   :members: 

.. autofunction:: ic2par

.. autofunction:: par2ic

.. autofunction:: ic2mee

.. autofunction:: mee2ic

.. autofunction:: mee2par

.. autofunction:: par2mee