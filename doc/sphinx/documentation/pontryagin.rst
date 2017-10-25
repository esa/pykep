.. _pontryagin:

=======================
The pontryagin module
=======================

This module contains the class ``PyKEP.pontryagin.leg``, which allows one to
efficiently construct low-thrust trajectories using the indirect optimal
control transcription alla moda di Pontryagin's maximum principle.
``PyKEP.pontryagin.leg`` transcribes the two-point boundary value problem
resulting from `Pontryagin's maximum principle <https://en.wikipedia.org/wiki/Pontryagin%27s_maximum_principle>`_ , in which one must correctly
choose (via and optimiser, e.g. `SNOPT <https://esa.github.io/pagmo_plugins_nonfree/py_snopt7.html>`_) the nondimensional costate
variables (nondimensional variables each corresponding to the respective
components of the state) to lead the dynamical system (i.e. spacecraft dynamics)
to chosen boundary conditions (within some tolerance). The user can choose
the parametres ``type(freemass) == bool`` and ``type(freetime) == bool``
to enforce transversality conditions on the arrival mass and time.
If ``freemass == True`` than the final mass may vary; if ``freetime == True``
than the final time may vary.

This configurability allows one to correctly construct any of the following:

- single shooting trajectory optimisation problems
- multiple shooting trajectory optimisation problem
- trajectory optimisation problems with flybys

The list of classes and the detailed documentation follows:

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`PyKEP.pontryagin.leg`                   class           represents one leg in the Pontryagin model
=========================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: PyKEP.pontryagin.leg(object)

  .. automethod:: PyKEP.pontryagin.leg.__init__(t0=None, x0=None, l0=None, tf=None, xf=None, sc=PyKEP.sims_flanagan.spacecraft(1000, 0.3, 2500), mu=PyKEP.MU_SUN, freemass=True, freetime=True, alpha=1, bound=True)

  .. automethod:: PyKEP.pontryagin.leg.set(t0, x0, l0, tf, xf)

  .. automethod:: PyKEP.pontryagin.leg.mismatch_constraints(atol=1e-5, rtol=1e-5)

  .. automethod:: PyKEP.pontryagin.leg.get_states(atol=1e-12, rtol=1e-12)

  .. automethod:: PyKEP.pontryagin.leg.plot_traj(axis, mark="k.-", atol=1e-11, rtol=1e-11, units=PyKEP.AU)

  .. automethod:: PyKEP.pontryagin.leg.plot(x, y, mark="k.-", atol=1e-12, rtol=1e-12, unitsx=1, unitsy=1, xlabel=False, ylabel=False)
