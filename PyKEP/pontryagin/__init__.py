"""
This module contains the class `PyKEP.pontryagin.leg`, which allows one to
efficiently construct low-thrust trajectories using the indirect optimal
control transcription alla moda di Pontryagin's maximum principle.
`PyKEP.pontryagin.leg` transcibes the two-point boundary value
problem resulting from Pontryagin's maximum principle, in which one must
correctly choose (via and optimiser, e.g. SNOPT) the nondimensional costate
variables (nondimensional variables each corresponding to the respective
components of the state) to lead the dynamical system (i.e. spacecraft dynamics)
to chosen boundary conditions (within some tolerance).
The user can choose the parametres `type(freemass) == bool` and `type(freetime) == bool`
to enforce transversality conditions on the arrival mass and time.
If `freemass == True` than the final mass may vary; if `freeimte == True`
than the final time may vary.

This configurability allows one to correctly construst any of the following:
- single shooting trajectory optimisation problems
- multiple shooting trajectory optimisation problem
- trajectory optimisation problems with flybys
"""

from PyKEP.pontryagin._leg import leg
