.. _propagation:

Numerical Propagation
========================

The backbone of numerical propagation in `pykep` is based on Lagrangian coefficients for 
Kepler's  dynamics and Taylor numerical integration, as implemented in the 
`Heyoka <https://bluescarni.github.io/heyoka.py/index.html>`_ :cite:p:`biscaniheyoka1` python package, for all
other cases. The state transition matrix is also available and provided, in the case of numerical integration,
seamlessly via variational equations.

The main routines are listed here:
 
Keplerian 
---------------
 
.. currentmodule:: pykep

.. autofunction:: propagate_lagrangian

.. autofunction:: propagate_lagrangian_grid

---------------------------------------------------------

Taylor adaptive propagators
-----------------------------

In addition to the above mentioned propagators,
Taylor adaptive integrators are offered in `pykep` wrapping some of the functionalities of
`Heyoka <https://bluescarni.github.io/heyoka.py/index.html>`_ :cite:p:`biscaniheyoka1` python package.
Their variational version is also offered (at order one) as to be able to produce stms and, where needed,
other useful quantities. Higher order variational equations can also be obtained directly using the available dynamics and 
using `Heyoka <https://bluescarni.github.io/heyoka.py/index.html>`_ :cite:p:`biscaniheyoka1` syntax.

Some of the Taylor integrators are "zero hold" versions of a given dynamics, meaning that they
consider a constant thrust vector in some frame and include the mass variation due to said thrust.
These are intended to be used to build "zero hold" trajectory legs and direct methods to optimally
control the spacecraft.

Some of the Taylor adaptive integrators are associated to OCPs (Optimal Control Problems) of
relevance to interplanetary flight. In particular they are born when applying Pontryagin 
principle to dynamics of interests and result in two point boundary value problems (TPBVP).
In these cases, a number of auxiliary functiond describing
the control, the Hamiltonian, the switching function, etc., are also provided.

--------------------------------------------------------

Two Body Problem (Kepler)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: pykep.ta

.. autofunction:: get_kep

.. autofunction:: get_kep_var

.. autofunction:: kep_dyn 


Circular Restricted Three Body Problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: pykep.ta

.. autofunction:: get_cr3bp

.. autofunction:: get_cr3bp_var

.. autofunction:: cr3bp_dyn 

.. autofunction:: cr3bp_jacobi_C 

.. autofunction:: cr3bp_effective_potential_U

Bicircular Problem
~~~~~~~~~~~~~~~~~~~~~

Introduced by Simo' et al. in his '97 paper :cite:p:`simo1995bicircular`.

.. autofunction:: get_bcp 
    
.. autofunction:: get_bcp_var

.. autofunction:: bcp_dyn 

Zero Hold Keplerian Propagator in Cartesian coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_zoh_kep

.. autofunction:: get_zoh_kep_var

.. autofunction:: zoh_kep_dyn 

Zero Hold Keplerian Propagator in Equinoctial elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_zoh_eq

.. autofunction:: get_zoh_eq_var

.. autofunction:: zoh_eq_dyn 

Zero Hold CR3BP Propagator in Cartesian coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_zoh_cr3bp

.. autofunction:: get_zoh_cr3bp_var

.. autofunction:: zoh_cr3bp_dyn 

Zero Hold Solar Sail Propagator in Cartesian coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_zoh_ss

.. autofunction:: get_zoh_ss_var

.. autofunction:: zoh_ss_dyn 

Low-thrust Pontryagin Cartesian TPBVP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_pc

.. autofunction:: get_pc_var

.. autofunction:: pc_dyn

Low-thrust Pontryagin Equinoctial TPBVP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_peq

.. autofunction:: get_peq_var

.. autofunction:: peq_dyn 
