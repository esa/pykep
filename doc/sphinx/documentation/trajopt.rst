.. _trajopt:

===================================
The trajopt module (requires PyGMO)
===================================

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`PyKEP.trajopt.mga_1dsm`                 class           A multiple Gravity Assist Trajectory with one deep space manouvre 
:class:`PyKEP.trajopt.mga_lt_nep`               class           A multiple Gravity Assist Trajectory low-thrust optimization problem
:class:`PyKEP.trajopt.mr_lt_nep`                class           A multiple ransezvous low-thrust optimization problem
=========================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: PyKEP.trajopt.mga_1dsm(*args)
 
   .. automethod:: PyKEP.trajopt.mga_1dsm.__init__(*args)
   
   .. automethod:: PyKEP.trajopt.mga_1dsm.set_tof(*args)
   
   .. automethod:: PyKEP.trajopt.mga_1dsm.set_launch_window(*args)
   
   .. automethod:: PyKEP.trajopt.mga_1dsm.set_vinf(*args)
   
   .. automethod:: PyKEP.trajopt.mga_1dsm.pretty(*args)
   
   .. automethod:: PyKEP.trajopt.mga_1dsm.plot(*args)
   
.. autoclass:: PyKEP.trajopt.mga_lt_nep(*args)
 
   .. automethod:: PyKEP.trajopt.mga_lt_nep.__init__(*args)
   
   .. automethod:: PyKEP.trajopt.mga_lt_nep.high_fidelity(*args)
   
   .. automethod:: PyKEP.trajopt.mga_lt_nep.ic_from_mga_1dsm(*args)
       
   .. automethod:: PyKEP.trajopt.mga_lt_nep.plot(*args)

.. autoclass:: PyKEP.trajopt.mr_lt_nep(*args)
 
   .. automethod:: PyKEP.trajopt.mr_lt_nep.__init__(*args)
       
   .. automethod:: PyKEP.trajopt.mr_lt_nep.plot(*args)