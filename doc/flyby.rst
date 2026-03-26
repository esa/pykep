.. _flyby:

Fly-by routines
########################
The gravity assist technique, widely used in the Pioneer and Voyager missions, was developed by Michael Andrew Minovitch when he was a UCLA 
graduate student working as intern at NASA's Jet Propulsion Laboratory. Many authors later refined the technique which, avoiding to solve
directly the full three body problem, is based on patching the planetocentric hyperbola with incoming and outgoing (elliptic) trajectories.
In `pykep` we offer the basic routines that allow to use this technique in the context of trajectory optimization and patched conics propagations.

.. currentmodule:: pykep

------------------------------------------------------

Minovitch fly-by
****************
.. autofunction:: fb_con

.. autofunction:: fb_dv

.. autofunction:: fb_vout 