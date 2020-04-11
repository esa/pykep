#include <string>

#include "planet_docstrings.hpp"

namespace pykep
{

std::string planet_spice_doc()
{
    return R"(

pykep.planet.spice(target, observer, ref_frame, aberrations, mu_central_body, mu_self, radius, self_radius)

- target: Target body (see `NAIF IDs <http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html#NAIF%20Object%20ID%20numbers>`_)
- observer:		Observer body (see `NAIF IDs <http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html#NAIF%20Object%20ID%20numbers>`_
- ref_frame:		The reference frame (see `SPICE supported frames <http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html#Frames%20Supported%20in%20SPICE>`_)
- aberrations: 	Aberration correction type (see spkezr_c docs <http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html>`_)
- mu_central_body: Gravity parameter of the central body (this is not necessary to compute the ephemerides).
- mu_self: Gravity parameter of the target (this is not necessary to compute the ephemerides).
- radius: Radius of target body (this is not necessary to compute the ephemerides).
- self_radius: Safe radius of target body (this is not necessary to compute the ephemerides).

.. note::
    The presence in memory of the needed SPICE kernel files are only checked upon call to the ephemerides method. 
    As a consequence this object can still be constructed with invalid arguments. Only later the ephemerides call
    will fail throwing an excpetion. To load SPICE kernels use the :py:func:`pykep.utils.load_spice_kernel` utility.

.. note::
    mu_central_body must be set to compute the period of the orbital elements.

Example::


    import pykep as pk
    pk.load_spice_kernel("de405.bsp")
    planet = pk.planet.spice('EARTH', 'SUN', 'ECLIPJ2000', 'NONE', MU_SUN, MU_EARTH, EARTH_R, EARTH_R * 1.05)
)";
}





} // namespace pykep
