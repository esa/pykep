import pykep as _pk


def planet_to_keplerian(pla, when: _pk.epoch, mu=None):
    """Transforms a planet to its Keplerian version

    The :class:`~pykep.planet` returned will be Keplerian, with elements identical to
    the osculating elements at *t0* of the input planet. Hence its ephemerides can be computed
    at any epoch, but will only be matching the original ones at *t0*.

    This utility is useful to extend the validity range of some :class:`~pykep.planet`'s ephemerides.

    Args:
        *pla* (:class:`~pykep.planet`): the input planet.
        
        *when* (:class:`~pykep.epoch`): the epoch to match the osculating elements.
        
        *mu* (:class:`float`, optional): the central body parameter. Defaults to the one computed from the input :class:`~pykep.planet`.

    Returns:
        :class:`~pykep.planet`: a Keplerian planet.
    """
    if mu is None:
        mu = pla.get_mu_central_body()
        if mu == -1:
            print(
                f"PyKEP ERROR: Cannot get the central body parameter from {pla.get_name()}, please define explicitly a mu."
            )
            raise ValueError
    posvel = pla.eph(when)
    udpla = _pk.udpla.keplerian(
        when=when,
        posvel=posvel,
        mu_central_body=mu,
        name=pla.get_name() + "(K)",
        added_params=[pla.get_mu_self(), pla.get_radius(), pla.get_safe_radius()],
    )
    return _pk.planet(udpla)
