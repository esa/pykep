# Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
#                          Advanced Concepts Team, European Space Agency (ESA)
#
# This file is part of the pykep library.
#
# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pykep as _pk

class cr3bp:
    """User-defined planet-like interface for a CR3BP reference trajectory.

    This UDPLA exposes a simple planet-like interface backed by a CR3BP
    integrator. It returns ephemerides in SI units consistent with the
    length/time scales supplied at construction, and implements the
    optional ``acc`` method required by gradient-enabled problems.

    Args:
        *when* (:class:`float` or :class:`pykep.epoch`): Reference epoch (MJD2000) for the provided
            reference state. If a :class:`pykep.epoch` is provided, its ``mjd2000`` value is used.

        *state_nd* (:class:`list` or :class:`numpy.ndarray`): Reference 6D state in non-dimensional
            CR3BP units used to initialize the internal integrator: ``[x,y,z,vx,vy,vz]``.

        *mu_cr3bp* (:class:`float`): The CR3BP mass parameter (non-dimensional).

        *TIME* (:class:`float`): Time scale (non-dimensional seconds) used to convert between
            MJD2000/seconds and the integrator's internal time units.

        *L* (:class:`float`): Length scale used to map integrator positions to SI units.

        *name* (:class:`str`, optional): Human-readable name for the UDPLA. If omitted a default
            descriptive name containing ``mu_cr3bp`` is assigned.

        *tol* (:class:`float`, optional): Tolerance passed to the internal Taylor adaptive integrator.

    """
    def __init__(self, when, state_nd, mu_cr3bp, TIME, L, name = None, tol = 1e-16):
        if name is None:
            name = "An unkown body in the CR3BP with mu = " + str(mu_cr3bp)
        import heyoka as hy
        # when can be a mjd2000 or a pykep epoch, consistently with other pykep planet interfaces.
        if isinstance(when, _pk.epoch):
            self.ref_mjd2000 = when.mjd2000
        else:
            self.ref_mjd2000 = when
        self.ref_state = state_nd
        self.mu_cr3bp = mu_cr3bp
        # These are the nd units declared by the user. 
        self.TIME = TIME
        self.L = L
        self.V = self.L / self.TIME
        self.ACC = self.V / self.TIME
        self.name = name
        # We use the nd numerical integrator of the cr3bp given by pykep.
        self.ta = _pk.ta.get_cr3bp(tol)
        # Set the parameters and the initial state of the integrator.
        self.ta.pars[0] = self.mu_cr3bp
        self.reset_ta()
        # The cfunc for the acceleration to implement the acc method of the pykep planet interface.
        dyn = _pk.ta.cr3bp_dyn()
        x,y,z,vx,vy,vz = hy.make_vars("x", "y", "z", "vx", "vy", "vz")
        self.acc_cfunc = hy.cfunc([e[1] for e in dyn[3:]], vars = [x,y,z,vx,vy,vz])

    def reset_ta(self):
        # This method allows to reset the internal integrator to the reference state and epoch,
        # which can be useful for debugging or to limit the possible loss of precision if the user calls 
        # repetedly eph with non-monotone epochs.
        self.ta.state[:] = self.ref_state
        self.ta.time = self.ref_mjd2000 * _pk.DAY2SEC / self.TIME
        
    def eph(self, mjd2000):
        """Return position and velocity at the requested epoch.

        Args:
            *mjd2000* (:class:`float`): Epoch in MJD2000 days.

        Returns:
            :class:`list` [:class:`list`, :class:`list`]: Position and velocity in SI units (meters, meters/second).
        """
        # We convert to nd units
        epoch_nd = mjd2000 * _pk.DAY2SEC / self.TIME
        # Propagate until the desired epoch. We use propagate until so that if eph is called multiple times we avoid propagating more than needed.
        # This may have ramifications on the accuracy of the results if the user calls eph with non-monotone epochs
        # but we assume that this is not an issue.
        self.ta.propagate_until(epoch_nd)
        return [self.ta.state[:3] * self.L, self.ta.state[3:] * self.V]
    
    def acc(self, mjd2000):
        """Return acceleration from the CR3BP dynamics at the requested epoch.

        Args:
            *mjd2000* (:class:`float`): Epoch in MJD2000 days.

        Returns:
            :class:`list` or :class:`numpy.ndarray`: Acceleration at the epoch, returned in SI units (meters/second^2).
        """
        epoch_nd = mjd2000 * _pk.DAY2SEC / self.TIME
        self.ta.propagate_until(epoch_nd)  # avoid double propagation when possible
        return self.acc_cfunc(self.ta.state, pars=self.ta.pars) * self.ACC
    
    def get_name(self):
        """Return the human-readable name for this UDPLA."""
        return self.name
    
    def get_extra_info(self):
        """Return a short string with extra information about the UDPLA.

        Typically includes the reference state and epoch used to initialize the integrator.
        """
        retval = "Reference state (nd): " + str(self.ref_state)
        retval += "\nReference epoch (mjd2000): " + str(self.ref_mjd2000)
        return retval
        