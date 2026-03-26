## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)
## This file is part of the kep3 library.
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

from .core import epoch
from datetime import datetime, timezone, timedelta
import pykep as pk


def _to_datetime(self):
    """
    Convert this ``epoch`` to a timezone-aware Python ``datetime`` in UTC.

    This method is added to the class via monkey-patching:

    - It interprets ``self.mjd2000`` as a day count relative to the MJD2000
      reference epoch (2000-01-01 00:00:00) in UTC.
    - It converts days to seconds by multiplying by 86400.
    - It returns an *aware* ``datetime`` (``tzinfo=timezone.utc``) computed as::

          datetime(2000, 1, 1, tzinfo=timezone.utc) + timedelta(seconds=self.mjd2000 * 86400)

    Notes:
        This is a civil-time conversion that assumes the MJD2000 epoch and the
        stored day count are meant to be interpreted in UTC. If your underlying
        time scale is TT/TAI (or you require leap-second-aware conversions), use
        an astronomy time library (e.g., Astropy) instead of naive arithmetic.

    Returns:
        datetime: A timezone-aware ``datetime`` in UTC corresponding to this epoch.
    """
    return datetime(2000, 1, 1, tzinfo=timezone.utc) + timedelta(
        seconds=self.mjd2000 * 86400.0
    )


# Do the actual patching.
setattr(epoch, "to_datetime", _to_datetime)
