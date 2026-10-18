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

import unittest as _ut


class plot_tests(_ut.TestCase):
    def test_add_solar_system_outer_planets(self):
        """We test that Uranus and Neptune can be added and use their own marker sizes."""
        import matplotlib.pyplot as plt

        # One distinct marker size per body, so we can tell which entry was used
        sizes = [15, 2, 3, 3, 2, 8, 8, 5, 6]
        ax = _pk.plot.make_3Daxis()
        _pk.plot.add_solar_system(ax, bodies=[7, 8], when=_pk.epoch(0), s=sizes)

        # The planet markers are the scatter collections carrying the planet name as label
        scatters = {c.get_label(): c for c in ax.collections}
        self.assertEqual(scatters["Uranus"].get_sizes()[0], sizes[7])
        self.assertEqual(scatters["Neptune"].get_sizes()[0], sizes[8])
        plt.close(ax.figure)
