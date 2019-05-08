# -*- coding: utf-8 -*-

# Copyright 2017-2018 PaGMO development team
#
# This file is part of the PaGMO library.
#
# The PaGMO library is free software; you can redistribute it and/or modify
# it under the terms of either:
#
#   * the GNU Lesser General Public License as published by the Free
#     Software Foundation; either version 3 of the License, or (at your
#     option) any later version.
#
# or
#
#   * the GNU General Public License as published by the Free Software
#     Foundation; either version 3 of the License, or (at your option) any
#     later version.
#
# or both in parallel, as here.
#
# The PaGMO library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received copies of the GNU General Public License and the
# GNU Lesser General Public License along with the PaGMO library.  If not,
# see https://www.gnu.org/licenses/.

from __future__ import absolute_import as _ai

import unittest as _ut

class core_functions_test_case(_ut.TestCase):
    """Test case for the core functions

    """

    def runTest(self):
        self.run_epoch_test()

    def run_epoch_test(self):
        from .core import epoch
        ep = epoch(julian_date=0, julian_date_type='mjd2000')

class lambert_test_case(_ut.TestCase):
    """Test case for the lambert problem class

    """

    def runTest(self):
        from .core import lambert_problem
        lp = lambert_problem(r1 = [1,1,0], r2 = [0,1,0], mu = 1., cw = False, max_revs = 0, tof = 0.3)


def run_test_suite(level=0):
    """Run the full test suite.

    This function will raise an exception if at least one test fails.

    Args:
        level(``int``): the test level (higher values run longer tests)

    """
    
    retval = 0
    suite = _ut.TestLoader().loadTestsFromTestCase(core_functions_test_case)
    suite.addTest(lambert_test_case())
    
    test_result = _ut.TextTestRunner(verbosity=2).run(suite)

    
    if len(test_result.failures) > 0 or len(test_result.errors) > 0:
        retval = 1
    if retval != 0:
        raise RuntimeError('One or more tests failed.')
