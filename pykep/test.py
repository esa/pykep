"""Test runner entrypoint for the pykep package.

This module intentionally contains only the `run_test_suite` function so that
`python -c "import pykep.test as t; t.run_test_suite()"` works both when used
from the repository and when the package is installed.

All actual TestCase classes must live in separate files named `test_*.py`
inside the `pykep` package (for example `test_core.py`,
`test_leg_sims_flanagan.py`, etc.).
"""

import os
import warnings
import unittest as _ut


def run_test_suite(verbosity=2, suppress_solver_warnings=True):
    """Discover and run tests located in `pykep/test_*.py`.

    Parameters
    - verbosity: passed to the test runner
    - suppress_solver_warnings: if True, filter generic Warning subclasses
    """
    import inspect

    # Discover tests inside the installed `pykep` package directory. This
    # works both for the installed package (site-packages) and for a
    # development copy when PYTHONPATH or editable install is used.
    import pykep as _pk

    tests_dir = _pk.__path__[0]
    package_root = os.path.dirname(tests_dir)

    suite = _ut.defaultTestLoader.discover(
        start_dir=tests_dir, pattern="test_*.py", top_level_dir=package_root
    )

    class _FilePrefixedTextResult(_ut.TextTestResult):
        def getDescription(self, test):
            import os

            cls = test.__class__
            cls_name = cls.__name__
            if cls_name.startswith("test_"):
                cls_name = cls_name[len("test_") :]

            try:
                file_path = inspect.getfile(cls)
                filename = os.path.basename(file_path)
            except Exception:
                filename = getattr(cls, "__module__", "")

            return f"{filename}: {cls_name}.{test._testMethodName}"

    if suppress_solver_warnings:
        warnings.filterwarnings("ignore", category=Warning)

    runner = _ut.TextTestRunner(verbosity=verbosity, resultclass=_FilePrefixedTextResult)
    return runner.run(suite)
