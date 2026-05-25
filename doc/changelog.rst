.. _changelog:

Changelog
#########

.. _changelog_3_1_0:

3.1.0 (unreleased)
==================

- Added CMake option ``kep3_BUILD_CPP_LIBRARY`` (default ``ON``). When set to
  ``OFF``, the C++ library is not built from source and CMake will instead
  locate an existing installation via ``find_package``. A descriptive
  ``FATAL_ERROR`` is emitted if the library cannot be found. This makes it
  possible to build only the Python bindings, tests, or benchmarks against an
  already-installed ``kep3``.

.. _changelog_3_0_0:

3.0.0 (23-05-2026)
==================

pykep 3 is a complete reimagination of the library: a new C++23 core, a clean
pybind11 Python interface, type-erased *user-defined* abstractions (``udpla``,
``udp``, ``uda``), native integration with pagmo and heyoka, and an expanded
set of low-thrust transcriptions and benchmark problems.
