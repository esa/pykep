.. _changelog:

Changelog
#########

.. _changelog_3_1_0:

3.1.0 (unreleased)
==================

Build system
------------

- Added CMake option ``kep3_BUILD_CPP_LIBRARY`` (default ``ON``). When set to
  ``OFF``, the C++ library is not built from source and CMake will instead
  locate an existing installation via ``find_package``. A descriptive
  ``FATAL_ERROR`` is emitted if the library cannot be found. This makes it
  possible to build only the Python bindings, tests, or benchmarks against an
  already-installed ``kep3``.

Bug fixes
---------

- Fixed pickling of TOPS gym classes. The lambdas passed as ``state2cart``
  and equivalent callbacks to the internal UDP prevented serialization. They were
  replaced with :func:`functools.partial` wrapping small module-level helper
  functions (``_cart_scale``, ``_mee_to_cart``, ``_cfunc_call``), which are
  picklable by name.

- Fixed pickling of :class:`~pykep.leg.zoh`. The C++ class is bound as
  ``_zoh_cpp`` in ``pykep.core``, so pickle looked for ``pykep.leg._zoh_cpp``
  when reconstructing objects. The alias was missing from ``pykep/leg/__init__.py``
  (only ``_zoh`` was defined). Added ``_zoh_cpp = _core._zoh_cpp`` to make the
  lookup succeed, consistent with the existing ``_sims_flanagan`` and
  ``_sims_flanagan_alpha`` aliases.

- Fixed several bugs in :class:`~pykep.trajopt.pl2pl_N_impulses`: a crash when
  ``phase_free=False`` and ``t0_bounds=None`` (``t0`` was assigned instead of
  ``t0_bounds``); a dead/unreachable ``N_max < 2`` guard; ``plot_primer_vector``
  using ``MU_SUN`` instead of ``self._common_mu``; alpha lower bounds of ``0.0``
  in phase-free mode (now ``1e-3``, consistent with the non-phase-free path, to
  avoid ``log(0)`` in ``alpha2direct``). The ``pretty`` method was also corrected
  to print times-of-flight in days (the internal ``decode`` stores durations in
  seconds, and the conversion was missing).

.. _changelog_3_0_0:

3.0.0 (23-05-2026)
==================

pykep 3 is a complete reimagination of the library: a new C++23 core, a clean
pybind11 Python interface, type-erased *user-defined* abstractions (``udpla``,
``udp``, ``uda``), native integration with pagmo and heyoka, and an expanded
set of low-thrust transcriptions and benchmark problems.
