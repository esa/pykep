.. _installation:

Installation
############

.. _installation_deps:

Dependencies
============

pykep has both C++ and Python dependencies.

On the **C++ side**, pykep depends on:

* the `heyoka C++ library <https://github.com/bluescarni/heyoka>`_ (version ≥ 7, mandatory),
* the `Boost <https://www.boost.org/>`_ C++ libraries (version ≥ 1.73, mandatory),
* the `{fmt} <https://fmt.dev/>`_ library (version ≥ 10, mandatory),
* the `spdlog <https://github.com/gabime/spdlog>`_ library (version ≥ 1.12, mandatory),
* the `xtensor <https://xtensor.readthedocs.io/>`_ and ``xtensor-blas`` libraries (version ≥ 0.26, mandatory),
* the `NLopt <https://nlopt.readthedocs.io/>`_ nonlinear optimization library (mandatory).

On the **Python side**, pykep requires Python ≥ 3.11 and depends on:

* `NumPy <https://numpy.org/>`_ (mandatory),
* `pygmo <https://esa.github.io/pygmo2/>`_ (mandatory, for trajectory optimization problem interfaces),
* `heyoka.py <https://bluescarni.github.io/heyoka.py/>`_ (version ≥ 7, mandatory, for Taylor integration),
* `spiceypy <https://spiceypy.readthedocs.io/>`_ (optional, for JPL SPICE ephemeris support),
* `sgp4 <https://pypi.org/project/sgp4/>`_ (optional, for SGP4 propagation of Earth-orbiting objects),
* `matplotlib <https://matplotlib.org/>`_ (optional, for the built-in plotting utilities),
* `scipy <https://scipy.org/>`_ (optional, for certain utility functions).

.. _installation_packages:

Packages
========

.. _installation_conda:

conda
-----

conda is the **recommended** installation method, as it provides a fully managed
scientific Python stack with all C++ dependencies resolved automatically.

Add ``conda-forge`` to your channels (if not already present) and install:

.. code-block:: console

   $ conda config --add channels conda-forge
   $ conda config --set channel_priority strict
   $ conda install pykep

The conda packages are maintained by the core development team and updated with
each new release.

.. _installation_pip:

pip
---

Binary wheels for pykep are also available on `PyPI <https://pypi.org/project/pykep/>`_:

.. code-block:: console

   $ pip install pykep

Prebuilt wheels are currently provided for Linux (x86-64 and arm64).

.. note::

   pykep bundles several C++ libraries inside its pip wheels. There is a
   non-negligible chance of conflicts with other packages that bundle the same
   C++ libraries. Installing via conda is strongly preferred whenever possible.

.. _installation_source:

Installation from source
========================

Building from source is the recommended path when working on the development
version of pykep or when you need a custom build configuration.

.. warning::

   When building from source, make sure the compiler used is ABI-compatible
   with the compilers used to build pykep's dependencies. Mixing incompatible
   compilers can cause hard-to-diagnose build or runtime errors.

pykep requires a compiler supporting **C++23** and CMake ≥ 3.28.

**Step 1** — Create and activate the development environment:

.. code-block:: console

   $ conda env create -f kep3_devel.yml
   $ conda activate kep3_devel

This installs all mandatory C++ and Python dependencies into an isolated conda
environment.

**Step 2** — Configure and build:

.. code-block:: console

   $ cmake -S . -B build -G Ninja              \
       -DCMAKE_EXPORT_COMPILE_COMMANDS=1        \
       -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"   \
       -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"      \
       -Dkep3_BUILD_PYTHON_BINDINGS=ON
   $ cmake --build build --target install --parallel

**Step 3** — Verify the installation:

.. code-block:: console

   $ python -c "import pykep; print(pykep.__version__)"

.. rubric:: CMake flags reference

+-------------------------------------------+----------------------------------------------------------+
| Flag                                      | Description                                              |
+===========================================+==========================================================+
| ``-DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"``| Install into the active conda environment.               |
+-------------------------------------------+----------------------------------------------------------+
| ``-DCMAKE_PREFIX_PATH="$CONDA_PREFIX"``   | Resolve dependencies from the active conda environment.  |
+-------------------------------------------+----------------------------------------------------------+
| ``-Dkep3_BUILD_CPP_LIBRARY=ON``           | Build the ``kep3`` C++ library from source (default).    |
|                                           | Set to ``OFF`` to use an already-installed ``kep3``      |
|                                           | located via ``find_package`` instead.                    |
+-------------------------------------------+----------------------------------------------------------+
| ``-Dkep3_BUILD_PYTHON_BINDINGS=ON``       | Build and install the ``pykep`` Python extension module. |
+-------------------------------------------+----------------------------------------------------------+
| ``-Dkep3_BUILD_TESTS=ON``                 | Build the C++ unit-test suite.                           |
+-------------------------------------------+----------------------------------------------------------+
| ``-Dkep3_BUILD_BENCHMARKS=ON``            | Build the benchmark executables.                         |
+-------------------------------------------+----------------------------------------------------------+
| ``-DKEP3_VERBOSE_CONFIGURE=ON``           | Emit additional configure-time diagnostics.              |
+-------------------------------------------+----------------------------------------------------------+
| ``-DCMAKE_BUILD_TYPE=Debug``              | Enable debug symbols (single-config generators only).    |
+-------------------------------------------+----------------------------------------------------------+

.. _installation_verify:

Verifying the installation
==========================

After installation, confirm pykep is working correctly:

.. code-block:: console

   $ python -c "import pykep; print(pykep.__version__)"

To run the Python test suite:

.. code-block:: console

   $ python -m pytest /path/to/pykep/tests

.. rubric:: Running the C++ test suite

To verify the C++ library itself, configure a build with ``kep3_BUILD_TESTS=ON``
and run the tests via CTest:

.. code-block:: console

   $ cmake -S . -B build -G Ninja              \
       -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"      \
       -Dkep3_BUILD_TESTS=ON
   $ cmake --build build --parallel
   $ ctest --test-dir build --output-on-failure

.. note::

   The ``-DCMAKE_INSTALL_PREFIX`` flag is not required when building tests
   only — nothing is installed. CTest discovers the test executables directly
   from the build tree.

Getting help
============

If you encounter problems during installation, please open an issue on the
`GitHub issue tracker <https://github.com/esa/pykep/issues>`_.
