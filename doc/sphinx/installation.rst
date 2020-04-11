.. _howtoinstall:

Installation
======================

Dependencies
------------

pykep has the following **mandatory** runtime dependencies:

* `Python <https://www.python.org/>`__ 3.4 or later (Python 2.x is
  **not** supported),
* the `keplerian_toolbox C++ library <https://esa.github.io/pykep/>`__, same version as pykep included in the pykep source,
* the `Boost serialization library <https://www.boost.org/doc/libs/release/libs/serialization/doc/index.html>`__,
  version 1.60 or later,

Additionally, pykep has the following **optional** runtime
dependencies:

* `scipy <https://www.scipy.org/>`__, which is used in many plotting and related utilities.
* `Matplotlib <https://matplotlib.org/>`__, which is used in many plotting and related utilities.
  plotting utilities,
* `scikit-learn <https://scikit-learn.org/stable//>`__, which is used, for example, in the phasing module
* `numba <http://numba.pydata.org/>`__, which is used to speed up some pure python code.


Packages
--------

pykep packages are available from a variety
of package managers on several platforms.

Conda
^^^^^

pykep is available via the `conda <https://conda.io/docs/>`__
package manager for Linux, OSX and Windows
thanks to the infrastructure provided by `conda-forge <https://conda-forge.org/>`__.
In order to install pykep via conda, you just need to add ``conda-forge``
to the channels, and then you can immediately install pykep:

.. code-block:: console

   $ conda config --add channels conda-forge
   $ conda config --set channel_priority strict
   $ conda install pykep

The conda packages for pykep are maintained by the core development team,
and they are regularly updated when new pykep versions are released.

Please refer to the `conda documentation <https://conda.io/docs/>`__
for instructions on how to setup and manage
your conda installation.

pip
^^^

pykep is also available on Linux via the `pip <https://pip.pypa.io/en/stable/>`__
package installer. The installation of pykep with pip is straightforward:

.. code-block:: console

   $ pip install pykep

If you want to install pykep for a single user instead of
system-wide, which is in general a good idea, you can do:

.. code-block:: console

   $ pip install --user pykep

An even better idea is to make use of Python's
`virtual environments <https://virtualenv.pypa.io/en/latest/>`__.

The pip packages for pykep are maintained by the core development team,
and they are regularly updated when new pykep versions are released.

.. warning::

   Note however that we **strongly** encourage users to install pykep with conda
   rather than with pip. The reason is that pykep is built on a
   moderately complicated
   stack of C++ libraries, which have to be bundled together with pykep
   in the pip package.
   This is a problem if one uses pykep together with other Python
   packages sharing dependencies with pykep, because multiple incompatible
   versions of the same C++ library might end up being loaded at the
   same time, leading to crashes and erratic runtime behaviour.
   The conda packages do not suffer from this issue.

Installation from source
------------------------

In order to install pygmo from source, you will need:

* a C++11 capable compiler (any recent version of GCC,
  Clang or MSVC should do),
* a `Python <https://www.python.org/>`__ installation,
* the `Boost libraries <https://www.boost.org/>`__,
* `CMake <https://cmake.org/>`__, version 3.3 or later.

After making sure the above dependencies are installed on your system, you can
download the pykep source code from the
`GitHub release page <https://github.com/esa/pykep/releases>`__. Alternatively,
and if you like living on the bleeding edge, you can get the very latest
version of pykep via ``git``:

.. code-block:: console

   $ git clone https://github.com/esa/pykep.git

We follow the usual PR-based development workflow, thus pykep's ``master``
branch is normally kept in a working state.

After downloading and/or unpacking pykep's
source code, go to pykep's
source tree, create a ``build`` directory and ``cd`` into it. E.g.,
on a Unix-like system:

.. code-block:: console

   $ cd /path/to/pykep
   $ mkdir build
   $ cd build

Once you are in the ``build`` directory, you must configure your build
using ``cmake``. There are various useful CMake variables you can set,
such as:

* ``CMAKE_BUILD_TYPE``: the build type (``Release``, ``Debug``, etc.),
  defaults to ``Release``.
* ``CMAKE_INSTALL_PREFIX``: the path into which pygmo will be installed
  (e.g., this defaults to ``/usr/local`` on Unix-like platforms).
* ``CMAKE_PREFIX_PATH``: additional paths that will be searched by CMake
  when looking for dependencies.

Please consult `CMake's documentation <https://cmake.org/cmake/help/latest/>`_
for more details about CMake's variables and options.

Before compiling pykep we need to install the keplerian toolbox, that is the C++
code at the heart of most computationally expensive routines. In a linux installation,
this would typically look like:

.. code-block:: console

    $ cmake -DBoost_NO_BOOST_CMAKE=ON \
        -DPYKEP_BUILD_KEP_TOOLBOX=yes \
        -DPYKEP_BUILD_PYKEP=no \
        -DPYKEP_BUILD_SPICE=yes \
        -DPYKEP_BUILD_TESTS=yes \
        -DCMAKE_INSTALL_PREFIX = ~/.local \
        -DCMAKE_PREFIX_PATH = ~/.local \
        -DCMAKE_BUILD_TYPE=Release \
        ../;
    $ cmake  --build . --target install
 
Here is a typical example of the output you should expect::

  [ 99%] Built target propagate_lagrangian_u_test
  [ 99%] Built target anomalies_test
  [ 99%] Built target lambert_test
  [ 99%] Built target propagate_taylor_test
  [ 99%] Built target load_spice_kernel_test
  [ 99%] Built target propagate_lagrangian_test
  [ 99%] Built target propagate_taylor_s_test 
  [100%] Built target spice_planet_test
  [100%] Built target sgp4_test
  Install the project...
  -- Install configuration: "Release"

You can also run the C++ tests typing:

.. code-block:: console

     $ make test

And see, hopefully, something like the following output::

  Running tests...
  Test project /home/dario/Develop/pykep/build
      Start  1: lambert_test
  1/12 Test  #1: lambert_test .....................   Passed    1.19 sec
      Start  2: propagate_lagrangian_test
  2/12 Test  #2: propagate_lagrangian_test ........   Passed    0.15 sec
      Start  3: propagate_lagrangian_u_test
  3/12 Test  #3: propagate_lagrangian_u_test ......   Passed    0.20 sec
      Start  4: propagate_taylor_test
  4/12 Test  #4: propagate_taylor_test ............   Passed    0.19 sec
      Start  5: propagate_taylor_J2_test
  5/12 Test  #5: propagate_taylor_J2_test .........   Passed    0.33 sec
      Start  6: propagate_taylor_jorba_test
  6/12 Test  #6: propagate_taylor_jorba_test ......   Passed    0.22 sec
      Start  7: propagate_taylor_s_test
  7/12 Test  #7: propagate_taylor_s_test ..........   Passed    0.22 sec
      Start  8: leg_s_test
  8/12 Test  #8: leg_s_test .......................   Passed    0.00 sec
      Start  9: sgp4_test
  9/12 Test  #9: sgp4_test ........................   Passed    0.12 sec
      Start 10: anomalies_test
  10/12 Test #10: anomalies_test ...................   Passed    0.56 sec
      Start 11: load_spice_kernel_test
  11/12 Test #11: load_spice_kernel_test ...........   Passed    0.04 sec
      Start 12: spice_planet_test
  12/12 Test #12: spice_planet_test ................   Passed    0.01 sec

  100% tests passed, 0 tests failed out of 12

  Total Test time (real) =   3.24 sec

You are now ready to compile and install pykep. Move back one level from the folder where you installed the 
keplerian_toolbox and create a new build directory:

.. code-block:: console

     $ cd ..
     $ mkdir build_pykep
     $ cd build_pykep

There we need to run cmake again with slightly different options:

.. code-block:: console

     $ cmake -DBoost_NO_BOOST_CMAKE=ON \
         -DPYKEP_BUILD_KEP_TOOLBOX=no \
         -DPYKEP_BUILD_PYKEP=yes \
         -DPYKEP_BUILD_TESTS=no \
         -DCMAKE_INSTALL_PREFIX = ~/.local \
         - DCMAKE_PREFIX_PATH = ~/.local \
         - DCMAKE_BUILD_TYPE=Release \
         ../;

Watch carefully the message in the terminal where the installation path is given to check
that the correct python dist-packages or site-packages directory has been located. If not, 
set explicitly the relevant cmake variables, for example ```Boost_PYTHON37_LIBRARY_RELEASE``` 
and ```PYTHON_EXECUTABLE```.

Verifying the installation
--------------------------

You can verify that pykep was successfully compiled and
installed by running the test suite. From a
Python session, run the following commands:

.. code-block:: python

   >>> import pykep
   >>> pykep.test.run_test_suite()

If these commands execute without any error, then
your pykep installation is ready for use.

