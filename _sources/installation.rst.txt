.. _howtoinstall:

Install pykep
======================

pykep supports python 3.X. Python 2.7 support was dropped as of version 2.4.
Both PyPi and conda package mangers contain the binaries for pykep, but only for a limited set of
architectures. In case yours is not include, you will have to compile the code from source.


Using Binaries 
--------------

The Python module called ``pykep`` can be installed from conda or pip only for some architectures / python version combinations.

Installing with conda
^^^^^^^^^^^^^^^^^^^^^
``pykep`` is available in the `conda <https://conda.io/en/latest/>`__ package manager
from the `conda-forge <https://conda-forge.org/>`__ channel. A single package is available:

* `pykep <https://anaconda.org/conda-forge/pykep>`__, which contains the ``pykep`` python module.

In order to install ``pykep`` via conda, you just need
to add ``conda-forge`` to the channels:

.. code-block:: console

   $ conda config --add channels conda-forge
   $ conda install pykep

Please refer to the `conda documentation <https://conda.io/en/latest/>`__ for instructions
on how to setup and manage your conda installation.

Installing with pip
^^^^^^^^^^^^^^^^^^^
We also provide the pip packages (mainly for linux 64 bit architectures). Check on the 
`PyPi pykep page <https://pypi.org/project/pykep/>`_ if the needed package is provided.

.. code-block:: console

   $ pip install pykep

Compiling and Installing under Linux
------------------------------------

Assuming you have prepared your system for compiling pykep (see :ref:`prepareyoursystem`) and that you have just downloaded the source code
following the instructions given, see :ref:`howtodownload`, you will have a directory pykep with the code, move there::

  cd pykep

You will now need to 1) build and install the keplerian toolbox (i.e. the c++ library) and then build and install pykep (the python module
that depends on the keplerian toolbox library). So lets start::

  mkdir build
  cd build

and, using cmake::

  cmake -DBoost_NO_BOOST_CMAKE=ON \
        -DPYKEP_BUILD_KEP_TOOLBOX=yes \
        -DPYKEP_BUILD_PYKEP=no \
        -DPYKEP_BUILD_SPICE=yes \
        -DPYKEP_BUILD_TESTS=no \
        -DCMAKE_BUILD_TYPE=Release ../;
  make install

Here is a typical example of the output obtained::

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

You can run the tests now typing::

  make test

And you should see something like::

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

Now you need to compile and install the pykep module::

  cd /pykep
  mkdir build_pykep
  cd build_pykep
  cmake -DBoost_NO_BOOST_CMAKE=ON \
        -DPYKEP_BUILD_KEP_TOOLBOX=no \
        -DPYKEP_BUILD_PYKEP=yes \
        -DPYKEP_BUILD_TESTS=no \
        -DCMAKE_BUILD_TYPE=Release \
        -DBoost_PYTHON${PYTHON_VERSION}_LIBRARY_RELEASE=/usr/local/lib/${BOOST_PYTHON_LIBRARY_NAME} \
        -DPYTHON_EXECUTABLE=/opt/python/${PYTHON_DIR}/bin/python ../;
  make -j2 install

Watch carefully the message in the terminal where the installation path is given to check
that the correct python dist-packages or site-packages directory has been located

Compiling and Installing under Windows 
------------------------------------------------------------------

We do not really support nor reccomend doing this, but in case you are really motivated, you can get inspired by our
`azure CI script <https://github.com/esa/pykep/blob/master/azure-pipelines.yml>`_ that works using the
`conda <https://conda.io/en/latest/>`_ package manager or our
`appveyor CI script <https://github.com/esa/pykep/blob/master/tools/install_appveyor_mingw.py>`_  which
makes use on minGW. 