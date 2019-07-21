.. _howtoinstall:

Install pykep
======================

pykep supports 32 and 64 bits architectures and python 3.X. Python 2.7 support was dropped as of version 2.4.
Both PyPi and conda package mangers contain the binaries for pykep latest code, but only for a limited set of
architectures. In case yours is not include, you will have to compile the code from source.


Using Binaries (encouraged whenever possible)
----------------------------------------------

If you have a compatible architecture you can install pykep via pip/conda typing::

  pip install pykep

or, if you are using conda::

  conda install pykep

In case both fail, you probably do not have an architecture which we support binaries for.

Compiling and Installing under Linux (degree of pain: low)
------------------------------------------------------------------

Assuming you have prepared your system for compiling pykep (see :ref:`prepareyoursystem`) and that you have just downloaded the source code
following the instructions given, see :ref:`howtodownload`, you will have created a directory keptoolbox in your current directory, 
move there::

  cd pykep

You will now need to create a build directory where to build the source code, so::

  mkdir build

You can now move there::

  cd build

and have ccmake help you select the options that are most suitable for you::

  ccmake ../

At this point (after pressing c once to configure and having selected the correct options) you should be seeing something like this on the screen:

.. image:: images/ccmake.png

Note that in the case above only the tests are built, not the python module nor the headers will be installed. Also note that the 
```CMAKE_INSTALL_PREFIX``` points to the ```.local``` folder of the user. 
This correspond to a (suggested) local installation of **pykep**.

You can now press 'g' to generate a make file and exit ccmake utility. 
You are back to the prompt where you can now type::

  make

and::

  sudo make install

Watch carefully the message in the terminal where the installation path is given to check
that the correct python dist-packages or site-packages directory has been located

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

You can now activate, in ccmake, the build option ```BUILD_PYKEP``` and compile/install the python module, or ```INSTALL_HEADER``` and install the headers.

.. note::

   Check carefully what boost python library is selected automatically by cmake, and if needed change it.

Compiling and Installing under Windows (degree of pain: high)
------------------------------------------------------------------

Unsing minGW things should be roughly the same as under Unix, just make sure that

* You have compiled the boost libraries correctly (i.e invoking bjam with the option toolset=gcc link=shared).
* Place the whole boost directory where the CMake script can find it (e.g. in C:/boost). This may also require renaming the folder from boost_x_xx_xx to boost)
* Check, when running CMake, that all libraries are found correctly
* When running a make install, Windows will probably put your pykep directory under Program Files/kep_toolbox, move it to the correct place (e.g. C:/PythonXX/Lib/site-packages/)
* Put all dll (boost and keplerian_toolbox) in pykep/core
* Hope for the best (kidding its super easy ...)
* No, really hope for the best

Systems with multiple python versions
-------------------------------------------------

If your system has several versions of python installed check the PYTHON_EXECUTABLE variable in cmake. The libraries, includes and site-packages directory are determined accordingly. If you want to change python version, just define explicitly such a variable. For example (assuming you are in a directory pykep/build)::

  cmake ../ -DBUILD_PYKEP="ON"  -DPYTHON_EXECUTABLE="/usr/bin/python3.3m"

It is always good practice to check what cmake has actually located by typing::

  cmake ../

which could look something like::

  -- OS detected: Darwin
  -- CXX Compiler detected: Clang
  -- CMake additional search path for libraries: /usr/local/lib
  -- Enabling '-ftemplate-depth=256' compiler flag required since boost 1.54.
  -- Enabling '-std=c++11' compiler flag
  -- CXX compilation flags:  -ftemplate-depth=256 -std=c++11
  -- Python interpreter: /usr/local/bin/python3
  -- Python interpreter verison: 3.4
  -- Python includes path: /usr/local/Cellar/python3/3.4.2_1/Frameworks/Python.framework/Versions/3.4/include/python3.4m
  -- Python modules install path: /usr/local/Cellar/python3/3.4.2_1/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages
  -- Python library name: /usr/local/Cellar/python3/3.4.2_1/Frameworks/Python.framework/Versions/3.4/lib/libpython3.4.dylib
  -- Required Boost libraries: serialization;date_time;python3
  -- Boost version: 1.57.0
  -- Found the following Boost libraries:
  --   serialization
  --   date_time
  --   python3
  -- Detected Boost version: 105700
  -- Boost include dirs: /usr/local/include
  -- Boost libraries: /usr/local/lib/libboost_serialization-mt.dylib;/usr/local/lib/libboost_date_time-mt.dylib;/usr/local/lib/libboost_python3-mt.dylib
  -- Configuring done
  -- Generating done
  -- Build files have been written to: /Users/darioizzo/Documents/pykep/build
