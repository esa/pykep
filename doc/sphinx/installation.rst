.. _howtoinstall:

Install pykep
======================

pykep supports 32 and 64 bits architectures and both python 2.7 and python 3.X We are unfortunately not in the position to provide
binaries for everybody, so it is likely you will need to compile the source code youself. This process, under unix architectures is
painless, while under windows architectures it can be troublesome, but we do provide windows 64 bits binaries via PyPi

Using Windows Binaries
----------------------

If you have a windows 64-bits based python (27, 34, 35 or 36) you can install pykep via pip typing::

  pip install pykep

Compiling and Installing under Unix
-----------------------------------

Assuming you have prepared your system for compiling pykep (see :ref:`prepareyoursystem`) and that you have just downloaded the source code following the instructions given, see :ref:`howtodownload`, you will have
created a directory keptoolbox in your current directory, move there::

  cd pykep

You will now need to create a build directory where to build the source code, so::

  mkdir build

You can now move there::

  cd build

and have ccmake help you select the options that are most suitable for you::

  ccmake ../

At this point (after pressing c once to configure) you should be seeing something like this on the screen:

.. image:: images/ccmake.*

You can now press 'g' to generate a make file and exit ccmake utility. You are back to the prompt where you can now type::

  make

and::

  sudo make install

Watch carefully the message in the terminal where the installation path is given to check
that the correct python dist-packages or site-packages directory has been located

Here is a typical example of the output obtained::

  [ 91%] Built target propagate_lagrangian_test
  [ 92%] Built target propagate_lagrangian_u_test
  [ 94%] Built target propagate_taylor_jorba_test
  [ 96%] Built target propagate_taylor_s_test
  [ 98%] Built target propagate_taylor_test
  [100%] Built target sgp4_test
  Install the project...
  -- Install configuration: ""
  -- Installing: /usr/local/lib/libkeplerian_toolbox.dylib
  -- Up-to-date: /usr/local/lib/python2.7/site-packages/pykep/__init__.py
  -- Installing: /usr/local/lib/python2.7/site-packages/pykep/core/_core.so

Compiling and Installing under Windows
--------------------------------------

Unsing minGW things will be the same as under Unix, just make sure that

* You have compiled the boost libraries correctly (i.e invoking bjam with the option toolset=gcc link=shared).
* Place the whole boost directory where the CMake script can find it (e.g. in C:/boost). This may also require renaming the folder from boost_x_xx_xx to boost)
* Check, when running CMake, that all libraries are found correctly
* When running a make install, Windows will probably put your pykep directory under Program Files/kep_toolbox, move it to the correct place (e.g. C:/Python27/Lib/site-packages/)
* Put all dll (boost and keplerian_toolbox) in pykep/core
* Hope for the best (kidding its super easy ...)

Systems with both Python 2 and Python 3 installed
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
