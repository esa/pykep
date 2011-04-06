.. _howtoinstall:

Install PyKEP
======================

Using Windows Binaries
----------------------

Download the `binaries <http://sourceforge.net/projects/keptoolbox/files/Windows/PyKEP-1.0_python27_mgw.rar/download>`_
and extract the whole folder as it is in your python site-packages directory 
(e.g. extract to C:/Python27/Lib/site-packages/)

Compiling and Installing under Unix
-----------------------------------

Assuming you have just downloaded the source code following the instructions given, see :ref:`howtodownload`, you will have 
created a directory keptoolbox in your current directory, move there::

  cd keptoolbox

You will now need to create a build directory where to build the source code, so::

  mkdir build

You can now move there::

  cd build

and have ccmake help you select the options that are most suitable for you::

  ccmake ../

At this point (after pressing c once to configure) you should be seeing something like this on the screen:

.. image:: images/ccmake.*

You need to activate the option BUILD_PYKEP. As for the other options, just ignore them, they are useful only for c++ developers / users. 
Upon pressing 'c' again, cmake tries to locate the boost python library necessary to build the code and, in case it is successfull,
displays something like

.. image:: images/ccmake2.*

You can now press 'g' to generate a make file and exit ccmake utility. You are back to the prompt where you can now type::

  make

and::

  sudo make install

Watch carefully the message in the terminal where the installation path is given to check 
that the correct python dist-packages or site-packages directory has been located

Here is a typical example of the output obtained (gentoo system)::
  [ 47%] Built target keplerian_toolbox

  [ 94%] Built target keplerian_toolbox_static

  [100%] Built target _PyKEP

  Install the project...

  -- Install configuration: "Release"

  -- Installing: /usr/local/lib/libkeplerian_toolbox.so

  -- Installing: /usr/local/lib/python2.6/site-packages/PyKEP/__init__.py

  -- Installing: /usr/local/lib/python2.6/site-packages/PyKEP/_PyKEP.so

  -- Removed runtime path from "/usr/local/lib/python2.6/site-packages/PyKEP/_PyKEP.so"

Compiling and Installing under Windows
--------------------------------------

Same as under Unix, just make sure that

* You have compiled the boost libraries correctly (i.e invoking bjam with the option toolset=gcc link=shared). 
* Place the whole boost directory where the CMake script can find it (e.g. in C:/boost). This may also require renaming the folder from boost_x_xx_xx to boost)
* Check, when running CMake, that all libraries are found correctly
* When running a make install, Windows will probably put your PyKEP directory under Program Files/kep_toolbox, move it to the correct place (e.g. C:/Python27/Lib/site-packages/)
* Put all dll (boost and keplerian_toolbox) in PyKEP/core
* Hope for the best (honestly, if it works for you just download the binaries.... it is easier)