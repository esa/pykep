.. _prepareyoursystem:

Preparing your system to compile pykep
======================================

Compiling under a Unix System
-----------------------------

As pykep has some dependencies, you need to install a few packages before being able to install and compile pykep on your machine. 
All are supported by most package managers.

* Install `git <http://git-scm.com>`_
* Install `CMake <http://www.cmake.org>`_ with its ccmake utility
* Install `boost libraries <http://www.boost.org>`_ both headers and compiled libraries

Compiling under Windows
-----------------------

Whoa, Whoa why would you want to do such a crazy crazy thing???? Anyway you need to:

* Install `Python 3.X <https://www.python.org/>`_ 
* Install `MinGW (gcc compiler for windows) <https://sourceforge.net/projects/mingw-w64/files/mingw-w64/mingw-w64-release/>`_ 
* Install `CMake <http://www.cmake.org>`_ 
* Install `boost libraries <http://www.boost.org>`_ both headers and compiled libraries

Alternatively to minGW also clang-cl is tested in our CI.