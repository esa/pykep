# Copyright (C) 2004-2009 The PaGMO development team,
# Advanced Concepts Team (ACT), European Space Agency (ESA)
# http://apps.sourceforge.net/mediawiki/pagmo
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits
# act@esa.int
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

INCLUDE(FindPythonLibs)
# We need the Python interpreter to figure out Python's version in SuckOSX.
INCLUDE(FindPythonInterp)

# Find Python libraries
FIND_PACKAGE(PythonLibs REQUIRED)
MESSAGE(STATUS "Python libraries: " "${PYTHON_LIBRARIES}")
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
MESSAGE(STATUS "Python library: " "${PYTHON_LIBRARY}")

# These flags are used to signal the need to override the default extension of the Python modules
# depending on the architecture. Under Windows, for instance, CMake produces shared objects as
# .dll files, but Python from 2.5 onwards requires .pyd files (hence the need to override).
# A similar thing happens in SuckOSX.
SET(PYDEXTENSION FALSE)
SET(SOEXTENSION FALSE)

IF(UNIX)
	# We need the Python interpreter in order to detect the appropriate directory of modules installation.
	IF(NOT PYTHONINTERP_FOUND)
		MESSAGE(FATAL_ERROR "Unable to locate Python interpreter.")
	ENDIF(NOT PYTHONINTERP_FOUND)
	MESSAGE(STATUS "Python interpreter is: ${PYTHON_EXECUTABLE}")
	# Now we must establish if the installation dir for Python modules is named 'site-packages' (as usual)
	# or 'dist-packages' (apparently Ubuntu 9.04 or maybe Python 2.6, it's not clear).
        EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/PyKEP/python_packages_dir.py
		OUTPUT_VARIABLE PY_PACKAGES_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
	MESSAGE(STATUS "Python packages dir is: ${PY_PACKAGES_DIR}")
	# SuckOSX suckages.
	IF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
		MESSAGE(STATUS "OSX system detected.")
		SET(SOEXTENSION TRUE)
		# Let's determine Python version by running the interpreter with the --version flag.
		EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} --version ERROR_VARIABLE PY_VERSION_OSX ERROR_STRIP_TRAILING_WHITESPACE)
		MESSAGE(STATUS "Python interpeter returns string: " ${PY_VERSION_OSX})
		STRING(REGEX MATCH [0-9]*\\.[0-9]* PYTHON_LIBRARY_VERSION_DOT ${PY_VERSION_OSX})
	ELSE(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
		# In sane Unix system we can fetch the Python version number directly from the library.
		STRING(REGEX MATCH libpython[0-9]*\\.[0-9]* PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY})
	ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	# Remove the dot from the Python version.
	STRING(REGEX REPLACE libpython "" PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY_VERSION_DOT})
	STRING(REGEX REPLACE \\. "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION_DOT})
	# Let's use CMAKE_INSTALL_PREFIX, so that if we specify a different install path it will be respected.
	SET(PYTHON_MODULES_PATH lib/python${PYTHON_LIBRARY_VERSION_DOT}/${PY_PACKAGES_DIR})
ELSE(UNIX)
	STRING(REGEX MATCH python[0-9]* PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY})
	STRING(REGEX REPLACE python "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION})
	SET(PYTHON_MODULES_PATH .)
	IF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
		MESSAGE(STATUS "Python >= 2.5 detected on WIN32 platform. Output extension for compiled modules will be '.pyd'.")
		SET(PYDEXTENSION TRUE)
	ENDIF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
ENDIF(UNIX)
MESSAGE(STATUS "Python library version: " ${PYTHON_LIBRARY_VERSION})
MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_PATH}")
