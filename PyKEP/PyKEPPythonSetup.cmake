# Lets find the python libraries (and Python.h)
INCLUDE(FindPythonLibs)
IF(PYTHONLIBS_FOUND)
	INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
    MESSAGE(STATUS "Path to the python library: ${PYTHON_LIBRARIES}")
	MESSAGE(STATUS "Path to where Python.h is found: ${PYTHON_INCLUDE_PATH}")
	MESSAGE(STATUS "Version detected for python libraries: ${PYTHONLIBS_VERSION_STRING}")
ELSE(PYTHONLIBS_FOUND)
	MESSAGE(FATAL_ERROR "Unable to locate Python libraries.")
ENDIF(PYTHONLIBS_FOUND)

# Lets find the python interpreter
INCLUDE(FindPythonInterp)
IF(PYTHONINTERP_FOUND)
    MESSAGE(STATUS "Python interpreter: ${PYTHON_EXECUTABLE}")
	MESSAGE(STATUS "Version detected for the python interpreter: ${PYTHON_VERSION_STRING}")
ELSE(PYTHONINTERP_FOUND)
	MESSAGE(FATAL_ERROR "Unable to locate Python Interpreter.")
ENDIF(PYTHONINTERP_FOUND)

# These flags are used to signal the need to override the default extension of the Python modules
# depending on the architecture. Under Windows, for instance, CMake produces shared objects as
# .dll files, but Python from 2.5 onwards requires .pyd files (hence the need to override).
# A similar thing happens in SuckOSX.
SET(PYDEXTENSION FALSE)
SET(SOEXTENSION FALSE)

IF(UNIX)
	# SuckOSX suckages.
	IF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
		SET(SOEXTENSION TRUE)
	ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	ELSE(UNIX)
	# Win32 suckages.
	SET(PYDEXTENSION TRUE)
ENDIF(UNIX)

# We now construct the path where we will install all files by determining the system wide path for site-packages,
# which is expected to be something like /usr/blah/blah/blah/lib/python2.7/site-packages. We then take the last 3
# subdir to create lib/python2.7/site-packages (or dist-utils or whatever python version is correct)
execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print('/'.join(get_python_lib().split('/')[-3:]))" OUTPUT_VARIABLE PYTHON_MODULES_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_PATH}")
