# Lets find the python interpreter (FindPythonInterp actually does a decent job)
INCLUDE(FindPythonInterp)
IF(PYTHONINTERP_FOUND)
    MESSAGE(STATUS "Python interpreter: ${PYTHON_EXECUTABLE}")
	MESSAGE(STATUS "Version detected for the python interpreter: ${PYTHON_VERSION_STRING}")
ELSE(PYTHONINTERP_FOUND)
	MESSAGE(FATAL_ERROR "Unable to locate Python Interpreter.")
ENDIF(PYTHONINTERP_FOUND)

# We now use the python executable to tell us (via the distutil package) all relevant libarry information
# As the CMake script FindPythonLibs sucks biiiig time
UNSET(PYTHON_INCLUDE_DIR CACHE)
UNSET(PYTHON_MODULES_DIR CACHE)
UNSET(PYTHON_LIBRARIES CACHE)

EXECUTE_PROCESS ( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())" OUTPUT_VARIABLE PYTHON_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
EXECUTE_PROCESS ( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())" OUTPUT_VARIABLE PYTHON_MODULES_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
EXECUTE_PROCESS ( COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_SOURCE_DIR}/cmake_scripts/find_python_library.py" OUTPUT_VARIABLE PYTHON_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})
MESSAGE(STATUS "Python includes path: "        "${PYTHON_INCLUDE_DIR}")
MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_DIR}")
MESSAGE(STATUS "Python library name: "         "${PYTHON_LIBRARIES}")

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

MARK_AS_ADVANCED(CLEAR PYTHON_EXECUTABLE)
# We make the variables cached so that we can change them in GUIs editors like ccmake 
SET(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE} CACHE PATH "The Python executable")
SET(PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR} CACHE PATH "Path to the python include files, where pyconfig.h can be found")
SET(PYTHON_MODULES_DIR ${PYTHON_MODULES_DIR} CACHE PATH "Path of site-packages, PyKEP will be installed in .../PyKEP")
SET(PYTHON_LIBRARIES ${PYTHON_LIBRARIES} CACHE PATH "Name of the python library")
