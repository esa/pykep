if(UNIX)
	# On unix platforms, we set the rpath of the pykep libraries
	# to the lib subdirectory of the install prefix. This allows
	# to install pykep in nonstandard locations (e.g., ~/.local)
	# without having to muck around with library search paths.
	set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# The name of the BP target to link to depends on the Boost and Python versions.
# NOTE: since Boost 1.67, the naming of the Boost.Python library has changed to include the
# major and minor python version as a suffix. See the release notes:
# https://www.boost.org/users/history/version_1_67_0.html
if(${Boost_MAJOR_VERSION} GREATER 1 OR (${Boost_MAJOR_VERSION} EQUAL 1 AND ${Boost_MINOR_VERSION} GREATER 66))
    set(PYKEP_BP_TARGET "Boost::python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
else()
    if(${PYTHON_VERSION_MAJOR} EQUAL 2)
        set(PYKEP_BP_TARGET "Boost::python")
    else()
        set(PYKEP_BP_TARGET "Boost::python3")
    endif()
endif()

set(PYKEP_BOOST_TARGET "keplerian_toolbox" "Boost::disable_autolinking" "Boost::date_time" "Boost::serialization")

# Boost.Python version 1.64 (and possibly later?) has issues with default hidden visibility.
if(${Boost_MAJOR_VERSION} GREATER 1 OR (${Boost_MAJOR_VERSION} EQUAL 1 AND ${Boost_MINOR_VERSION} GREATER 63))
	list(FIND PYKEP_CXX_FLAGS_DEBUG "-fvisibility=hidden" _PYKEP_HAS_VHIDDEN)
	if (NOT ${_PYKEP_HAS_VHIDDEN} EQUAL -1)
		message(STATUS "Removing the '-fvisibility=hidden' flag for Boost.Python > 1.63.")
		list(REMOVE_ITEM PYKEP_CXX_FLAGS_DEBUG "-fvisibility=hidden")
		list(REMOVE_ITEM PYKEP_CXX_FLAGS_RELEASE "-fvisibility=hidden")
	endif()
endif()

# Setup the installation path.
set(PYKEP_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/pykep")


ADD_SUBDIRECTORY(core)
ADD_SUBDIRECTORY(planet)
ADD_SUBDIRECTORY(sims_flanagan)
ADD_SUBDIRECTORY(util)
ADD_SUBDIRECTORY(pontryagin)
ADD_SUBDIRECTORY(orbit_plots)
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(trajopt)
ADD_SUBDIRECTORY(phasing)
INSTALL(FILES __init__.py test.py DESTINATION ${PYKEP_INSTALL_PATH})

# Definitions needed for windows ... (using dlls)
ADD_DEFINITIONS(-DKEP_TOOL_DLL_EXPORT_API)




