if(UNIX)
# On unix platforms, we set the rpath of the dcgpy libraries
# to the lib subdirectory of the install prefix. This allows
# to install dcgpy in nonstandard locations (e.g., ~/.local)
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

# Boost.Python version 1.64 (and possibly later?) has issues with default hidden visibility.
if(${Boost_MAJOR_VERSION} GREATER 1 OR (${Boost_MAJOR_VERSION} EQUAL 1 AND ${Boost_MINOR_VERSION} GREATER 63))
list(FIND KEP_TOOLBOX_CXX_FLAGS_DEBUG "-fvisibility=hidden" _PYKEP_HAS_VHIDDEN)
if (NOT ${_PYKEP_HAS_VHIDDEN} EQUAL -1)
	message(STATUS "Removing the '-fvisibility=hidden' flag for Boost.Python > 1.63.")
	list(REMOVE_ITEM KEP_TOOLBOX_CXX_FLAGS_DEBUG "-fvisibility=hidden")
	list(REMOVE_ITEM KEP_TOOLBOX_CXX_FLAGS_RELEASE "-fvisibility=hidden")
endif()
endif()

# Setup of the header-only pykep library.
add_library(pykep INTERFACE)
target_link_libraries(pykep INTERFACE ${PYKEP_BP_TARGET} Boost::disable_autolinking Keplerian_toolbox::keplerian_toolbox)
target_include_directories(pykep INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/pykep> $<INSTALL_INTERFACE:include>)

# Extra flags to be appended in CI experiments
set(PYKEP_CXX_FLAGS_EXTRA)
mark_as_advanced(PYKEP_CXX_FLAGS_EXTRA)
list(APPEND DCGP_CXX_FLAGS_RELEASE ${PYKEP_CXX_FLAGS_EXTRA})

# Setting the install path
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





