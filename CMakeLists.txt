PROJECT(Keplerian_Toolbox)

ENABLE_TESTING()

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0)
# Setup and include CMake modules.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../../cmake_modules)

# Configuration for GCC.
IF(CMAKE_COMPILER_IS_GNUCXX)
	# Add to the base flags extra warnings, and set up message reporting from the compiler so that it is Eclipse-friendly.
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wdisabled-optimization -fmessage-length=0")
	# Suggested for multithreaded code.
	ADD_DEFINITIONS(-D_REENTRANT)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# Initialise (empty) list of libraries to link against.
SET(MANDATORY_LIBRARIES "")

#

# Provides build options to CMake
# Build Option: when active the file main.cpp is built and linked to the PaGMO static library
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Build Option: when active builds Python bindings called PyKEP
OPTION(BUILD_PYKEP "Build PyKEP." OFF)

# Build Option: when active the examples in the folder example are compiled into executables
OPTION(BUILD_TESTS "Build tests'." ON)

# Build Option: when active the serialization module is activated
OPTION(ENABLE_SERIALIZATION "Support for serialization" OFF)

# Keplerian_toolbox lib source files.
SET(SRC_LIST
        ${CMAKE_CURRENT_SOURCE_DIR}/src/epoch.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/planet.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/lambert_problem.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/asteroid_gtoc2.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/asteroid_gtoc5.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/planet_mpcorb.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/planet_ss.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sims_flanagan/leg.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sims_flanagan/fb_traj.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/core_functions/jorba.c
	)

# If serialization support is requested look for the boost library and link to it
IF(ENABLE_SERIALIZATION)
        FIND_LIBRARY(SERIALIZATION_LIBRARY NAMES boost_serialization)
        IF(NOT SERIALIZATION_LIBRARY)
                MESSAGE(FATAL_ERROR "Serialization support was enabled but no libboost_serialization found.")
        ENDIF(NOT SERIALIZATION_LIBRARY)
        MESSAGE(STATUS "Boost serialization library: ${SERIALIZATION_LIBRARY}")
        SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${SERIALIZATION_LIBRARY})
        ADD_DEFINITIONS(-DKEP_TOOLBOX_ENABLE_SERIALIZATION)
ENDIF(ENABLE_SERIALIZATION)

FIND_LIBRARY(BOOST_DATE_TIME_LIBRARY NAMES boost_date_time)
IF(NOT BOOST_DATE_TIME_LIBRARY)
                MESSAGE(FATAL_ERROR "Boost date_time library is required but libboost_date_time not found.")
ENDIF(NOT BOOST_DATE_TIME_LIBRARY)
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${BOOST_DATE_TIME_LIBRARY})

ADD_LIBRARY(keplerian_toolbox_static STATIC ${SRC_LIST})

# Link main to pagmo_static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
        TARGET_LINK_LIBRARIES(main keplerian_toolbox_static ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN)

IF(BUILD_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
ENDIF(BUILD_TESTS)

IF(BUILD_PYKEP)
        SET(LIB_INSTALL_PATH "lib")
        INCLUDE(FindPythonLibs)
        FIND_LIBRARY(BOOST_PYTHON_LIBRARY NAMES boost_python)
        SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${BOOST_PYTHON_LIBRARY})
        ADD_LIBRARY(keplerian_toolbox SHARED ${SRC_LIST})
        # Request linking against the Python libs, in addition to our mandatory libs.
        TARGET_LINK_LIBRARIES(keplerian_toolbox ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARY})
        # Add the definitions needed when building the PaGMO DLL.
        INSTALL(TARGETS keplerian_toolbox
                RUNTIME DESTINATION ${LIB_INSTALL_PATH}
                LIBRARY DESTINATION ${LIB_INSTALL_PATH}
        )
        ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/PyKEP")
ENDIF(BUILD_PYKEP)

