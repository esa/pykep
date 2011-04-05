PROJECT(Keplerian_Toolbox)

ENABLE_TESTING()

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0)
# Setup and include CMake modules.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../../cmake_modules)

 IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING 
			"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Configuration for GCC.
IF(CMAKE_COMPILER_IS_GNUCXX)
	# Add to the base flags extra warnings, and set up message reporting from the compiler so that it is Eclipse-friendly.
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wdisabled-optimization -fmessage-length=0 -ffast-math")
	# This flag is suggested for performance. One can check with the tests what difference it really makes
	#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
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

SET(DYNAMIC_LIB_KEP_TOOL_USE_FLAGS "-DBOOST_THREAD_USE_DLL -DBOOST_SERIALIZATION_DYN_LINK=1")

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
	${CMAKE_CURRENT_SOURCE_DIR}/src/sims_flanagan/spacecraft.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/core_functions/jorba.c
	)

# Finding the boost libraries needed for the keplerian_toolbox
set(Boost_ADDITIONAL_VERSIONS "1.46.1")

SET(REQUIRED_BOOST_LIBS serialization date_time)
IF(BUILD_PYKEP)
	SET(REQUIRED_BOOST_LIBS ${REQUIRED_BOOST_LIBS} python)
ENDIF(BUILD_PYKEP)
MESSAGE(STATUS "Required Boost libraries: ${REQUIRED_BOOST_LIBS}")

FIND_PACKAGE(Boost 1.42.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
MESSAGE(STATUS "Detected Boost version: ${Boost_VERSION}")
MESSAGE(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

SET(MANDATORY_BOOST_LIBS ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})
IF(BUILD_PYKEP)
	SET(MANDATORY_BOOST_LIBS ${MANDATORY_BOOST_LIBS} ${Boost_PYTHON_LIBRARY})
ENDIF(BUILD_PYKEP)
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${MANDATORY_BOOST_LIBS})
	
#Build Static Library
ADD_LIBRARY(keplerian_toolbox_static STATIC ${SRC_LIST})

# Build main and link it to pagmo_static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
        TARGET_LINK_LIBRARIES(main keplerian_toolbox_static ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN)

# Build Tests
IF(BUILD_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
ENDIF(BUILD_TESTS)

#Build PyKEP
IF(BUILD_PYKEP)
        SET(LIB_INSTALL_PATH "lib")
        INCLUDE(FindPythonLibs)
        ADD_LIBRARY(keplerian_toolbox SHARED ${SRC_LIST})
        # Request linking against the Python libs, in addition to our mandatory libs.
        TARGET_LINK_LIBRARIES(keplerian_toolbox ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARY})
        # Add the definitions needed when building the PaGMO DLL.
		SET_TARGET_PROPERTIES(keplerian_toolbox PROPERTIES COMPILE_FLAGS "-DKEP_TOOL_DLL_EXPORT_API")
        INSTALL(TARGETS keplerian_toolbox
                RUNTIME DESTINATION ${LIB_INSTALL_PATH}
                LIBRARY DESTINATION ${LIB_INSTALL_PATH}
        )
        ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/PyKEP")
ENDIF(BUILD_PYKEP)

#ADD_LIBRARY(keplerian_toolbox SHARED ${SRC_LIST})
 #       # Request linking against the Python libs, in addition to our mandatory libs.
 #       TARGET_LINK_LIBRARIES(keplerian_toolbox ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARY})
  #      # Add the definitions needed when building the PaGMO DLL.
	#SET_TARGET_PROPERTIES(keplerian_toolbox PROPERTIES COMPILE_FLAGS "-DKEP_TOOL_DLL_EXPORT_API")

#ADD_EXECUTABLE(mainlol main.cpp)
#SET_TARGET_PROPERTIES(mainlol PROPERTIES COMPILE_FLAGS "-DKEP_TOOL_DLL_IMPORT_API ${DYNAMIC_LIB_KEP_TOOL_USE_FLAGS}")
#TARGET_LINK_LIBRARIES(mainlol keplerian_toolbox ${MANDATORY_LIBRARIES})
