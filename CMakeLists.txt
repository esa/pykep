PROJECT(Keplerian_Toolbox)

# List of directories to search for CMake modules.
# Commands like include() and find_package() search for files in directories listed
# by this variable before checking the default modules that come with CMake.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)

# We enable cmake test capability
ENABLE_TESTING()

# We set the minimum required cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12.2)

#Some screen output
MESSAGE(STATUS "OS detected: ${CMAKE_SYSTEM_NAME}")
MESSAGE(STATUS "CXX Compiler detected: ${CMAKE_CXX_COMPILER_ID}")

# Setting OS specific cmake variables
# 1 - OSX
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	SET(CMAKE_MACOSX_RPATH OFF)
	SET(CMAKE_OSX_ROOT "")
	SET(CMAKE_OSX_DEPLOYMENT_TARGET "")
	# Homebrew directory added to path
	SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# 2 - NONE

# Some screen output
MESSAGE(STATUS "CMake Library additional search path: ${CMAKE_LIBRARY_PATH}")

# We set compiler specific flags
INCLUDE(cmake_scripts/set_compiler_flags.cmake)

# Provides build options to CMake
# Build Option: when active the file main.cpp is built and linked to the PaGMO static library
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Build Option: when active builds Python bindings called PyKEP
OPTION(BUILD_PYKEP "Build PyKEP." OFF)

# Build Option: when active the examples in the folder example are compiled into executables
OPTION(BUILD_TESTS "Build tests'." ON)

# Build Option: when active libsgp4 is also built and related classes activated
OPTION(ENABLE_SGP4 "Build SGP4 related classes'." OFF)
IF(ENABLE_SGP4)
		ADD_DEFINITIONS(-DPYGMO_ENABLE_SGP4)
ENDIF(ENABLE_SGP4)

SET(DYNAMIC_LIB_KEP_TOOL_USE_FLAGS "-DBOOST_THREAD_USE_DLL -DBOOST_SERIALIZATION_DYN_LINK=1")

# Finding the boost libraries needed for the keplerian_toolbox
SET(REQUIRED_BOOST_LIBS serialization date_time)
IF(BUILD_PYKEP)
	INCLUDE(cmake_scripts/python_setup.cmake)
	IF(PYTHON_VERSION_MAJOR LESS 3)
		SET(REQUIRED_BOOST_LIBS ${REQUIRED_BOOST_LIBS} python)
	ELSE(PYTHON_VERSION_MAJOR LESS 3)
		SET(REQUIRED_BOOST_LIBS ${REQUIRED_BOOST_LIBS} python3)
	ENDIF(PYTHON_VERSION_MAJOR LESS 3)
ENDIF(BUILD_PYKEP)
MESSAGE(STATUS "Required Boost libraries: ${REQUIRED_BOOST_LIBS}")
FIND_PACKAGE(Boost 1.48.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
MESSAGE(STATUS "Detected Boost version: ${Boost_VERSION}")
MESSAGE(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# Initialise (empty) list of libraries to link against.
SET(LIBRARIES_4_STATIC "")
SET(LIBRARIES_4_DYNAMIC "")
SET(LIBRARIES_4_PYTHON_MODULES "")

# Define the libraries to link against.
SET(LIBRARIES_4_STATIC ${LIBRARIES_4_STATIC} ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})
SET(LIBRARIES_4_DYNAMIC ${LIBRARIES_4_DYNAMIC} ${LIBRARIES_4_STATIC})
IF(PYTHON_VERSION_MAJOR LESS 3)
	SET(LIBRARIES_4_PYTHON_MODULES ${LIBRARIES_4_PYTHON_MODULES} ${Boost_PYTHON_LIBRARY}  ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})
ELSE(PYTHON_VERSION_MAJOR LESS 3)
	SET(LIBRARIES_4_PYTHON_MODULES ${LIBRARIES_4_PYTHON_MODULES} ${Boost_PYTHON3_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})
ENDIF(PYTHON_VERSION_MAJOR LESS 3)

# Here we define the file list and make the static and dynamic library
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

# Build main and link it to static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main keplerian_toolbox_static ${LIBRARIES_4_STATIC})
ENDIF(BUILD_MAIN)

# Build PyKEP and link it to dynamic library.
IF(BUILD_PYKEP)
	SET(LIBRARIES_4_PYTHON_MODULES keplerian_toolbox ${LIBRARIES_4_PYTHON_MODULES})
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/PyKEP")
ENDIF(BUILD_PYKEP)

# Build Tests and link them to static library.
IF(BUILD_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
	IF(ENABLE_SGP4)
		FILE(COPY "${CMAKE_SOURCE_DIR}/tests/sgp4_test.txt" DESTINATION "${CMAKE_BINARY_DIR}/tests")
	ENDIF(ENABLE_SGP4)
ENDIF(BUILD_TESTS)

