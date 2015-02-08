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

#Some screen output
MESSAGE(STATUS "CMake Library additional search path: ${CMAKE_LIBRARY_PATH}")

# Setting compiler secific flags
INCLUDE(CheckCXXCompilerFlag)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

	CHECK_CXX_COMPILER_FLAG(-ftemplate-depth=256 CLANG_TEMPLATE_DEPTH)
	IF(CLANG_TEMPLATE_DEPTH)
		MESSAGE(STATUS "Enabling '-ftemplate-depth=256' compiler flag required since boost 1.54.")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=256")
	ENDIF(CLANG_TEMPLATE_DEPTH)

ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

	CHECK_CXX_COMPILER_FLAG(fno-strict-aliasing GNUCXX_NO_STRICT_ALIASING)
	IF(GNUCXX_NO_STRICT_ALIASING)
		MESSAGE(STATUS "Enabling '-fno-strict-aliasing' compiler flag.")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing")
	ENDIF(GNUCXX_NO_STRICT_ALIASING)

	CHECK_CXX_COMPILER_FLAG(-Wnon-virtual-dtor GNUCXX_W_ALL)
	IF(GNUCXX_W_ALL)
		MESSAGE(STATUS "Enabling '-Wall' compiler flag.")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
	ENDIF(GNUCXX_W_ALL)

	CHECK_CXX_COMPILER_FLAG(-Wno-deprecated GNUCXX_W_EXTRA)
	IF(GNUCXX_W_EXTRA)
		MESSAGE(STATUS "Enabling '-Wextra' compiler flag.")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
	ENDIF(GNUCXX_W_EXTRA)

	CHECK_CXX_COMPILER_FLAG(-Wno-ignored-qualifiers GNUCXX_W_NOEXCEPT)
	IF(GNUCXX_W_NOEXCEPT)
		MESSAGE(STATUS "Enabling '-Wnoexcept' compiler flag.")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnoexcept")
	ENDIF(GNUCXX_W_NOEXCEPT)

	CHECK_CXX_COMPILER_FLAG(-Wno-ignored-qualifiers GNUCXX_W_DISOPT)
	IF(GNUCXX_W_DISOPT)
		MESSAGE(STATUS "Enabling '-Wdisabled-optimization' compiler flag.")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wdisabled-optimization")
	ENDIF(GNUCXX_W_DISOPT)

	CHECK_CXX_COMPILER_FLAG(-Wno-ignored-qualifiers GNUCXX_FASTMATH)
	IF(GNUCXX_FASTMATH)
		MESSAGE(STATUS "Enabling '-ffast-math' compiler flag.")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
	ENDIF(GNUCXX_FASTMATH)

ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
	CHECK_CXX_COMPILER_FLAG(-mieee-fp INTEL_IEEE_COMPLIANT)
	IF(INTEL_IEEE_COMPLIANT)
			MESSAGE(STATUS "Enabling '-mieee-fp' compiler flag to get IEEE compliant code")
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mieee-fp")
	ENDIF(INTEL_IEEE_COMPLIANT)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
ENDIF()

# The following flags will be tried on all compilers
# If a C++11 able compiler is detected, then the flag is added
# otherwise an error is thrown
CHECK_CXX_COMPILER_FLAG(-std=c++11 ALL_C11)
IF(ALL_C11)
	MESSAGE(STATUS "Enabling '-std=c++11' compiler flag")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF(ALL_C11)

IF(NOT ALL_C11)
	CHECK_CXX_COMPILER_FLAG(-std=c++0x ALL_C0X)
	IF(ALL_C0X)
		MESSAGE(STATUS "Enabling '-std=c++0x' compiler flag")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	ELSE(ALL_C0X)
	    MESSAGE(FATAL_ERROR "Unable to locate a c++11 compatible compiler")
	ENDIF(ALL_C0X)
ENDIF(NOT ALL_C11)

# We finally set and log the compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "CXX compilation flags: ${CMAKE_CXX_FLAGS}")

# Initialise (empty) list of libraries to link against.
SET(LIBRARIES_4_KEP_TOOLBOX "")
SET(LIBRARIES_4_PYKEP "")

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
	SET(REQUIRED_BOOST_LIBS ${REQUIRED_BOOST_LIBS} python)
ENDIF(BUILD_PYKEP)
MESSAGE(STATUS "Required Boost libraries: ${REQUIRED_BOOST_LIBS}")
FIND_PACKAGE(Boost 1.42.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
MESSAGE(STATUS "Detected Boost version: ${Boost_VERSION}")
MESSAGE(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
SET(LIBRARIES_4_KEP_TOOLBOX ${LIBRARIES_4_KEP_TOOLBOX} ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})
SET(LIBRARIES_4_PYKEP ${LIBRARIES_4_PYKEP} ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} ${Boost_PYTHON_LIBRARY})

# Here we define the file list
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

# Build main and link it to static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main keplerian_toolbox_static ${MANDATORY_LIBRARIES_4_KEP_TOOLBOX})
ENDIF(BUILD_MAIN)

# Build PyKEP and link it to dynamic library.
IF(BUILD_PYKEP)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/PyKEP")
ENDIF(BUILD_PYKEP)

# Build Tests and link them to static library.
IF(BUILD_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
	IF(ENABLE_SGP4)
		FILE(COPY "${CMAKE_SOURCE_DIR}/tests/sgp4_test.txt" DESTINATION "${CMAKE_BINARY_DIR}/tests")
	ENDIF(ENABLE_SGP4)
ENDIF(BUILD_TESTS)

