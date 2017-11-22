# CMake version check.
cmake_minimum_required(VERSION 3.2)

# Main pykep project version.
set(PYKEP_PROJECT_VERSION 2.1)
project(pykep VERSION ${PYKEP_PROJECT_VERSION})
enable_testing()

# Set default build type to "Release".
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

#Some information to screen
message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
message(STATUS "CXX Compiler detected: ${CMAKE_CXX_COMPILER_ID}")

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
MESSAGE(STATUS "CMake additional search path for libraries: ${CMAKE_LIBRARY_PATH}")

# Writing some configuration files with version update ...
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pykep/__init__.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/pykep/__init__.py" @ONLY)

# We set compiler specific flags
INCLUDE(cmake_scripts/set_compiler_flags.cmake)

# Provides build options to CMake
# Build Option: when active the file main.cpp is built and linked to the pagmo static library
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Build Option: when active builds Python bindings called pykep
OPTION(BUILD_PYKEP "Build pykep." OFF)

# Build Option: when active the examples in the folder example are compiled into executables
OPTION(BUILD_TESTS "Build tests'." ON)

# Build Option: when active the JPL SPICE toolbox is used to extend pykep. You can choose
# to build the library or use a pre-build library
OPTION(USE_SPICE "Interface pykep to JPL SPICE toolbox'." OFF)
OPTION(BUILD_SPICE "Interface pykep to JPL SPICE toolbox'." OFF)

# Check consistency.
if(USE_SPICE AND BUILD_SPICE)
message(FATAL_ERROR "Please select whether to build spice or just try to find a library")
endif()

# Build Option: when active headers will be installed in the CMAKE_INSTALL_PREFIX/include directory
OPTION(INSTALL_HEADERS "Installs the header files" OFF)

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
# Minimum version must have std::array serialization
FIND_PACKAGE(Boost 1.56.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
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

# We build cspice static library and link it in pykep libraries
IF(BUILD_SPICE)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src/third_party/cspice")
	SET(LIBRARIES_4_STATIC ${LIBRARIES_4_STATIC} cspice)
    SET(LIBRARIES_4_DYNAMIC ${LIBRARIES_4_DYNAMIC} cspice)
    ADD_DEFINITIONS(-DPYKEP_USING_SPICE)
ENDIF(BUILD_SPICE)

# Here we define the file list and make the static and dynamic library
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

# Build main and link it to static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main keplerian_toolbox_static ${LIBRARIES_4_STATIC})
ENDIF(BUILD_MAIN)

# Build pykep and link it to dynamic library.
IF(BUILD_PYKEP)
	SET(LIBRARIES_4_PYTHON_MODULES keplerian_toolbox ${LIBRARIES_4_PYTHON_MODULES})
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/pykep")
	if(MINGW)
		message(STATUS "Creating the files for the generation of a binary wheel for MinGW.")
		configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py")
		configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt")
	else(MINGW)
        message(STATUS "Creating the files for the generation of a manylinux binary wheel.")
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/manylinux_wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py")
    endif()
ENDIF(BUILD_PYKEP)

# Build Tests and link them to static library.
IF(BUILD_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
	FILE(COPY "${CMAKE_SOURCE_DIR}/tests/data/sgp4_test.txt" DESTINATION "${CMAKE_BINARY_DIR}/tests")
	FILE(COPY "${CMAKE_SOURCE_DIR}/tests/data/C_G_1000012_2012_2017.bsp" DESTINATION "${CMAKE_BINARY_DIR}/tests")
	FILE(COPY "${CMAKE_SOURCE_DIR}/tests/data/pck00010.tpc" DESTINATION "${CMAKE_BINARY_DIR}/tests")
	FILE(COPY "${CMAKE_SOURCE_DIR}/tests/data/gm_de431.tpc" DESTINATION "${CMAKE_BINARY_DIR}/tests")
ENDIF(BUILD_TESTS)
