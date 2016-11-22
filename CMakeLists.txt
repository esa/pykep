cmake_minimum_required(VERSION 3.0.0)

PROJECT(Keplerian_Toolbox VERSION 1.3.0)
# We enable cmake test capability
ENABLE_TESTING()

# We set the minimum required cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Set default build type to "Release", change it in the GUI if you need to build with debug.
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

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
MESSAGE(STATUS "CMake additional search path for libraries: ${CMAKE_LIBRARY_PATH}")

# We set compiler specific flags
INCLUDE(cmake_scripts/set_compiler_flags.cmake)

# Provides build options to CMake
# Build Option: when active the file main.cpp is built and linked to the PaGMO static library
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Build Option: when active builds Python bindings called PyKEP
OPTION(BUILD_PYKEP "Build PyKEP." OFF)

# Build Option: when active the examples in the folder example are compiled into executables
OPTION(BUILD_TESTS "Build tests'." ON)

# Build Option: when active the JPL SPICE toolbox is used to extend PyKEP
OPTION(BUILD_SPICE "Interface PyKEP to JPL SPICE toolbox'." OFF)

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
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
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

# We build cspice static library and link it in PyKEP libraries
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

# Build PyKEP and link it to dynamic library.
IF(BUILD_PYKEP)
	SET(LIBRARIES_4_PYTHON_MODULES keplerian_toolbox ${LIBRARIES_4_PYTHON_MODULES})
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/PyKEP")
	if(MINGW)
		message(STATUS "Creating the files for the generation of a binary wheel for MinGW.")
		configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py")
		configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt")
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
