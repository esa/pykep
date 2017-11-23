# CMake version check.
cmake_minimum_required(VERSION 3.2)

# Main pykep project version.
set(PYKEP_PROJECT_VERSION 2.1)
project(pykep VERSION ${PYKEP_PROJECT_VERSION})
enable_testing()

# Module path setup.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules" "${CMAKE_SOURCE_DIR}/cmake_modules/yacma")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Set default build type to "Release".
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# Writing some configuration files with version update ...
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pykep/__init__.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/pykep/__init__.py" @ONLY)

# Initial setup of compiler flags.
include(YACMACompilerLinkerSettings)

# Provides build options to CMake
# Build Option: when active the file main.cpp is built and linked to the PYKEP static library
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Build Option: when active builds Python bindings called pykep
OPTION(BUILD_PYKEP "Build pykep." OFF)

# Build Option: when active the examples in the folder example are compiled into executables
OPTION(BUILD_TESTS "Build tests'." ON)

# Build Option: when active the JPL SPICE toolbox is used to extend pykep. You can choose
# to build the library or use a pre-build library
OPTION(BUILD_SPICE "Interface pykep to JPL SPICE toolbox'." OFF)

# Build Option: when active headers will be installed in the CMAKE_INSTALL_PREFIX/include directory
OPTION(INSTALL_HEADERS "Installs the header files" OFF)

# Assemble the flags.
set(PYKEP_CXX_FLAGS_DEBUG ${YACMA_CXX_FLAGS} ${YACMA_CXX_FLAGS_DEBUG} ${YACMA_THREADING_CXX_FLAGS})
set(PYKEP_CXX_FLAGS_RELEASE ${YACMA_CXX_FLAGS} ${YACMA_THREADING_CXX_FLAGS})
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND YACMA_COMPILER_IS_CLANGXX)
  message(STATUS "Clang compiler on OSX detected, setting the standard library to 'libc++'.")
  list(APPEND PYKEP_CXX_FLAGS_DEBUG "-stdlib=libc++")
  list(APPEND PYKEP_CXX_FLAGS_RELEASE "-stdlib=libc++")
endif()
if(YACMA_COMPILER_IS_MSVC)
  # Disable the idiotic minmax macros on MSVC, some annoying warnings,
  # and enable the bigobj option.
  list(APPEND PYKEP_CXX_FLAGS_DEBUG "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj")
  list(APPEND PYKEP_CXX_FLAGS_RELEASE "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj")
endif()
if(YACMA_COMPILER_IS_INTELXX)
  # NOTE: on MSVC we use the push/pop pragmas, but they do not seem to work on Intel (the pragmas
  # in icc influence the behaviour at instantiation point, not at definition point).
  # These warnings are useful in principle, but they are generated a lot from cereal and we have no
  # way of disabling them selectively. Just rely on the other compilers to provde good diagnostic.
  list(APPEND PYKEP_CXX_FLAGS_DEBUG "-diag-disable" "2259,1682,68")
  list(APPEND PYKEP_CXX_FLAGS_RELEASE "-diag-disable" "2259,1682,68")
endif()
if(MINGW)
	# Flag needed to deal with big binaries in MinGW.
	message(STATUS "Enabling the '-Wa,-mbig-obj' flag in MinGW builds.")
	list(APPEND PYKEP_CXX_FLAGS_DEBUG "-Wa,-mbig-obj")
	list(APPEND PYKEP_CXX_FLAGS_RELEASE "-Wa,-mbig-obj")
endif()

if(BUILD_PYKEP)
# pygmo dependencies.
include(YACMAPythonSetup)

# Python version check.
if(${PYTHON_VERSION_MAJOR} LESS 2 OR (${PYTHON_VERSION_MAJOR} EQUAL 2 AND ${PYTHON_VERSION_MINOR} LESS 7))
	message(FATAL_ERROR "Minimum supported Python version is 2.7.")
endif()

# Internal variable that will be used to tell PYKEPFindBoost to locate Boost.Python.
set(_PYKEP_FIND_BOOST_PYTHON TRUE)
endif()

# Boost setup. Will allow using Boost::library_name for the needed components
include(PykepFindBoost)

# We build cspice static library and link it in pykep libraries
if(BUILD_SPICE)
add_subdirectory("${CMAKE_SOURCE_DIR}/src/third_party/cspice")
    add_definitions(-DPYKEP_USING_SPICE)
endif(BUILD_SPICE)

# Here we define the file list and make the static and dynamic library
add_subdirectory("${CMAKE_SOURCE_DIR}/src")

# Build main and link it to static library.
if(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main keplerian_toolbox_static Boost::boost Boost::serialization Boost::date_time)
	endif(BUILD_MAIN)

# Build Tests and link them to static library.
if(BUILD_TESTS)
	add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
	file(COPY "${CMAKE_SOURCE_DIR}/tests/data/sgp4_test.txt" DESTINATION "${CMAKE_BINARY_DIR}/tests")
	file(COPY "${CMAKE_SOURCE_DIR}/tests/data/C_G_1000012_2012_2017.bsp" DESTINATION "${CMAKE_BINARY_DIR}/tests")
	file(COPY "${CMAKE_SOURCE_DIR}/tests/data/pck00010.tpc" DESTINATION "${CMAKE_BINARY_DIR}/tests")
	FILE(COPY "${CMAKE_SOURCE_DIR}/tests/data/gm_de431.tpc" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endif(BUILD_TESTS)

# Build pykep and link it to dynamic library.
if(BUILD_PYKEP)
	add_subdirectory("${CMAKE_SOURCE_DIR}/pykep")
	if(MINGW)
		message(STATUS "Creating the files for the generation of a binary wheel for MinGW.")
		configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py")
		configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt")
	else(MINGW)
        message(STATUS "Creating the files for the generation of a manylinux binary wheel.")
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/manylinux_wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py")
	endif()
endif()



