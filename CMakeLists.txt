# Main pykep project version.
set(PYKEP_PROJECT_VERSION 2.4)

# CMake version check.
cmake_minimum_required(VERSION 3.2)

# Module path setup.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules" "${CMAKE_SOURCE_DIR}/cmake_modules/yacma")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Set default build type to "Release".
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# Main build options: build keplerian_toolbox or pykep. They cannot be on at the same time,
# and only one must be chosen.
option(PYKEP_BUILD_KEP_TOOLBOX "Build the cpp project." ON)
option(PYKEP_BUILD_PYKEP "Build pykep." OFF)

# Check consistency.
if(PYKEP_BUILD_KEP_TOOLBOX AND PYKEP_BUILD_PYKEP)
    message(FATAL_ERROR "Please select whether to build keplerian_toolbox or pykep: you cannot build them both at the same time.")
endif()

if((NOT PYKEP_BUILD_KEP_TOOLBOX) AND (NOT PYKEP_BUILD_PYKEP))
    message(FATAL_ERROR "Please select if you want to build the cpp project or pykep")
endif()

if(PYKEP_BUILD_KEP_TOOLBOX)
    # Initial setup of a keplerian_toolbox build.
    project(keplerian_toolbox VERSION ${PYKEP_PROJECT_VERSION})
    enable_testing()

    # Build option: enable test set.
    option(PYKEP_BUILD_TESTS "Build test set." ON)

    # Build option: enable examples
    option(PYKEP_BUILD_SPICE "Build also the spice tools." ON)
else()
    # Initial setup of a pykep build.
    project(pykep VERSION ${PYKEP_PROJECT_VERSION})
endif()

# Common general bits.

# Initial setup of compiler flags.
include(YACMACompilerLinkerSettings)

# Threading setup.
include(YACMAThreadingSetup)

# Assemble the flags.
set(KEP_TOOLBOX_CXX_FLAGS_DEBUG ${YACMA_CXX_FLAGS} ${YACMA_CXX_FLAGS_DEBUG} ${YACMA_THREADING_CXX_FLAGS})
set(KEP_TOOLBOX_CXX_FLAGS_RELEASE ${YACMA_CXX_FLAGS} ${YACMA_THREADING_CXX_FLAGS})
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND YACMA_COMPILER_IS_CLANGXX)
  message(STATUS "Clang compiler on OSX detected, setting the standard library to 'libc++'.")
  list(APPEND KEP_TOOLBOX_CXX_FLAGS_DEBUG "-stdlib=libc++")
  list(APPEND KEP_TOOLBOX_CXX_FLAGS_RELEASE "-stdlib=libc++")
endif()
if(YACMA_COMPILER_IS_MSVC)
  include(CheckCXXCompilerFlag)
  # Disable the idiotic minmax macros on MSVC, some annoying warnings,
  # and enable the bigobj option.
  # Also, enable the WIN32_LEAN_AND_MEAN definition.
  list(APPEND KEP_TOOLBOX_CXX_FLAGS_DEBUG "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj" "-DWIN32_LEAN_AND_MEAN")
  list(APPEND KEP_TOOLBOX_CXX_FLAGS_RELEASE "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj" "-DWIN32_LEAN_AND_MEAN")
  # Enable strict conformance mode, if supported.
  set(CMAKE_REQUIRED_QUIET TRUE)
  check_cxx_compiler_flag("/permissive-" _KEP_TOOLBOX_MSVC_SUPPORTS_STRICT_CONFORMANCE)
  unset(CMAKE_REQUIRED_QUIET)
  if(_KEP_TOOLBOX_MSVC_SUPPORTS_STRICT_CONFORMANCE)
    message(STATUS "The '/permissive-' flag is supported, enabling it.")
    list(APPEND KEP_TOOLBOX_CXX_FLAGS_DEBUG "/permissive-")
    list(APPEND KEP_TOOLBOX_CXX_FLAGS_RELEASE "/permissive-")
  endif()
  unset(_KEP_TOOLBOX_MSVC_SUPPORTS_STRICT_CONFORMANCE)
  if(YACMA_COMPILER_IS_CLANGXX)
    # clang-cl emits various warnings from GMP/MPFR, let's just silence them.
    # NOTE: at one point in the recent past, MSVC added an options similar to GCC's isystem:
    # https://blogs.msdn.microsoft.com/vcblog/2017/12/13/broken-warnings-theory/
    # We probably just need to wait for this to be picked up by CMake/clang-cl. Let's
    # revisit the issue in the future.
    list(APPEND _KEP_TOOLBOX_CLANG_CL_DISABLED_WARNINGS
      "-Wno-unused-variable"
      "-Wno-inconsistent-dllimport"
      "-Wno-unknown-pragmas"
      "-Wno-unused-parameter"
      "-Wno-sign-compare"
      "-Wno-deprecated-declarations"
      "-Wno-deprecated-dynamic-exception-spec"
      "-Wno-old-style-cast"
      "-Wno-sign-conversion"
      "-Wno-non-virtual-dtor"
      "-Wno-deprecated"
      "-Wno-shadow"
      "-Wno-shorten-64-to-32"
      "-Wno-reserved-id-macro"
      "-Wno-undef"
      "-Wno-c++98-compat-pedantic"
      "-Wno-documentation-unknown-command"
      "-Wno-zero-as-null-pointer-constant"
      "-Wno-language-extension-token"
      "-Wno-gnu-anonymous-struct"
      "-Wno-nested-anon-types"
      "-Wno-documentation"
      "-Wno-comma"
      "-Wno-nonportable-system-include-path"
      "-Wno-global-constructors"
      "-Wno-redundant-parens"
      "-Wno-exit-time-destructors"
      "-Wno-missing-noreturn"
      "-Wno-switch-enum"
      "-Wno-covered-switch-default"
      "-Wno-float-equal"
      "-Wno-double-promotion"
      "-Wno-microsoft-enum-value"
      "-Wno-missing-prototypes"
      "-Wno-implicit-fallthrough"
      "-Wno-format-nonliteral"
      "-Wno-cast-qual"
      "-Wno-disabled-macro-expansion"
      "-Wno-unused-private-field"
      "-Wno-unused-template"
      "-Wno-unused-macros"
      "-Wno-extra-semi-stmt"
      "-Wno-c++98-compat")
    list(APPEND KEP_TOOLBOX_CXX_FLAGS_DEBUG ${_KEP_TOOLBOX_CLANG_CL_DISABLED_WARNINGS})
    list(APPEND KEP_TOOLBOX_CXX_FLAGS_RELEASE ${_KEP_TOOLBOX_CLANG_CL_DISABLED_WARNINGS})
    unset(_KEP_TOOLBOX_CLANG_CL_DISABLED_WARNINGS)
  endif()
endif()
if(YACMA_COMPILER_IS_INTELXX)
  # NOTE: on MSVC we use the push/pop pragmas, but they do not seem to work on Intel (the pragmas
  # in icc influence the behaviour at instantiation point, not at definition point).
  # These warnings are useful in principle, but they are generated a lot from cereal and we have no
  # way of disabling them selectively. Just rely on the other compilers to provde good diagnostic.
  list(APPEND KEP_TOOLBOX_CXX_FLAGS_DEBUG "-diag-disable" "2259,1682,68")
  list(APPEND KEP_TOOLBOX_CXX_FLAGS_RELEASE "-diag-disable" "2259,1682,68")
endif()
if(MINGW)
	# Flag needed to deal with big binaries in MinGW.
	message(STATUS "Enabling the '-Wa,-mbig-obj' flag in MinGW builds.")
	list(APPEND KEP_TOOLBOX_CXX_FLAGS_DEBUG "-Wa,-mbig-obj")
	list(APPEND KEP_TOOLBOX_CXX_FLAGS_RELEASE "-Wa,-mbig-obj")
endif()

# Some flags that generate warnings due to Eigen obsolete versions
list(REMOVE_ITEM KEP_TOOLBOX_CXX_FLAGS_DEBUG "-Wduplicated-branches")
list(REMOVE_ITEM KEP_TOOLBOX_CXX_FLAGS_DEBUG "-Wold-style-cast")

if(PYKEP_BUILD_PYKEP)
    include(YACMAPythonSetup)
    # Python version check.
    if(${PYTHON_VERSION_MAJOR} LESS 2 OR (${PYTHON_VERSION_MAJOR} EQUAL 2 AND ${PYTHON_VERSION_MINOR} LESS 7))
        message(FATAL_ERROR "Minimum supported Python version is 2.7.")
    endif()

    # NOTE: for the time being, require that dcgp/dcgpy versions are matching exactly.
    find_package(keplerian_toolbox ${PYKEP_PROJECT_VERSION} EXACT REQUIRED)

    # Internal variable that will be used to tell PykepFindBoost to locate Boost.Python.
    set(_PYKEP_FIND_BOOST_PYTHON TRUE)
endif()

# Boost setup. Will allow using Boost::library_name for the needed components
include(PykepFindBoost)

if(PYKEP_BUILD_KEP_TOOLBOX)
    # List of source files.
    SET(KEP_TOOLBOX_SRC_FILES
        # Core
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epoch.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lambert_problem.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/sims_flanagan/leg.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/sims_flanagan/leg_s.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/sims_flanagan/spacecraft.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/core_functions/jorba.c"
        # Planet
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/base.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/keplerian.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/j2.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/jpl_low_precision.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/mpcorb.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/gtoc2.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/gtoc5.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/gtoc6.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/gtoc7.cpp"
    )
    # We keep these in a separate list as to be able to have different compile flags
    SET(LIBSGP4_SRC_FILES
        # PYKEP FILE
        "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/tle.cpp"
        # SGP4 FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/libsgp4/Util.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/libsgp4/Tle.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/libsgp4/SolarPosition.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/libsgp4/SGP4.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/libsgp4/OrbitalElements.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/libsgp4/Observer.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/libsgp4/Eci.cpp"
    )
    # We cannot use c++ strict compiling flags on a c file. The following ad hoc property is thus set to 
    # silence resulting warnings
    if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
        set_source_files_properties(${KEP_TOOLBOX_SRC_FILES} PROPERTIES COMPILE_OPTIONS "${KEP_TOOLBOX_CXX_FLAGS_RELEASE}")
    elseif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        set_source_files_properties(${KEP_TOOLBOX_SRC_FILES} PROPERTIES COMPILE_OPTIONS "${KEP_TOOLBOX_CXX_FLAGS_DEBUG}")
    endif()
    set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/src/core_functions/jorba.c" PROPERTIES COMPILE_OPTIONS -w)

    if(PYKEP_BUILD_SPICE)
        # Add cpp files to the keplerian toolbox
        set(KEP_TOOLBOX_SRC_FILES
            "${CMAKE_CURRENT_SOURCE_DIR}/src/planet/spice.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/util/spice_utils.cpp"
            "${KEP_TOOLBOX_SRC_FILES}")
        # Add c files to the keplerian toolbox
        file(GLOB CSPICE_LIB_SRC_LIST "src/third_party/cspice/*.c")
        # Silence warnings
        set_source_files_properties(${CSPICE_LIB_SRC_LIST} PROPERTIES COMPILE_OPTIONS -w)
        set(PYKEP_USING_SPICE "#define PYKEP_USING_SPICE")
    endif(PYKEP_BUILD_SPICE)

    # Creates the file config.hpp
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/include/keplerian_toolbox/config.hpp" @ONLY)

    # Setup of the keplerian_toolbox library.
    add_library(keplerian_toolbox SHARED "${CSPICE_LIB_SRC_LIST}" "${KEP_TOOLBOX_SRC_FILES}" "${LIBSGP4_SRC_FILES}" )
    set_property(TARGET keplerian_toolbox PROPERTY VERSION "1.0")
    set_property(TARGET keplerian_toolbox PROPERTY SOVERSION 1)
    set_target_properties(keplerian_toolbox PROPERTIES CXX_VISIBILITY_PRESET hidden)
    set_target_properties(keplerian_toolbox PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)
    # Let's setup the target C++ standard, but only if the user did not provide it manually.
    if(NOT CMAKE_CXX_STANDARD)
        set_property(TARGET keplerian_toolbox PROPERTY CXX_STANDARD 11)
    endif()
    set_property(TARGET keplerian_toolbox PROPERTY CXX_STANDARD_REQUIRED YES)
    set_property(TARGET keplerian_toolbox PROPERTY CXX_EXTENSIONS NO)
    # NOTE: make sure the include directories from the current build
    # are included first, so that if there is already a pagmo installation
    # in the prefix path we don't risk including the headers from that
    # one instead.
    target_include_directories(keplerian_toolbox PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>)

     # Boost.
    target_link_libraries(keplerian_toolbox PUBLIC Boost::boost Boost::serialization Boost::date_time)

    # Build Tests and link them to static library.
    if(PYKEP_BUILD_TESTS)
        add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
        file(COPY "${CMAKE_SOURCE_DIR}/tests/data/sgp4_test.txt" DESTINATION "${CMAKE_BINARY_DIR}/tests")
        file(COPY "${CMAKE_SOURCE_DIR}/tests/data/C_G_1000012_2012_2017.bsp" DESTINATION "${CMAKE_BINARY_DIR}/tests")
        file(COPY "${CMAKE_SOURCE_DIR}/tests/data/pck00010.tpc" DESTINATION "${CMAKE_BINARY_DIR}/tests")
        FILE(COPY "${CMAKE_SOURCE_DIR}/tests/data/gm_de431.tpc" DESTINATION "${CMAKE_BINARY_DIR}/tests")
    endif(PYKEP_BUILD_TESTS)

    # Configure config.hpp.
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/include/keplerian_toolbox/config.hpp" @ONLY)
    # Configure the doc files.
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py" @ONLY)
endif()

# Build pykep and link it to dynamic library.
if(PYKEP_BUILD_PYKEP)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pykep/__init__.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/pykep/__init__.py" @ONLY)
	add_subdirectory("${CMAKE_SOURCE_DIR}/pykep")
    if(MINGW OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        message(STATUS "Creating the files for the generation of a binary wheel.")
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py" @ONLY)
        if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
            # NOTE: this is necessary on linux but harmful on mingw.
            configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.cfg" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.cfg" @ONLY)
        endif()
        if(MINGW)
            configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" @ONLY)
        endif()
    endif()
endif()


# Library installation.
if(PYKEP_BUILD_KEP_TOOLBOX)
    # Installation of the header files.
    install(DIRECTORY include/ DESTINATION include)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/include/keplerian_toolbox/config.hpp" DESTINATION include/keplerian_toolbox)
    # Installation of the library.
    install(TARGETS keplerian_toolbox
        EXPORT keplerian_toolbox_export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
    )
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/keplerian_toolbox-config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/keplerian_toolbox-config.cmake" @ONLY)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/keplerian_toolbox-config.cmake" DESTINATION "lib/cmake/keplerian_toolbox")
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/PykepFindBoost.cmake" DESTINATION "lib/cmake/keplerian_toolbox")
    install(EXPORT keplerian_toolbox_export NAMESPACE Keplerian_toolbox:: DESTINATION lib/cmake/keplerian_toolbox)
    # Take care of versioning.
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/keplerian_toolbox-config-version.cmake" VERSION ${keplerian_toolbox_VERSION}
        COMPATIBILITY ExactVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/keplerian_toolbox-config-version.cmake" DESTINATION "lib/cmake/keplerian_toolbox")
endif()

# Uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)
    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()







