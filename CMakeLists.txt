# NOTE: current minimum kep3 requirement.
cmake_minimum_required(VERSION 3.28.0)
# This effectively requests behavior preferred as of a given CMake version and tells newer CMake versions to warn about their new policies.
cmake_policy(VERSION 3.30)

# This removes a strange behaviour on osx CI machines where ninja would error out with 
# '/bin/sh: CMAKE_CXX_COMPILER_CLANG_SCAN_DEPS-NOTFOUND: command not found' even if clang-scan-deps was installed
# Since we here do not use modules we deactivate their search.
set(CMAKE_CXX_SCAN_FOR_MODULES 0)

# Set default build type to "Release".
# NOTE: this should be done before the project command since the latter can set
# CMAKE_BUILD_TYPE itself (it does so for nmake).
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
endif()

project(kep3 VERSION 3.0.0 LANGUAGES CXX C)

option(KEP3_VERBOSE_CONFIGURE "Enable detailed configure-time logging." OFF)
if(KEP3_VERBOSE_CONFIGURE)
    # Propagate verbose messages without requiring extra CLI flags.
    set(CMAKE_MESSAGE_LOG_LEVEL VERBOSE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/yacma")

# Run the YACMA compiler setup.
include(YACMACompilerLinkerSettings)

# Build options.
option(kep3_BUILD_TESTS "Build unit tests." OFF)
option(kep3_BUILD_BENCHMARKS "Build benchmarks." OFF)
option(kep3_BUILD_PYTHON_BINDINGS "Build Python bindings." OFF)

# NOTE: on Unix systems, the correct library installation path
# could be something other than just "lib", such as "lib64",
# "lib32", etc., depending on platform/configuration. Apparently,
# CMake provides this information via the GNUInstallDirs module.
# Let's enable this for now on all Unixes except OSX.
# NOTE: potentially, this could be applicable to Cygwin as well.
#
# https://cmake.org/cmake/help/v3.15/module/GNUInstallDirs.html
# https://cmake.org/pipermail/cmake/2013-July/055375.html
if(UNIX AND NOT APPLE)
    include(GNUInstallDirs)
    set(_kep3_INSTALL_LIBDIR_DEFAULT "${CMAKE_INSTALL_LIBDIR}")
else()
    set(_kep3_INSTALL_LIBDIR_DEFAULT "lib")
endif()
if(NOT kep3_INSTALL_LIBDIR)
    set(kep3_INSTALL_LIBDIR "${_kep3_INSTALL_LIBDIR_DEFAULT}" CACHE STRING
        "Library installation directory." FORCE)
endif()
mark_as_advanced(kep3_INSTALL_LIBDIR)

# Assemble the flags.
set(kep3_CXX_FLAGS_DEBUG ${YACMA_CXX_FLAGS} ${YACMA_CXX_FLAGS_DEBUG})
set(kep3_CXX_FLAGS_RELEASE ${YACMA_CXX_FLAGS})
if(YACMA_COMPILER_IS_MSVC)
  # On both cl and clang-cl, disable the idiotic minmax macros and enable the bigobj option.
  # Also, enable the WIN32_LEAN_AND_MEAN definition:
  # https://stackoverflow.com/questions/11040133/what-does-defining-win32-lean-and-mean-exclude-exactly
  list(APPEND kep3_CXX_FLAGS_DEBUG "-DNOMINMAX" "/bigobj" "-DWIN32_LEAN_AND_MEAN")
  list(APPEND kep3_CXX_FLAGS_RELEASE "-DNOMINMAX" "/bigobj" "-DWIN32_LEAN_AND_MEAN")
  if(YACMA_COMPILER_IS_CLANGXX)
    # clang-cl emits various warnings, let's just silence them.
    # NOTE: at one point in the recent past, MSVC added an options similar to GCC's isystem:
    # https://blogs.msdn.microsoft.com/vcblog/2017/12/13/broken-warnings-theory/
    # We probably just need to wait for this to be picked up by CMake/clang-cl. Let's
    # revisit the issue in the future.
    list(APPEND _kep3_CLANG_CL_DISABLED_WARNINGS
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
    list(APPEND kep3_CXX_FLAGS_DEBUG ${_kep3_CLANG_CL_DISABLED_WARNINGS})
    list(APPEND kep3_CXX_FLAGS_RELEASE ${_kep3_CLANG_CL_DISABLED_WARNINGS})
    unset(_kep3_CLANG_CL_DISABLED_WARNINGS)
  else()
    # Same as above, disable some cl warnings.
    list(APPEND kep3_CXX_FLAGS_DEBUG "/wd4459" "/wd4127" "/wd4251")
    list(APPEND kep3_CXX_FLAGS_RELEASE "/wd4459" "/wd4127" "/wd4251")
  endif()
  # Enable strict conformance mode, if supported.
  set(CMAKE_REQUIRED_QUIET TRUE)
  check_cxx_compiler_flag("/permissive-" _kep3_MSVC_SUPPORTS_STRICT_CONFORMANCE)
  unset(CMAKE_REQUIRED_QUIET)
  if(_kep3_MSVC_SUPPORTS_STRICT_CONFORMANCE)
        message(VERBOSE "MSVC: enabling '/permissive-' strict conformance mode.")
    list(APPEND kep3_CXX_FLAGS_DEBUG "/permissive-")
    list(APPEND kep3_CXX_FLAGS_RELEASE "/permissive-")
  endif()
  unset(_kep3_MSVC_SUPPORTS_STRICT_CONFORMANCE)
endif()

# List of source files.
set(kep3_SRC_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epoch.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/planet.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lambert_problem.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/linalg.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/udpla/keplerian.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/udpla/jpl_lp.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/udpla/vsop2013.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/leg/sims_flanagan.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/leg/sims_flanagan_alpha.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/leg/zoh.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/leg/zoh_ss.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/leg/sf_checks.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core_astro/flyby.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core_astro/ic2par2ic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core_astro/ic2mee2ic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core_astro/mima.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core_astro/mee2par2mee.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core_astro/stm.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core_astro/propagate_lagrangian.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core_astro/encodings.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/core_astro/basic_transfers.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ta/kep.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ta/zoh_kep.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ta/zoh_eq.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ta/zoh_cr3bp.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ta/zoh_ss.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ta/cr3bp.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ta/bcp.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ta/pontryagin_cartesian.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ta/pontryagin_equinoctial.cpp"
)

# Setup of the kep3 shared library.
add_library(kep3 SHARED "${kep3_SRC_FILES}")
set_property(TARGET kep3 PROPERTY VERSION "1.0")
set_property(TARGET kep3 PROPERTY SOVERSION 1)
set_target_properties(kep3 PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_target_properties(kep3 PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)

# Adding this to shut down xtensor warnings (as of v0.26 still there)
if (YACMA_COMPILER_IS_GNUCXX)
    target_compile_options(kep3 PRIVATE -Wno-array-bounds)
    target_compile_options(kep3 PRIVATE -Wno-stringop-overflow)
endif()


target_compile_options(kep3 PRIVATE
    "$<$<CONFIG:Debug>:${kep3_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${kep3_CXX_FLAGS_RELEASE}>"
)

# Ensure that C++20 is employed when both compiling and consuming kep3.
target_compile_features(kep3 PUBLIC cxx_std_23)
# Enforce vanilla C++20 when compiling kep3.
set_property(TARGET kep3 PROPERTY CXX_EXTENSIONS NO)

target_include_directories(kep3 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
# Boost.
# NOTE: need 1.69 for safe numerics.
set(_kep3_MIN_BOOST_VERSION "1.69")
set(_kep3_MIN_FMT_VERSION "10")
set(_kep3_MIN_HEYOKA_VERSION "7")

find_package(Boost ${_kep3_MIN_BOOST_VERSION} CONFIG REQUIRED COMPONENTS serialization)
target_link_libraries(kep3 PUBLIC Boost::boost Boost::serialization)
# NOTE: quench warnings from Boost when building the library.
target_compile_definitions(kep3 PRIVATE BOOST_ALLOW_DEPRECATED_HEADERS)

# fmt.
find_package(fmt ${_kep3_MIN_FMT_VERSION} CONFIG REQUIRED)
target_link_libraries(kep3 PUBLIC fmt::fmt)

# heyoka.
find_package(heyoka ${_kep3_MIN_HEYOKA_VERSION} CONFIG REQUIRED)
target_link_libraries(kep3 PUBLIC heyoka::heyoka)

# spdlog.
find_package(spdlog CONFIG REQUIRED)
target_link_libraries(kep3 PRIVATE spdlog::spdlog)

# xtensor.
find_package(xtensor CONFIG REQUIRED)
target_link_libraries(kep3 PRIVATE xtensor)

# xtensor.
find_package(xtensor-blas CONFIG REQUIRED)
target_link_libraries(kep3 PRIVATE xtensor-blas)

# Configure config.hpp.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/include/kep3/config.hpp" @ONLY)

# Installation of the header files.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/kep3" DESTINATION include)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/kep3/config.hpp" DESTINATION include/kep3)

# Installation of the library.
install(TARGETS kep3
    EXPORT kep3_export
    LIBRARY DESTINATION "${kep3_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${kep3_INSTALL_LIBDIR}"
    RUNTIME DESTINATION bin
)

include(CMakePackageConfigHelpers)

# Setup of the CMake config file.
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/kep3-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/kep3-config.cmake"
    INSTALL_DESTINATION "${kep3_INSTALL_LIBDIR}/cmake/kep3"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kep3-config.cmake"
    DESTINATION "${kep3_INSTALL_LIBDIR}/cmake/kep3")
install(EXPORT kep3_export NAMESPACE kep3:: DESTINATION "${kep3_INSTALL_LIBDIR}/cmake/kep3")
# Take care of versioning.
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/kep3-config-version.cmake" COMPATIBILITY SameMinorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kep3-config-version.cmake" DESTINATION "${kep3_INSTALL_LIBDIR}/cmake/kep3")


if(kep3_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

if(kep3_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)

endif()

if(kep3_BUILD_PYTHON_BINDINGS)
    # Find python
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
    message(VERBOSE "Python3 interpreter: ${Python3_EXECUTABLE}")
    message(VERBOSE "Python3 installation directory: ${Python3_SITEARCH}")
    message(VERBOSE "Python3 include directories: ${Python3_INCLUDE_DIRS}")

    set(PYKEP_INSTALL_PATH "" CACHE STRING "pykep module installation path")
    mark_as_advanced(PYKEP_INSTALL_PATH)

    # pybind11.
    find_package(pybind11 REQUIRED)

    # Build directory
    add_subdirectory("${CMAKE_SOURCE_DIR}/pykep")
endif()

message(STATUS "Config summary: kep3 ${kep3_VERSION} on ${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR} with ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Config summary: build type=${CMAKE_BUILD_TYPE}, tests=${kep3_BUILD_TESTS}, benchmarks=${kep3_BUILD_BENCHMARKS}, python=${kep3_BUILD_PYTHON_BINDINGS}, verbose=${KEP3_VERBOSE_CONFIGURE}")
message(STATUS "Config summary: install libdir=${kep3_INSTALL_LIBDIR}, deps Boost>=${_kep3_MIN_BOOST_VERSION}, fmt>=${_kep3_MIN_FMT_VERSION}, heyoka>=${_kep3_MIN_HEYOKA_VERSION}")

