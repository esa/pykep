# Configure the version file.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/_version.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/_version.py" @ONLY)   

# The list of pykep's Python files.
set(PYKEP_PYTHON_FILES 
    __init__.py 
    _patch_planet.py 
    _patch_epoch.py 
    _version.py 
    test.py 
    test_ta.py 
    test_trajopt.py 
    test_utils.py 
    test_leg_sims_flanagan.py
    test_leg_zoh.py)

# Core module.
Python3_add_library(core MODULE WITH_SOABI 
    core.cpp
    docstrings.cpp
    common_utils.cpp
    python_udpla.cpp
    expose_udplas.cpp
)

# This is to silence an aggressive compiler setting creating false positives in pybind11 stuff.
if (YACMA_COMPILER_IS_GNUCXX)
    set_source_files_properties(core.cpp PROPERTIES COMPILE_FLAGS "-Wno-null-dereference")
endif()

target_link_libraries(core PRIVATE kep3)
target_link_libraries(core PRIVATE "${pybind11_LIBRARIES}")

target_include_directories(core SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}" "${Python3_INCLUDE_DIRS}")
target_compile_definitions(core PRIVATE "${pybind11_DEFINITIONS}")
target_compile_options(core PRIVATE
    "$<$<CONFIG:Debug>:${kep3_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${kep3_CXX_FLAGS_RELEASE}>"
    
)
target_include_directories(core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
set_target_properties(core PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_target_properties(core PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)
target_compile_features(core PRIVATE cxx_std_17)
set_property(TARGET core PROPERTY CXX_EXTENSIONS NO)

# Installation setup.
if(PYKEP_INSTALL_PATH STREQUAL "")
    message(VERBOSE "pykep install location (default): ${Python3_SITEARCH}")
    set(_PYKEP_INSTALL_DIR "${Python3_SITEARCH}/pykep")
else()
    message(VERBOSE "pykep install location (custom): ${PYKEP_INSTALL_PATH}")
    set(_PYKEP_INSTALL_DIR "${PYKEP_INSTALL_PATH}/pykep")
endif()

# Install the core module.
install(TARGETS core
    RUNTIME DESTINATION ${_PYKEP_INSTALL_DIR}
    LIBRARY DESTINATION ${_PYKEP_INSTALL_DIR}
)

# Install the Python files.
install(FILES ${PYKEP_PYTHON_FILES} DESTINATION ${_PYKEP_INSTALL_DIR})
ADD_SUBDIRECTORY(udpla)
ADD_SUBDIRECTORY(leg)
ADD_SUBDIRECTORY(ta)
ADD_SUBDIRECTORY(trajopt)
ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(data)
ADD_SUBDIRECTORY(plot)


unset(_PYKEP_INSTALL_DIR)