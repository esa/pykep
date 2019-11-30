# Setup of the dcgpy core module.
YACMA_PYTHON_MODULE(core 
    core.cpp 
    docstrings.cpp
)
target_link_libraries(core PRIVATE ${PYKEP_BP_TARGET} pykep)
target_compile_options(core PRIVATE "$<$<CONFIG:DEBUG>:${KEP_TOOLBOX_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${KEP_TOOLBOX_CXX_FLAGS_RELEASE}>")
set_property(TARGET core PROPERTY CXX_STANDARD 11)
set_property(TARGET core PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET core PROPERTY CXX_EXTENSIONS NO)

# Setup the installation path.
set(PYKEP_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/pykep")
install(TARGETS core
RUNTIME DESTINATION ${PYKEP_INSTALL_PATH}/core
LIBRARY DESTINATION ${PYKEP_INSTALL_PATH}/core
)

# Add the Python files.
install(FILES __init__.py DESTINATION ${PYKEP_INSTALL_PATH}/core)
