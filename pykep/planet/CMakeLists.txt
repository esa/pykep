# Setup of the pykep planet module.
YACMA_PYTHON_MODULE(planet 
    planet.cpp 
)
target_link_libraries(planet PRIVATE ${PYKEP_BP_TARGET} pykep)
target_compile_options(planet PRIVATE "$<$<CONFIG:DEBUG>:${KEP_TOOLBOX_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${KEP_TOOLBOX_CXX_FLAGS_RELEASE}>")
set_property(TARGET planet PROPERTY CXX_STANDARD 11)
set_property(TARGET planet PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET planet PROPERTY CXX_EXTENSIONS NO)

# Setup the installation path.
set(PYKEP_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/pykep")
install(TARGETS planet
RUNTIME DESTINATION ${PYKEP_INSTALL_PATH}/planet
LIBRARY DESTINATION ${PYKEP_INSTALL_PATH}/planet
)

install(FILES __init__.py DESTINATION ${PYKEP_INSTALL_PATH}/planet)
