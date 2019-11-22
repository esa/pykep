# Setup of the pykep sims_flanagan module.
YACMA_PYTHON_MODULE(util 
util.cpp 
)
target_link_libraries(util PRIVATE ${PYKEP_BP_TARGET} pykep)
target_compile_options(util PRIVATE "$<$<CONFIG:DEBUG>:${KEP_TOOLBOX_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${KEP_TOOLBOX_CXX_FLAGS_RELEASE}>")
set_property(TARGET core PROPERTY CXX_STANDARD 11)
set_property(TARGET core PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET core PROPERTY CXX_EXTENSIONS NO)

# Setup the installation path.
set(PYKEP_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/pykep")
install(TARGETS util
RUNTIME DESTINATION ${PYKEP_INSTALL_PATH}/util
LIBRARY DESTINATION ${PYKEP_INSTALL_PATH}/util
)

install(FILES __init__.py DESTINATION ${PYKEP_INSTALL_PATH}/util)
install(FILES gravity_spherical_harmonic.py DESTINATION ${PYKEP_INSTALL_PATH}/util)
install(FILES load_gravity_model.py DESTINATION ${PYKEP_INSTALL_PATH}/util)
install(DIRECTORY gravity_models DESTINATION ${PYKEP_INSTALL_PATH}/util)
