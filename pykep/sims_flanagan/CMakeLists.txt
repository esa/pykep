# Setup of the pykep sims_flanagan module.
YACMA_PYTHON_MODULE(sims_flanagan 
    sims_flanagan.cpp 
)
target_link_libraries(sims_flanagan PRIVATE ${PYKEP_BP_TARGET} pykep)
target_compile_options(sims_flanagan PRIVATE "$<$<CONFIG:DEBUG>:${KEP_TOOLBOX_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${KEP_TOOLBOX_CXX_FLAGS_RELEASE}>")
set_property(TARGET sims_flanagan PROPERTY CXX_STANDARD 14)
set_property(TARGET sims_flanagan PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET sims_flanagan PROPERTY CXX_EXTENSIONS NO)

# Setup the installation path.
set(PYKEP_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/pykep")
install(TARGETS sims_flanagan
RUNTIME DESTINATION ${PYKEP_INSTALL_PATH}/sims_flanagan
LIBRARY DESTINATION ${PYKEP_INSTALL_PATH}/sims_flanagan
)

install(FILES __init__.py DESTINATION ${PYKEP_INSTALL_PATH}/sims_flanagan)
