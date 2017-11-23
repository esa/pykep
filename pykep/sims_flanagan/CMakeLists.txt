# Setup of the pygmo core module.
YACMA_PYTHON_MODULE(_sims_flanagan sims_flanagan.cpp)
target_link_libraries(_sims_flanagan PRIVATE ${PYKEP_BP_TARGET} ${PYKEP_BOOST_TARGET})
target_compile_options(_sims_flanagan PRIVATE "$<$<CONFIG:DEBUG>:${PAGMO_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${PAGMO_CXX_FLAGS_RELEASE}>")
set_property(TARGET _sims_flanagan PROPERTY CXX_STANDARD 11)
set_property(TARGET _sims_flanagan PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET _sims_flanagan PROPERTY CXX_EXTENSIONS NO)

install(TARGETS _sims_flanagan
RUNTIME DESTINATION ${PYKEP_INSTALL_PATH}/sims_flanagan
LIBRARY DESTINATION ${PYKEP_INSTALL_PATH}/sims_flanagan
)

install(FILES __init__.py DESTINATION ${PYKEP_INSTALL_PATH}/sims_flanagan)
