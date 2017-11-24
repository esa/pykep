# Setup of the pygmo core module.
YACMA_PYTHON_MODULE(_planet planet.cpp)
target_link_libraries(_planet PRIVATE ${PYKEP_BP_TARGET} ${PYKEP_BOOST_TARGET})
target_compile_options(_planet PRIVATE "$<$<CONFIG:DEBUG>:${PAGMO_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${PAGMO_CXX_FLAGS_RELEASE}>")
set_property(TARGET _planet PROPERTY CXX_STANDARD 11)
set_property(TARGET _planet PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET _planet PROPERTY CXX_EXTENSIONS NO)

install(TARGETS _planet
RUNTIME DESTINATION ${PYKEP_INSTALL_PATH}/planet
LIBRARY DESTINATION ${PYKEP_INSTALL_PATH}/planet
)

install(FILES __init__.py DESTINATION ${PYKEP_INSTALL_PATH}/planet)
