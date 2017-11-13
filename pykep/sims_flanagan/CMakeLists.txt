ADD_LIBRARY(_sims_flanagan SHARED sims_flanagan.cpp)
PYTHON_CPPMODULE(_sims_flanagan)
INSTALL(FILES __init__.py DESTINATION ${PYKEP_INSTALL_PATH}/sims_flanagan)
