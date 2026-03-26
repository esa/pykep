set(PYKEP_PLOT_PYTHON_FILES 
    __init__.py 
    _planet.py 
    _lambert.py 
    _ballistic.py 
    _sf_leg.py 
    _mit.py)

install(FILES ${PYKEP_PLOT_PYTHON_FILES} DESTINATION ${_PYKEP_INSTALL_DIR}/plot)
