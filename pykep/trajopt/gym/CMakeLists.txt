set(PYKEP_GYM_PYTHON_FILES 
    __init__.py 
    _cassini_mga.py 
    _cassini_mga_1dsm.py 
    _em_Nimp.py
    _eve_mga_1dsm.py 
    _juice_mga_1dsm.py 
    _messenger.py
    _rosetta_mga_1dsm.py
    )
    
install(FILES ${PYKEP_GYM_PYTHON_FILES} DESTINATION ${_PYKEP_INSTALL_DIR}/trajopt/gym)
