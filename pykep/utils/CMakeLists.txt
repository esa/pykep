set(PYKEP_UTILS_PYTHON_FILES 
__init__.py 
_spice_utils.py 
_planet_to_keplerian.py 
_encoding_conversions.py
_knn.py)

install(FILES ${PYKEP_UTILS_PYTHON_FILES} DESTINATION ${_PYKEP_INSTALL_DIR}/utils)
