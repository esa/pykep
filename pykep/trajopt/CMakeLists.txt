set(PYKEP_TRAJOPT_PYTHON_FILES 
    __init__.py 
    _sf_point2point.py 
    _sf_pl2pl.py 
    _sf_pl2pl_alpha.py 
    _zoh_point2point.py
    _zoh_ss_point2point.py
    _zoh_pl2pl.py
    _mga.py _mga_1dsm.py
    _min_Bu_bu.py 
    _launchers.py 
    _mim.py
    _pl2pl_N_impulses.py
    _pontryagin_cartesian.py
    _pontryagin_equinoctial.py
    _primer_vector.py)
    
install(FILES ${PYKEP_TRAJOPT_PYTHON_FILES} DESTINATION ${_PYKEP_INSTALL_DIR}/trajopt)

ADD_SUBDIRECTORY(gym)