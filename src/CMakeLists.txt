# Keplerian_toolbox lib source files.
SET(PYKEP_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/epoch.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/lambert_problem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sims_flanagan/leg.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/sims_flanagan/leg_s.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sims_flanagan/spacecraft.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/core_functions/jorba.c

	${CMAKE_CURRENT_SOURCE_DIR}/planets/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/keplerian.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/jpl_low_precision.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/mpcorb.cpp


	${CMAKE_CURRENT_SOURCE_DIR}/planets/planet.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/asteroid_gtoc2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/asteroid_gtoc5.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/asteroid_gtoc7.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/planet_mpcorb.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/planet_ss.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/planet_js.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planets/planet_tle.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/Vector.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/Util.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/Tle.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/TimeSpan.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/SolarPosition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/SGP4.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/OrbitalElements.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/Observer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/Eci.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/DateTime.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/CoordTopocentric.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/CoordGeodetic.cpp
)

#Build Static Library
ADD_LIBRARY(keplerian_toolbox_static STATIC ${PYKEP_LIB_SRC_LIST})

#Build Dynamic Library (only if needed by PyKEP)
IF(BUILD_PYKEP)
	SET(LIB_INSTALL_PATH "lib")
	#INCLUDE(FindPythonLibs)
	ADD_LIBRARY(keplerian_toolbox SHARED ${PYKEP_LIB_SRC_LIST})
	TARGET_LINK_LIBRARIES(keplerian_toolbox ${LIBRARIES_4_DYNAMIC})
	SET_TARGET_PROPERTIES(keplerian_toolbox PROPERTIES COMPILE_FLAGS "-DKEP_TOOL_DLL_EXPORT_API")
	INSTALL(TARGETS keplerian_toolbox
		RUNTIME DESTINATION ${LIB_INSTALL_PATH}
		LIBRARY DESTINATION ${LIB_INSTALL_PATH}
	)
ENDIF(BUILD_PYKEP)
