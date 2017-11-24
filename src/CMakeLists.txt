# Keplerian_toolbox lib source files.
SET(PYKEP_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/epoch.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/lambert_problem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sims_flanagan/leg.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sims_flanagan/leg_s.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sims_flanagan/spacecraft.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/core_functions/jorba.c

	${CMAKE_CURRENT_SOURCE_DIR}/planet/base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet/keplerian.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet/j2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet/jpl_low_precision.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet/mpcorb.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet/tle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet/gtoc2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet/gtoc5.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet/gtoc6.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet/gtoc7.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/Util.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/Tle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/SolarPosition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/SGP4.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/OrbitalElements.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/Observer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsgp4/Eci.cpp
)

IF(BUILD_SPICE)
	SET(PYKEP_LIB_SRC_LIST ${PYKEP_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/planet/spice.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/util/spice_utils.cpp
	)
ENDIF()

if(BUILD_MAIN OR BUILD_TESTS)
	ADD_LIBRARY(keplerian_toolbox_static STATIC ${PYKEP_LIB_SRC_LIST})
    set_property(TARGET keplerian_toolbox_static PROPERTY CXX_STANDARD 11)
    set_property(TARGET keplerian_toolbox_static PROPERTY CXX_STANDARD_REQUIRED YES)
	set_property(TARGET keplerian_toolbox_static PROPERTY CXX_EXTENSIONS NO)
	target_link_libraries(keplerian_toolbox "Boost::boost")	
endif()

#Build Dynamic Library
if(BUILD_PYKEP)
	set(LIB_INSTALL_PATH "lib")
	set(KEP_TOOL_NEEDED "Boost::boost" "Boost::date_time" "Boost::serialization")
	if(BUILD_SPICE)
		list(APPEND KEP_TOOL_NEEDED "cspice")
	endif()
	add_library(keplerian_toolbox SHARED ${PYKEP_LIB_SRC_LIST})
	target_link_libraries(keplerian_toolbox ${KEP_TOOL_NEEDED})
	set_target_properties(keplerian_toolbox PROPERTIES COMPILE_FLAGS "-DKEP_TOOL_DLL_EXPORT_API")
    set_property(TARGET keplerian_toolbox PROPERTY CXX_STANDARD 11)
    set_property(TARGET keplerian_toolbox PROPERTY CXX_STANDARD_REQUIRED YES)
    set_property(TARGET keplerian_toolbox PROPERTY CXX_EXTENSIONS NO)
	install(TARGETS keplerian_toolbox
		RUNTIME DESTINATION ${LIB_INSTALL_PATH}
		LIBRARY DESTINATION ${LIB_INSTALL_PATH}
	)
endif()

#Installs the headers
IF(INSTALL_HEADERS)
	FILE(GLOB HEADERS_IN_SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
	INSTALL(FILES ${HEADERS_IN_SRC} DESTINATION "include/keplerian_toolbox")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/core_functions" DESTINATION "include/keplerian_toolbox" FILES_MATCHING PATTERN "*.h")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/numerics" DESTINATION "include/keplerian_toolbox" FILES_MATCHING PATTERN "*.h")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/planet" DESTINATION "include/keplerian_toolbox" FILES_MATCHING PATTERN "*.h")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/sims_flanagan" DESTINATION "include/keplerian_toolbox" FILES_MATCHING PATTERN "*.h")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/third_party" DESTINATION "include/keplerian_toolbox" FILES_MATCHING PATTERN "*.h")
	INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/util" DESTINATION "include/keplerian_toolbox" FILES_MATCHING PATTERN "*.h")
ENDIF(INSTALL_HEADERS)
