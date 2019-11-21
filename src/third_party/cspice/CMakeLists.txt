# Add all c files to the library
file(GLOB CSPICE_LIB_SRC_LIST "*.c")

#Build Static Library
add_library(cspice STATIC ${CSPICE_LIB_SRC_LIST})
set_property(TARGET cspice PROPERTY POSITION_INDEPENDENT_CODE ON)

# Disable all warnings for SPICE compilation
target_compile_options(cspice PUBLIC "-w")



