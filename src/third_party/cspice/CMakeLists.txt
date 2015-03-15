# Add all c files to the library
FILE(GLOB CSPICE_LIB_SRC_LIST "*.c")

#Build Static Library
ADD_LIBRARY(cspice STATIC ${CSPICE_LIB_SRC_LIST})

# Disable all warnings for SPICE compilation
TARGET_COMPILE_OPTIONS(cspice PUBLIC "-w")



