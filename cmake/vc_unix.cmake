#
# Unix like systems (Linux, OSX, etc.) specific build configurations go here.
#
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

include(vc_functions)

if(NOT (UNIX OR APPLE))
    vc_message("Unix/Darwin NOT detected")
    vc_message("COMPLETE: ${CMAKE_CURRENT_SOURCE_DIR}/vc_unix.cmake")
    vc_message(" ")
    return()
endif()


# Unix and Apple Systems (We assume gcc is being used)
vc_message("System Detected: *nix or Darwin")

# Configure paths for TPLS
set(TPL_LIBRARY_DIR ${TPL_INSTALL_ROOT_DIR}/lib)
set(TPL_INCLUDE_DIR ${TPL_INSTALL_ROOT_DIR}/include)


# Enable C++ flags if the compiler supports them and if
# VOROCRUST_ENABLE_WARNINGS is enabled.
# Add  -DVOROCRUST_ENABLE_WARNINGS:BOOL=ON  to CMake command line to enable
if(VOROCRUST_ENABLE_WARNINGS)
    enable_flag_if_supported_unique( CMAKE_CXX_FLAGS "-Wall" )
    enable_flag_if_supported_unique( CMAKE_CXX_FLAGS "-Wextra" )
endif()

enable_flag_if_supported_unique( CMAKE_CXX_FLAGS "-Wno-unknown-pragmas")
enable_flag_if_supported_unique( CMAKE_CXX_FLAGS "-fPIC")

# SEACAS + TPLs may need this
set(VC_LINKFLAG_LIBDL "-ldl")

# STATIC linking flags (problematic on OSX, not tested on *nix)
#set(BUILD_SHARED_LIBS OFF CACHE STRING "Set by CMake")
#enable_cxxflag_if_supported(CMAKE_EXE_LINKER_FLAGS "-static")
#enable_cxxflag_if_supported(CMAKE_STATIC_LINKER_FLAGS "-static")


# If VOROCRUST_ENABLE_COVERAGE_TESTING is enabled, gcc needs the -gc parameter
if(VOROCRUST_ENABLE_COVERAGE_TESTING)
endif()

vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
