#
# CMake system file for Windows systems
#
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

if(NOT WIN32)
    vc_message("Windows (WIN32) was NOT detected")
    vc_message("COMPLETE: ${CMAKE_CURRENT_SOURCE_DIR}/vc_windows.cmake")
    vc_message(" ")
    return()
endif()

# Windows Specific flags (MSVC)
vc_message("System Detected: WIN32")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/usr/local/include"
    CACHE STRING "omp include")

set(CMAKE_EXE_LINKER_FLAGS "-lomp -L/usr/local/lib/omp" CACHE STRING "omp lib")

# SEACAS + TPLs need -ldl (TODO: what is the windows equivalent?)
set(VC_LINKFLAG_LIBDL "")


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
