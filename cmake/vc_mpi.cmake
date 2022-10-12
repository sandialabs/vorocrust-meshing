#
# MPI Configuration
#
#
include_guard()
vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(VOROCRUST_ENABLE_MPI)

set(VOROCRUST_HAVE_MPICXX_VAR "NO_MPI")

if(VOROCRUST_ENABLE_MPI)

    find_package(MPI REQUIRED)

    # Sets:
    #       (3.0 compatibility)
    vc_message_var(MPI_FOUND)
    vc_message_var(MPI_CXX_FOUND)
    vc_message_var(MPI_CXX_VERSION)
    vc_message_var(MPI_CXX_COMPILER)
    vc_message_var(MPI_CXX_COMPILER_FLAGS)
    vc_message_var(MPI_CXX_INCLUDE_PATH)
    vc_message_var(MPI_CXX_LIBRARIES)
    vc_message_var(MPI_CXX_LINK_FLAGS)
    vc_message_var(MPIEXEC_EXECUTABLE)
    vc_message_var(MPIEXEC_NUMPROC_FLAGS)
    vc_message_var(MPIEXEC_MAX_NUMPROCS)

    if(MPI_CXX_FOUND)
        set(VOROCRUST_HAVE_MPICXX_VAR "USE_MPI")
        include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
    else()
        # MPI was requested, but we didn't find it.
        message(WARNING "${TagPREFIX} MPI was requested but was not found.")
    endif()

endif()

vc_message_var(VOROCRUST_HAVE_MPICXX_VAR)

# ------------------------------------------------
# Check for configuration failure conditions
# ------------------------------------------------
if(VOROCRUST_ENABLE_MPI)
    if(NOT MPI_CXX_FOUND)
        # User wants MPI enabled, but cmake can't find an MPI C++ compiler.
        message(FATAL_ERROR "Error: No MPI C++ compiler found, but configuration requested it.")
    endif()
endif()


# MPI settings and dependencies
macro(VOROCRUST_TPLADD_MPI TARGET)
    if(VOROCRUST_ENABLE_MPI)
        target_include_directories(${TARGET} PUBLIC ${MPI_CXX_INCLUDE_PATH})
        target_compile_options(${TARGET}     PUBLIC ${MPI_CXX_COMPILE_FLAGS})
        target_link_libraries(${TARGET}      PUBLIC ${MPI_CXX_LIBRARIES})
    endif()
endmacro()


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
