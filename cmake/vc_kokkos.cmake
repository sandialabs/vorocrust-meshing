# Configuration file for TPL: Kokkos
#
# URL: https://github.com/kokkos/kokkos/archive/3.1.01.zip
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

include(FetchContent)
vc_message_var(VOROCRUST_TPL_ENABLE_KOKKOS)

if(VOROCRUST_TPL_ENABLE_KOKKOS)

    # This is set in the top-level CMakeLists.txt file and shows
    # up as an _advanced_ option.
    vc_message_var(VOROCRUST_TPL_KOKKOS_VERSION)

    # Kokkos requires a CXX standard of at least c++14
    set(Kokkos_CXX_STANDARD "14" CACHE STRING "Set by VoroCrust")

    FetchContent_Declare(
        Kokkos
        GIT_REPOSITORY  "https://github.com/kokkos/kokkos.git"
        GIT_TAG         ${VOROCRUST_TPL_KOKKOS_VERSION}
        #GIT_SHALLOW     ON
        GIT_PROGRESS    ON
    )
    # Note: Git shallow only works if we have a real branch to check out.
    #       It won't work if things are in a headless state, etc.

    FetchContent_GetProperties(Kokkos)
    if(NOT kokkos_POPULATED)
        vc_message("BEGIN   : Download Kokkos")
        FetchContent_Populate( Kokkos )
        vc_message("COMPLETE: Download Kokkos")
    endif()

    if(WIN32)
        vc_message("Disabling some Kokkos parameters for Windows compatibility")

        # There is an error when configuring on Windows with LIBDL. Disabling
        # this in Kokkos lets us configure and build but we won't be able to
        # use profiling (for now).
        set(Kokkos_ENABLE_PROFILING OFF
            CACHE BOOL
            "Disabled by VoroCrust due to Windows compatibility issues" FORCE)
        set(Kokkos_ENABLE_LIBDL OFF
            CACHE BOOL
            "Disabled by VoroCrust due to Windows compatibility issues" FORCE)
    endif(WIN32)


    if(NOT WIN32)
        # Visual Studio does not currently fully support OpenMP and there
        # are features in Kokkos that aren't available to Visual Studio
        # compilers yet.
        if(VOROCRUST_ENABLE_OPENMP)
            set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Enabled by VOROCRUST_ENABLE_OPENMP" FORCE)
        endif(VOROCRUST_ENABLE_OPENMP)
    endif(NOT WIN32)

    # Add in the Kokkos subdirectory.
    add_subdirectory(${kokkos_SOURCE_DIR} ${kokkos_BINARY_DIR})

    # Set the flag that we have Kokkos available
    set(VOROCRUST_HAVE_KOKKOS_VAR "USE_KOKKOS")

    vc_message_var(VOROCRUST_TPL_KOKKOS_ROOT_DIR)
    vc_message_var(CMAKE_BINARY_DIR)
    vc_message_var(CMAKE_CURRENT_BINARY_DIR)
    vc_message_var(kokkos_SOURCE_DIR)
    vc_message_var(kokkos_BINARY_DIR)

else(VOROCRUST_TPL_ENABLE_KOKKOS)

    set(VOROCRUST_HAVE_KOKKOS_VAR "NO_KOKKOS")

endif(VOROCRUST_TPL_ENABLE_KOKKOS)

vc_message_var(VOROCRUST_HAVE_KOKKOS_VAR)


# Kokkos settings and dependencies
macro(VOROCRUST_TPLADD_KOKKOS TARGET)
    if(VOROCRUST_TPL_ENABLE_KOKKOS)
        target_link_libraries(${TARGET} PUBLIC kokkos)
        target_link_libraries(${TARGET} PUBLIC Kokkos::kokkos)
    endif(VOROCRUST_TPL_ENABLE_KOKKOS)
endmacro()


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
