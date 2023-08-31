# Configuration file for TPL: Kokkos-Kernels
#
# URL: https://github.com/Kokkos/kokkos-kernels
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

include(FetchContent)
vc_message_var(VOROCRUST_TPL_ENABLE_KOKKOSKERNELS)


if(VOROCRUST_TPL_ENABLE_KOKKOSKERNELS)

    # This is set in the top-level CMakeLists.txt file and shows
    # up as an _advanced_ option.
    set(VOROCRUST_TPL_KOKKOSKERNELS_VERSION ${VOROCRUST_TPL_KOKKOS_VERSION})

    vc_message_var(VOROCRUST_TPL_KOKKOSKERNELS_VERSION)

    # Kokkos Kernels has some build issues on Windows platforms that need to be resolved.
    if(WIN32 AND (NOT VOROCRUST_ENABLE_DEVTEST))
        message(STATUS "${Red}ERROR ${ColorReset}Currently Kokkos-Kernels + Windows is not supported")
        message(FATAL_ERROR "To override this, enable ${Magenta}VOROCRUST_ENABLE_DEVTEST${ColorReset}")
    endif()

    FetchContent_Declare(
        KokkosKernels
        GIT_REPOSITORY  "https://github.com/kokkos/kokkos-kernels.git"
        GIT_TAG         ${VOROCRUST_TPL_KOKKOSKERNELS_VERSION}
        #GIT_PROGRESS    ON
        #GIT_SHALLOW     ON
    )
    # Note: Git shallow only works if we have a real branch to check out.
    #       It won't work if things are in a headless state, etc.

    FetchContent_GetProperties( KokkosKernels )

    if(NOT kokkoskernels_POPULATED)
        vc_message("BEGIN   : Download Kokkos-Kernels")
        FetchContent_Populate( KokkosKernels )
        vc_message("COMPLETE: Download Kokkos-Kernels")
    endif()

    # placeholder for now -- these settings are handled in Kokkos (I think)
    #if(WIN32)
    #    vc_message("Disabling some Kokkos-Kernels parameters for Windows compatibility")
    #endif()

    # Add in the Kokkos-Kernels subdirectory.
    vc_message_var(kokkoskernels_SOURCE_DIR)
    vc_message_var(kokkoskernels_BINARY_DIR)
    add_subdirectory(${kokkoskernels_SOURCE_DIR} ${kokkoskernels_BINARY_DIR})

    # Set the flag that we have Kokkos-Kernels available
    set(VOROCRUST_HAVE_KOKKOSKERNELS_VAR "USE_KOKKOSKERNELS")

    vc_message_var(VOROCRUST_TPL_KOKKOSKERNELS_ROOT_DIR)
    vc_message_var(CMAKE_BINARY_DIR)
    vc_message_var(CMAKE_CURRENT_BINARY_DIR)
    vc_message_var(kokkos_SOURCE_DIR)
    vc_message_var(kokkos_BINARY_DIR)

else(VOROCRUST_TPL_ENABLE_KOKKOS)

    set(VOROCRUST_HAVE_KOKKOSKERNELS_VAR "NO_KOKKOSKERNELS")

endif(VOROCRUST_TPL_ENABLE_KOKKOSKERNELS)

vc_message_var(VOROCRUST_HAVE_KOKKOSKERNELS_VAR)


# Kokkos-Kernels settings and dependencies
macro(VOROCRUST_TPLADD_KOKKOSKERNELS TARGET)
    if(VOROCRUST_TPL_ENABLE_KOKKOSKERNELS)
        target_link_libraries(${TARGET} PUBLIC kokkoskernels)
        target_link_libraries(${TARGET} PUBLIC Kokkos::kokkoskernels)
    endif()
endmacro()


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
