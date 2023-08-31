#
# External Project for Exodus support
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(VOROCRUST_ENABLE_EXODUS)

if(VOROCRUST_ENABLE_EXODUS)

    if(WIN32)
        if(NOT VOROCRUST_ENABLE_DEVTEST)
            # Guard against enabling Exodus on Windows
            vc_message("VOROCRUST_ENABLE_EXODUS does not work on Windows (yet, but we're working on it!).")
            vc_message("- Set VOROCRUST_ENABLE_DEVTEST:BOOL=ON to disable this error (developers only!)")
            message(FATAL_ERROR "${TagPREFIX} ERROR: Building VoroCrust with EXODUS enabled on Windows is not supported.")
        endif()
    endif()

    # Set the flag that we have Exodus available
    set(VOROCRUST_USE_EXODUS_VAR "USE_EXODUS")

    else()

    set(VOROCRUST_USE_EXODUS_VAR "NO_EXODUS")

endif(VOROCRUST_ENABLE_EXODUS)

vc_message_var(VOROCRUST_USE_EXODUS_VAR)

vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
