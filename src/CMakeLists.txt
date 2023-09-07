#
# VoroCrust src CMake
#
vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

# Force CMake to re-run the vc_git_tools.cmake script every
# time we build VoroCrust to ensure the SHA1 being set is accurate.
# Note: This won't work if this is not a git clone -- this will be replaced
#       when VoroCrust starts having numbered releases.
if(VOROCRUST_VERSION_USES_GIT_SHA1)
    add_custom_target(update_git_sha1 COMMAND cmake -P "${CMAKE_SOURCE_DIR}/cmake/vc_git_tools.cmake")
else()
    vc_message("Skip version configuration via Git SHA1")
endif()

# Add Subdirectory
add_subdirectory(Meshing) 

vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")

