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

# Add Subdirectories
# add_subdirectory(Common)
add_subdirectory(Meshing) # JHS 08/16/2022 Since everything is in Meshing, swap Common for Meshing
#add_subdirectory(Common)
#add_optional_subdirectory_if_exists("${CMAKE_CURRENT_SOURCE_DIR}/Meshing")
#add_optional_subdirectory_if_exists("${CMAKE_CURRENT_SOURCE_DIR}/Common") # JHS Also swap Meshing for Common here since DataCompression and Experimental still need VoroSpokes.*
#add_optional_subdirectory_if_exists("${CMAKE_CURRENT_SOURCE_DIR}/DataCompression")
#add_optional_subdirectory_if_exists("${CMAKE_CURRENT_SOURCE_DIR}/Experimental")
#add_optional_subdirectory_if_exists("${CMAKE_CURRENT_SOURCE_DIR}/AI")

vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")

