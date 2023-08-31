#
# VoroCrust Testing CMake
#
vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(CMAKE_SOURCE_DIR)
vc_message_var(CMAKE_CURRENT_SOURCE_DIR)
vc_message_var(CMAKE_BINARY_DIR)
vc_message_var(CMAKE_CURRENT_BINARY_DIR)

#
# BEGIN: Test Mesh/OptionParserVCMesh
#
add_executable(test_Mesh_OptionParserVCMesh test_Mesh_OptionParserVCMesh.cpp)
target_link_libraries(test_Mesh_OptionParserVCMesh PUBLIC libVCMesh) # libVCCommon libVCMesh) # JHS 081122 removal of libVCCommon

set_target_properties(test_Mesh_OptionParserVCMesh
                      PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
                      )

add_test(NAME test_Mesh_OptionParserVCMesh_help  COMMAND test_Mesh_OptionParserVCMesh -h)
add_test(NAME test_Mesh_OptionParserVCMesh_vc    COMMAND test_Mesh_OptionParserVCMesh -vc vc.in)
#
# END: Test Mesh/OptionParserVCMesh
#

vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")

