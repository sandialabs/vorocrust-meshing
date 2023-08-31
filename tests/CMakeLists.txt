include(CTest)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")


#==============================================
# Testing
#==============================================
vc_message_var(PROJECT_BINARY_DIR)
vc_message_var(PROJECT_SOURCE_DIR)


#==============================================
# Functions & Helpers
#==============================================



#==============================================
# Test Definitions
#==============================================


#====================================
# Sample Dataset Tests
#====================================


#---------------------------
# BEGIN Tests for `vc_mesh`
#---------------------------
if(VOROCRUST_BUILD_MESH)

    get_property(VOROCRUST_APPNAME_VCMESH GLOBAL PROPERTY "VOROCRUST_APPNAME_vc_mesh")
    vc_message("--")
    vc_message("-- Adding tests for ${VOROCRUST_APPNAME_VCMESH}")
    vc_message("--")

    #-----------
    # OneBox
    #-----------
    # A basic test using the OneBox example.
    set(OneBoxFiles OneBox.obj)

    vorocrust_add_simple_test("OneBox" ${VOROCRUST_APPNAME_VCMESH} "OneBox" "vc.in" "${OneBoxFiles}")
    set_tests_properties("OneBox" PROPERTIES TIMEOUT 1800)

    vorocrust_add_simple_test("OneBox_vcg" ${VOROCRUST_APPNAME_VCMESH} "OneBox" "vc_vcg.in" "${OneBoxFiles}")
    set_tests_properties("OneBox_vcg" PROPERTIES TIMEOUT 7200)

    if(VOROCRUST_ENABLE_EXODUS)
        vorocrust_add_simple_test("OneBox_exodus" ${VOROCRUST_APPNAME_VCMESH} "OneBox" "vc_exodus.in" "${OneBoxFiles}")
        set_tests_properties("OneBox_exodus" PROPERTIES TIMEOUT 7200)
    endif()

    if(OPENMP_FOUND)
        vorocrust_add_simple_test("OneBox_omp_4t" ${VOROCRUST_APPNAME_VCMESH} "OneBox" "vc_omp_4t.in" "${OneBoxFiles}")
        set_tests_properties("OneBox_omp_4t" PROPERTIES TIMEOUT 3600)
    endif()

    #-----------
    # Bird
    #-----------
    if(VOROCRUST_ENABLE_LONG_TESTS)
        set(testBirdFiles bird.obj)
        vorocrust_add_simple_test("Bird" ${VOROCRUST_APPNAME_VCMESH} "Bird" "vc.in" "${testBirdFiles}")
        set_tests_properties("Bird" PROPERTIES TIMEOUT 3600)
    endif()

    #-----------
    # Bracket
    #-----------
    if(VOROCRUST_ENABLE_LONG_TESTS)
        set(testBracketFiles bracket.obj)
        vorocrust_add_simple_test("Bracket" ${VOROCRUST_APPNAME_VCMESH} "Bracket" "vc.in" "${testBracketFiles}")
        set_tests_properties("Bracket" PROPERTIES TIMEOUT 7200)
    endif()

    #-----------
    # Non Manifold Box
    #-----------
    set(testNonManifoldBoxFiles non_manifold_box.obj)
    vorocrust_add_simple_test("Non_Manifold_Box" ${VOROCRUST_APPNAME_VCMESH} "Non_Manifold_Box" "vc.in" "${testNonManifoldBoxFiles}")
    set_tests_properties("Non_Manifold_Box" PROPERTIES TIMEOUT 3600)

    if(OPENMP_FOUND)
        vorocrust_add_simple_test("Non_Manifold_Box_omp_4t" ${VOROCRUST_APPNAME_VCMESH} "Non_Manifold_Box" "vc_omp_4t.in" "${testNonManifoldBoxFiles}")
        set_tests_properties("Non_Manifold_Box_omp_4t" PROPERTIES TIMEOUT 3600)
    endif()

    #-----------
    # Torus
    #-----------
    set(testTorusFiles torus.obj)
    # TEST_NAME DATASET_NAME VC_IN FILE_LIST
    vorocrust_add_simple_test("Torus" ${VOROCRUST_APPNAME_VCMESH} "Torus" "vc.in" "${testTorusFiles}")
    set_tests_properties("Torus" PROPERTIES TIMEOUT 3600) #3600

    #-----------
    # Monitoring Points (combined running w/ imposing pts and python script)
    #-----------
    set(OneBoxMPFiles OneBox.obj monitoring_points.csv check_points_in_mesh.py run_series.sh)
    vorocrust_add_monitoring_points_test("OneBox_monitoring_points" ${VOROCRUST_APPNAME_VCMESH} "OneBox_monitoring_points" "vc.in" "${OneBoxMPFiles}")
    set_tests_properties("OneBox_monitoring_points" PROPERTIES TIMEOUT 7200)

endif(VOROCRUST_BUILD_MESH)
#---------------------------
# END Tests for `vc_mesh`
#---------------------------


add_subdirectory(src)

vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
