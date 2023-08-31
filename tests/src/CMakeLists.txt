#
# VoroCrust Testing CMake
#
vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(CMAKE_SOURCE_DIR)
vc_message_var(CMAKE_CURRENT_SOURCE_DIR)
vc_message_var(CMAKE_BINARY_DIR)
vc_message_var(CMAKE_CURRENT_BINARY_DIR)


# If MPI is turned on, enable MPI specific tests.
if(VOROCRUST_ENABLE_MPI)

    #
    # Simple MPI hello-world app to make sure we're building with MPI
    #
    add_executable(test_mpi_hello test_mpi_hello.cpp)
    set_target_properties(test_mpi_hello
                          PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                          )
    target_include_directories(test_mpi_hello PRIVATE ${MPI_CXX_INCLUDE_PATH}
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_BINARY_DIR}/src
        ${CMAKE_BINARY_DIR}/src/Meshing)
    target_compile_options(test_mpi_hello PRIVATE ${MPI_CXX_COMPILE_FLAGS})
    target_link_libraries(test_mpi_hello ${MPI_CXX_LIBRARIES} libVCMesh)


    add_test(NAME test_mpi_hello
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 test_mpi_hello
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             )
endif(VOROCRUST_ENABLE_MPI)


# If Exodus is turned on, enable Exodus specific tests.
if(VOROCRUST_ENABLE_EXODUS)
        add_executable(test_exodus_simple test_exodus_simple.cpp)
        set_target_properties(test_exodus_simple
                              PROPERTIES
                              RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                              )
        target_include_directories(test_exodus_simple PRIVATE
            ${TPL_SEACAS_INCLUDES}
            ${CMAKE_BINARY_DIR}/src
            #${CMAKE_BINARY_DIR}/src/Common
            ${CMAKE_BINARY_DIR}/src/Meshing
            )

        target_link_libraries(test_exodus_simple
            ExodusIO
            ${TPL_SEACAS_LIBRARIES}
            ${TPL_NETCDF_LIBRARIES}
            ${TPL_HDF5_LIBRARIES}
            ${TPL_ZLIB_LIBRARIES}
            ${CMAKE_DL_LIBS}
            libVCMesh
            )

        add_test(NAME test_exodus_simple COMMAND test_exodus_simple)
endif(VOROCRUST_ENABLE_EXODUS)



# If Kokkos is turned on, enable Kokkos specific tests.
if(VOROCRUST_TPL_ENABLE_KOKKOS)

    #
    # Simple Kokkos hello-world app to make sure everything is working.
    #
    add_executable(test_kokkos_hello test_kokkos_hello.cpp)
    set_target_properties(test_kokkos_hello
                          PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                          )
    vorocrust_tpladd_kokkos(test_kokkos_hello)
    vorocrust_add_dependency(test_kokkos_hello libVCMesh)
    #target_link_libraries(test_kokkos_hello PUBLIC libVCMesh)
    #target_link_libraries(test_kokkos_hello Kokkos::kokkos libVCMesh)

    add_test(NAME test_kokkos_hello
             COMMAND test_kokkos_hello
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             )
endif()



# If KokkosKernels is turned on, enable Kokkos-kernels specific tests.
if(VOROCRUST_TPL_ENABLE_KOKKOSKERNELS)

    #
    # Simple Kokkos-Kernels example app to make sure everything is working.
    #
    add_executable(test_kokkoskernels_blas_abs test_kokkoskernels_blas_abs.cpp)
    set_target_properties(test_kokkoskernels_blas_abs
                          PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                          )
    vorocrust_tpladd_kokkos(test_kokkoskernels_blas_abs)
    vorocrust_tpladd_kokkoskernels(test_kokkoskernels_blas_abs)
    vorocrust_add_dependency(test_kokkoskernels_blas_abs libVCMesh)

    #target_link_libraries(test_kokkoskernels_hello PUBLIC libVCMesh)

    add_test(NAME test_kokkoskernels_blas_abs
             COMMAND test_kokkoskernels_blas_abs
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             )
endif()



#
# If LAPACK is turned on, enable LAPACK specific tests.
#
if(VOROCRUST_TPL_USE_LAPACK OR VOROCRUST_TPL_BUILD_OPENBLAS)

    add_executable(test_lapack test_lapack.cpp)
    set_target_properties(test_lapack
                          PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                          )

    vorocrust_tpladd_lapack(test_lapack)
    vorocrust_add_dependency(test_lapack libVCMesh)

    add_test(NAME test_lapack
             COMMAND test_lapack
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             )
endif()

vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
