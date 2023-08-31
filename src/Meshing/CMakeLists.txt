#
# VoroCrust Mesh CMake Script
#
vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(VOROCRUST_BUILD_MESH)

if(VOROCRUST_BUILD_MESH)

    set(target libVCMesh)

    # Sources
    file(GLOB VOROCRUST_MESH_SRCS *.cpp)
    file(GLOB VOROCRUST_MESH_HDRS *.h)
    set(VOROCRUST_MESH_EXTRA_HEADERS "")

    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    # OB VOROCRUST_COMMON_HDRS *.h)
    #set(VOROCRUST_COMMON_EXTRA_HEADERS "")

    # include_directories(${CMAKE_CURRENT_SOURCE_DIR})

    # Generate Version.h from Version.h.in
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Version.h.in" "${CMAKE_CURRENT_BINARY_DIR}/Version.h" @ONLY)
    list(APPEND VOROCRUST_MESH_EXTRA_HEADERS "${CMAKE_CURRENT_BINARY_DIR}/Version.h")

    # Generate Config.h from Config.h.in
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/Config.h" @ONLY)
    list(APPEND VOROCRUST_MESH_EXTRA_HEADERS "${CMAKE_CURRENT_BINARY_DIR}/Config.h")


    vc_message_var(VOROCRUST_MESH_EXTRA_HEADERS)
    vc_message_var(CMAKE_CURRENT_BINARY_DIR)
    vc_message_var(CMAKE_BINARY_DIR)

    # Include the binary dir to get the generated file(s).
    include_directories("${CMAKE_CURRENT_BINARY_DIR}")

    add_library(${target}
        ${VOROCRUST_MESH_HDRS}
        ${VOROCRUST_MESH_SRCS}
        ${VOROCRUST_MESH_EXTRA_HEADERS}
    )
    target_include_directories(${target}
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})
    # Add dependencies



    if(VOROCRUST_TPL_BUILD_ZLIB)
        add_dependencies(${target} zlib)
        target_link_libraries(${target} PUBLIC ${TPL_ZLIB_LIBRARIES})
    endif()
    if(VOROCRUST_TPL_BUILD_HDF5)
        add_dependencies(${target} hdf5)
    endif()
    if(VOROCRUST_TPL_BUILD_NETCDF)
        add_dependencies(${target} netcdf)
    endif()
    if(VOROCRUST_TPL_BUILD_SEACAS OR VOROCRUST_ENABLE_EXODUS)
        add_dependencies(${target} seacas)
    endif()

    # Exodus settings and dependencies
    if(VOROCRUST_ENABLE_EXODUS)
        add_library(ExodusIO OBJECT MeshingExodusIO.cpp MeshingExodusIO.h )

        target_include_directories(ExodusIO PUBLIC ${TPL_SEACAS_INCLUDE_DIR})
        target_compile_options(ExodusIO PUBLIC -c)
        target_link_libraries(ExodusIO PUBLIC
                ${TPL_SEACAS_LIBRARIES}
                ${TPL_NETCDF_LIBRARIES}
                ${TPL_HDF5_LIBRARIES}
                ${TPL_ZLIB_LIBRARIES}
                ${CMAKE_DL_LIBS}
                )

        target_include_directories(${target} PUBLIC ${TPL_SEACAS_INCLUDE_DIR})
        target_link_libraries(${target} PUBLIC
                ${TPL_SEACAS_LIBRARIES}
                ${TPL_NETCDF_LIBRARIES}
                ${TPL_HDF5_LIBRARIES}
                ${TPL_ZLIB_LIBRARIES}
                ${CMAKE_DL_LIBS}
                )

        link_directories(${TPL_SEACAS_LIB_DIR})

        add_dependencies(ExodusIO seacas)
        add_dependencies(${target} ExodusIO)
    endif(VOROCRUST_ENABLE_EXODUS)

    # vorocrust_add_dependency(${target} libVCCommon)

    # MPI settings and dependencies
    vorocrust_tpladd_mpi(${target})

    # Kokkos settings and dependencies
    vorocrust_tpladd_kokkos(${target})
    vorocrust_tpladd_kokkoskernels(${target})

    # LAPACK settings and dependencies
    vorocrust_tpladd_lapack(${target})


    # vc_mesh APP
    set(VOROCRUST_APPNAME_VCMESH "vc_mesh")
    set(VOROMESH_SOURCES "vc_mesh.cxx")
    set(VOROMESH_INCLUDES "") # "${CMAKE_CURRENT_BINARY_DIR}/Common")
    set(VOROMESH_LIBS  libVCMesh)
    vorocrust_add_app("${VOROCRUST_APPNAME_VCMESH}"
            VOROMESH_SOURCES
            VOROMESH_INCLUDES
            VOROMESH_LIBS
            )

    add_optional_subdirectory_if_exists("${CMAKE_CURRENT_SOURCE_DIR}/tests")

endif(VOROCRUST_BUILD_MESH)


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
