#
# External Project for SEACAS TPL
#
# URL: https://github.com/gsjaardema/seacas
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(VOROCRUST_TPL_BUILD_SEACAS)

if(VOROCRUST_TPL_BUILD_SEACAS)

    include(ExternalProject)

    # Copy the TPL zip file over from the TPLs subdir
    # This is useful when the SNL proxy issues come up
    if(VOROCRUST_TPL_NO_DOWNLOAD)
        configure_file(${CMAKE_SOURCE_DIR}/tpls/seacas-2021-05-12.tar.gz
                       ${TPL_BUILD_ROOT_DIR}/seacas/src/seacas-2021-05-12.tar.gz
                       COPYONLY)
    endif()

    # Configure options if MPI is on or off.
    # (currently disabled because we need pnetcdf for this)
    #if(VOROCRUST_ENABLE_MPI)
    #    vc_message("Building SEACAS with MPI support.")
    #    set(SEACAS_ENABLE_MPI "ON")
    #else()
    #    vc_message("Building SEACAS without MPI support.")
    #    set(SEACAS_ENABLE_MPI "OFF")
    #endif()

    set(TPL_SEACAS_EXTRA_FLAGS "")
    if(APPLE)
        set(TPL_SEACAS_EXTRA_FLAGS "-DCMAKE_MACOSX_RPATH:BOOL=ON")
    else()
        set(TPL_SEACAS_EXTRA_FLAGS "")
    endif(APPLE)

    set(TPL_SEACAS_INSTALL_DIR ${TPL_INSTALL_ROOT_DIR})

    set(TPL_SEACAS_ENABLE_TESTS ${VOROCRUST_ENABLE_TPL_TESTS})

    ExternalProject_Add(seacas
        URL     https://github.com/gsjaardema/seacas/archive/refs/tags/v2021-05-12.tar.gz
        URL_MD5 f8391413fca7cf20721a16fe1de2b218

        DOWNLOAD_NAME seacas-2021-05-12.tar.gz

        TEST_BEFORE_INSTALL ON

        PREFIX ${TPL_BUILD_ROOT_DIR}/seacas

        CMAKE_CACHE_ARGS
            -DCMAKE_INSTALL_RPATH:PATH=${TPL_SEACAS_INSTALL_DIR}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPL_SEACAS_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_C_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DBUILD_SHARED_LIBS:BOOL=OFF
            -DSEACASProj_ENABLE_SEACASExodus:BOOL=ON
            -DSEACASProj_ENABLE_SEACASExodus_for:BOOL=OFF
            -DSEACASProj_ENABLE_SEACASExoIIv2for32:BOOL=OFF
            -DSEACASProj_ENABLE_TESTS:BOOL=${TPL_SEACAS_ENABLE_TESTS}
            -DSEACASExodus_ENABLE_STATIC:BOOL=ON
            -DSEACASExodus_ENABLE_SHARED:BOOL=OFF
            -DSEACASProj_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=ON
            -DSEACASProj_HIDE_DEPRECATED_CODE:BOOL=OFF
            -DSEACASProj_ENABLE_Fortran:BOOL=OFF
            -DTPL_ENABLE_Netcdf:BOOL=ON
            -DTPL_ENABLE_MPI:BOOL=OFF
            -DTPL_ENABLE_Pthread:BOOL=ON
            -DSEACASExodus_ENABLE_THREADSAFE:BOOL=ON
            -DNetCDF_ROOT:PATH=${TPL_NETCDF_INSTALL_DIR}
            -DPNetCDF_ROOT:PATH=${TPL_NETCDF_INSTALL_DIR}
            -DHDF5_DIR:PATH=${TPL_HDF5_CONFIG_DIR}
            #-DHDF5_ROOT:PATH=${TPL_INSTALL_ROOT_DIR}
            -DHDF5_NO_SYSTEM_PATHS:BOOL=ON
            -DTPL_ENABLE_DLlib:BOOL=OFF
            -DTPL_TENTATIVE_ENABLE_DLlib:BOOL=OFF
            ${TPL_SEACAS_EXTRA_FLAGS}

        TEST_BEFORE_INSTALL ${TPL_SEACAS_ENABLE_TESTS}

        # If logging is enabled then the output goes to file only.
        LOG_DOWNLOAD  ${VOROCRUST_ENABLE_TPL_LOG_DOWNLOAD}
        LOG_CONFIGURE ${VOROCRUST_ENABLE_TPL_LOG_CONF}
        LOG_BUILD     ${VOROCRUST_ENABLE_TPL_LOG_BUILD}
        LOG_TEST      ${VOROCRUST_ENABLE_TPL_LOG_TEST}
        LOG_INSTALL   ${VOROCRUST_ENABLE_TPL_LOG_INST}
    )


    ExternalProject_Add_StepDependencies(seacas build netcdf)

    set(TPL_SEACAS_INCLUDE_DIR ${TPL_SEACAS_INSTALL_DIR}/include CACHE FILEPATH "Path to SEACAS headers")
    set(TPL_SEACAS_LIBRARY_DIR ${TPL_SEACAS_INSTALL_DIR}/lib     CACHE FILEPATH "Path to SEACAS library")

    include_directories(${TPL_SEACAS_INCLUDE_DIR})

    # configure link libraries
    if(WIN32)
        if(NOT VOROCRUST_ENABLE_DEVTEST)
            # Guard against enabling SEACAS / Exodus on Windows
            vc_message("SEACAS currently can't be built on windows.")
            vc_message("- Set VOROCRUST_ENABLE_DEVTEST:BOOL=ON to disable this error (developers only!)")
            message(FATAL_ERROR "${TagPREFIX} ERROR: Building VoroCrust with SEACAS enabled on Windows is not supported.")
        endif()
        set(TPL_SEACAS_LIBRARIES ${TPL_SEACAS_LIBRARY_DIR}/libexodus.lib)
    elseif(APPLE)
       set(TPL_SEACAS_LIBRARIES ${TPL_SEACAS_LIBRARY_DIR}/libexodus.a)
    else()
        set(TPL_SEACAS_LIBRARIES ${TPL_SEACAS_LIBRARY_DIR}/libexodus.a)
    endif()

    mark_as_advanced(FORCE TPL_SEACAS_INCLUDE_DIR)
    mark_as_advanced(FORCE TPL_SEACAS_LIBRARY_DIR)

    vc_message("TPL_SEACAS_INSTALL_DIR  = ${TPL_SEACAS_INSTALL_DIR}")
    vc_message("TPL_SEACAS_LIBRARY_DIR  = ${TPL_SEACAS_LIBRARY_DIR}")
    vc_message("TPL_SEACAS_INCLUDE_DIR  = ${TPL_SEACAS_INCLUDE_DIR}")
    vc_message("TPL_SEACAS_LIBRARIES    = ${TPL_SEACAS_LIBRARIES}")

endif()


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
