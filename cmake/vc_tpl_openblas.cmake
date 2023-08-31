include_guard(GLOBAL)
# External Project for OpenBLAS
#
# URL: https://github.com/xianyi/OpenBLAS/
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(VOROCRUST_TPL_BUILD_OPENBLAS)
vc_message_var(VOROCRUST_USE_BUILT_LAPACK)

if(VOROCRUST_USE_BUILT_LAPACK OR VOROCRUST_TPL_BUILD_OPENBLAS)

    # This banner should be removed once we figure out how to get OPENBLAS working
    # or abandon our attempt to pull it in as a TPL.
    vc_message("${Red}===============================================================${ColoReset}")
    vc_message("${Red}WARNING - OpenBLAS Integration is still a Work-In-Progress${ColoReset}")
    vc_message("${Red}          It is unlikely to work at this time${ColoReset}")
    vc_message("${Red}===============================================================${ColoReset}")

    include(ExternalProject)

    set(OPENBLAS_CONFIGURE_EXTRA_ARGS "")
    if(APPLE)
        set(OPENBLAS_CONFIGURE_EXTRA_ARGS "-DCMAKE_MACOSX_RPATH:STRING=1 ${OPENBLAS_CONFIGURE_EXTRA_ARGS}")
    endif()
    if(WIN32)
        #set(OPENBLAS_CONFIGURE_EXTRA_ARGS "-DNOFORTRAN=1 -DDYNAMIC_ARCH:BOOL=OFF -DVS_WINRT_COMPONENT:BOOL=ON")
        set(OPENBLAS_CONFIGURE_EXTRA_ARGS "${OPENBLAS_CONFIGURE_EXTRA_ARGS}")
    endif()
    vc_message_var(OPENBLAS_CONFIGURE_EXTRA_ARGS)

    # 'installation' directory for OpenBLAS & its TPLs within the build directory
    set(TPL_OPENBLAS_INSTALL_DIR ${TPL_INSTALL_ROOT_DIR})
    vc_message_var(TPL_OPENBLAS_INSTALL_DIR)

    # Toggle testing of OpenBLAS TPL
    set(TPL_OPENBLAS_ENABLE_TESTS ${VOROCRUST_ENABLE_TPL_TESTS})

    ExternalProject_Add(openblas
        URL     https://github.com/xianyi/OpenBLAS/releases/download/v0.3.18/OpenBLAS-0.3.18.zip

        URL_MD5 0ebf2e1ddc491f37be26bea4e0d1239a

        PREFIX  ${TPL_BUILD_ROOT_DIR}/OpenBLAS

        CMAKE_CACHE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=${TPL_OPENBLAS_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DBUILD_TESTING:BOOL=OFF
            -DCMAKE_C_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_INSTALL_LIBDIR:STRING=lib
            -DCMAKE_INSTALL_LIBEXEC:STRING=libexec
            -DCMAKE_INSTALL_INCLUDEDIR:STRING=include
            -DCMAKE_INSTALL_BINDIR:STRING=bin
            -DCMAKE_INSTALL_DATAROOTDIR:STRING=share
            ${OPENBLAS_CONFIGURE_EXTRA_ARGS}

        TEST_BEFORE_INSTALL ${TPL_OPENBLAS_ENABLE_TESTS}

        # If logging is enabled then the output goes to file only.
        LOG_DOWNLOAD  ${VOROCRUST_ENABLE_TPL_LOG_DOWNLOAD}
        LOG_CONFIGURE ${VOROCRUST_ENABLE_TPL_LOG_CONF}
        LOG_BUILD     ${VOROCRUST_ENABLE_TPL_LOG_BUILD}
        LOG_TEST      ${VOROCRUST_ENABLE_TPL_LOG_TEST}
        LOG_INSTALL   ${VOROCRUST_ENABLE_TPL_LOG_INST}
    )

    set(TPL_OPENBLAS_INCLUDE_DIR ${TPL_OPENBLAS_INSTALL_DIR}/include CACHE FILEPATH "Path to OpenBLAS headers")
    set(TPL_OPENBLAS_LIBRARY_DIR ${TPL_OPENBLAS_INSTALL_DIR}/lib     CACHE FILEPATH "Path to OpenBLAS library")

    include_directories(${TPL_OPENBLAS_INCLUDE_DIR})

    # configure link libraries
    if(WIN32)
        if(NOT VOROCRUST_ENABLE_DEVTEST)
            # Guard against enabling OpenBLAS on Windows
            vc_message("OpenBLAS currently can't be built on windows.")
            vc_message("- Set VOROCRUST_ENABLE_DEVTEST:BOOL=ON to disable this error (developers only!)")
            message(FATAL_ERROR "${TagPREFIX} ERROR: Building VoroCrust with OpenBLAS enabled on Windows is not supported.")
        endif()
        set(TPL_OPENBLAS_LIBRARIES ${TPL_OPENBLAS_LIBRARY_DIR}/openblas.lib)
    elseif(APPLE)
       set(TPL_OPENBLAS_LIBRARIES ${TPL_OPENBLAS_LIBRARY_DIR}/libopenblas.a)
    else()
        set(TPL_OPENBLAS_LIBRARIES ${TPL_OPENBLAS_LIBRARY_DIR}/libopenblas.a)
    endif()

    mark_as_advanced(FORCE TPL_OPENBLAS_INCLUDE_DIR)
    mark_as_advanced(FORCE TPL_OPENBLAS_LIBRARY_DIR)

    vc_message_var(TPL_OPENBLAS_INSTALL_DIR)
    vc_message_var(TPL_OPENBLAS_LIBRARY_DIR)
    vc_message_var(TPL_OPENBLAS_INCLUDE_DIR)
    vc_message_var(TPL_OPENBLAS_LIBRARIES)

endif()


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
