# External Project for TuckerMPI
#
# URL: https://gitlab.com/tensors/TuckerMPI
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(VOROCRUST_ENABLE_TUCKERMPI)

if(VOROCRUST_ENABLE_TUCKERMPI)

    # This banner should be removed once we figure out how to get TuckerMPI working
    # or abandon our attempt to pull it in as a TPL.
    vc_message("${Red}===============================================================${ColoReset}")
    vc_message("${Red}WARNING - TUCKERMPI Integration is still a Work-In-Progress${ColoReset}")
    vc_message("${Red}          It is unlikely to work at this time${ColoReset}")
    vc_message("${Red}===============================================================${ColoReset}")

    include(ExternalProject)

    if(APPLE)
        set(TUCKERMPI_CONFIGURE_EXTRA_ARGS -DCMAKE_MACOSX_RPATH:STRING=1)
        vc_message("-- ${Magenta}TUCKERMPI_CONFIGURE_EXTRA_ARGS = ${TUCKERMPI_CONFIGURE_EXTRA_ARGS}${ColorReset}")
    endif(APPLE)

    # 'installation' directory for SEACAS & its TPLs within the build directory
    set(TPL_TUCKERMPI_INSTALL_DIR ${TPL_INSTALL_ROOT_DIR})
    vc_message("-- ${Magenta}TPL_TUCKERMPI_INSTALL_DIR = ${TPL_TUCKERMPI_INSTALL_DIR}${ColorReset}")

    # Toggle testing of TuckerMPI TPL
    set(TPL_TUCKERMPI_ENABLE_TESTS ${VOROCRUST_ENABLE_TPL_TESTS})

    ExternalProject_Add(tuckermpi
        URL     https://gitlab.com/tensors/TuckerMPI/-/archive/master/TuckerMPI-master.zip

        #URL_MD5 1c9f62f0778697a09d36121ead88e08e

        #DOWNLOAD_NAME TuckerMPI-master.zip

        PREFIX  ${TPL_BUILD_ROOT_DIR}/TuckerMPI

        SOURCE_SUBDIR src

        CMAKE_CACHE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=${TPL_TUCKERMPI_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DBUILD_TESTING:BOOL=OFF
            -DCMAKE_C_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            #-DCMAKE_EXE_LINKER_FLAGS:STRING=${VC_LINKFLAG_LIBDL}
            #-DCMAKE_SKIP_RPATH:BOOL=ON
            #-DCMAKE_SKIP_INSTALL_RPATH:BOOL=ON
            -DINSTALL_BIN_DIR:PATH=${TPL_TUCKERMPI_INSTALL_DIR}/bin
            -DINSTALL_INC_DIR:PATH=${TPL_TUCKERMPI_INSTALL_DIR}/include
            -DINSTALL_LIB_DIR:PATH=${TPL_TUCKERMPI_INSTALL_DIR}/lib
            -DINSTALL_MAN_DIR:PATH=${TPL_TUCKERMPI_INSTALL_DIR}/share/man
            -DINSTALL_PKGCONFIG_DIR:PATH=${TPL_TUCKERMPI_INSTALL_DIR}/share/pkgconfig
            ${TUCKERMPI_CONFIGURE_EXTRA_ARGS}

        TEST_BEFORE_INSTALL ${TPL_TUCKERMPI_ENABLE_TESTS}

        # If logging is enabled then the output goes to file only.
        LOG_DOWNLOAD  ${VOROCRUST_ENABLE_TPL_LOG_DOWNLOAD}
        LOG_CONFIGURE ${VOROCRUST_ENABLE_TPL_LOG_CONF}
        LOG_BUILD     ${VOROCRUST_ENABLE_TPL_LOG_BUILD}
        LOG_TEST      ${VOROCRUST_ENABLE_TPL_LOG_TEST}
        LOG_INSTALL   ${VOROCRUST_ENABLE_TPL_LOG_INST}
    )

endif()

vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
