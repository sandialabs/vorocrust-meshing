# External Project for Exodus via SEACAS
#
# URL: https://github.com/gsjaardema/seacas
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(VOROCRUST_TPL_BUILD_ZLIB)
if(VOROCRUST_TPL_BUILD_ZLIB)

    include(ExternalProject)

    # Copy the TPL zip file over from the TPLs subdir
    # This is useful when the SNL proxy issues come up
    if(VOROCRUST_TPL_NO_DOWNLOAD)
        configure_file(${CMAKE_SOURCE_DIR}/tpls/zlib131.zip
                       ${TPL_BUILD_ROOT_DIR}/zlib/src/zlib131.zip
                       COPYONLY)
    endif()

    if(APPLE)
        set(ZLIB_CONFIGURE_EXTRA_ARGS -DCMAKE_MACOSX_RPATH:STRING=1)
        vc_message("ZLIB_CONFIGURE_EXTRA_ARGS = ${ZLIB_CONFIGURE_EXTRA_ARGS}")
    endif(APPLE)

    # 'installation' directory for SEACAS & its TPLs within the build directory
    set(TPL_ZLIB_INSTALL_DIR ${TPL_INSTALL_ROOT_DIR})
    vc_message("TPL_ZLIB_INSTALL_DIR = ${TPL_ZLIB_INSTALL_DIR}")

    # Toggle testing of ZLIB TPL
    set(TPL_ZLIB_ENABLE_TESTS ${VOROCRUST_ENABLE_TPL_TESTS})

    ExternalProject_Add(zlib
	URL https://www.zlib.net/zlib131.zip
	URL_MD5	ef76f7ebfd97778a6653b1d8413541c0

        PREFIX  ${TPL_BUILD_ROOT_DIR}/zlib

        CMAKE_CACHE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=${TPL_ZLIB_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_C_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            #-DCMAKE_EXE_LINKER_FLAGS:STRING=${VC_LINKFLAG_LIBDL}
            #-DCMAKE_SKIP_RPATH:BOOL=ON
            #-DCMAKE_SKIP_INSTALL_RPATH:BOOL=ON
            -DINSTALL_BIN_DIR:PATH=${TPL_ZLIB_INSTALL_DIR}/bin
            -DINSTALL_INC_DIR:PATH=${TPL_ZLIB_INSTALL_DIR}/include
            -DINSTALL_LIB_DIR:PATH=${TPL_ZLIB_INSTALL_DIR}/lib
            -DINSTALL_MAN_DIR:PATH=${TPL_ZLIB_INSTALL_DIR}/share/man
            -DINSTALL_PKGCONFIG_DIR:PATH=${TPL_ZLIB_INSTALL_DIR}/share/pkgconfig
            ${ZLIB_CONFIGURE_EXTRA_ARGS}

        TEST_BEFORE_INSTALL ${TPL_ZLIB_ENABLE_TESTS}

        # If logging is enabled then the output goes to file only.
        LOG_DOWNLOAD  ${VOROCRUST_ENABLE_TPL_LOG_DOWNLOAD}
        LOG_CONFIGURE ${VOROCRUST_ENABLE_TPL_LOG_CONF}
        LOG_BUILD     ${VOROCRUST_ENABLE_TPL_LOG_BUILD}
        LOG_TEST      ${VOROCRUST_ENABLE_TPL_LOG_TEST}
        LOG_INSTALL   ${VOROCRUST_ENABLE_TPL_LOG_INST}
    )

    set(TPL_ZLIB_INCLUDE_DIR ${TPL_ZLIB_INSTALL_DIR}/include CACHE FILEPATH "Path to ZLIB headers")
    set(TPL_ZLIB_LIBRARY_DIR ${TPL_ZLIB_INSTALL_DIR}/lib     CACHE FILEPATH "Path to ZLIB library")

    include_directories(${TPL_ZLIB_INCLUDE_DIR})

    set(TPL_ZLIB_PKGCONFIG "")
    if(WIN32)
        set(TPL_ZLIB_LIBRARIES ${TPL_ZLIB_LIBRARY_DIR}/zlibstaticd.lib)
    elseif(APPLE)
        set(TPL_ZLIB_LIBRARIES ${TPL_ZLIB_LIBRARY_DIR}/libz.a)
        set(TPL_ZLIB_PKGCONFIG ${TPL_ZLIB_INSTALL_DIR}/share/pkgconfig/zlib.pc)
    else()
        set(TPL_ZLIB_LIBRARIES ${TPL_ZLIB_LIBRARY_DIR}/libz.a)
        set(TPL_ZLIB_PKGCONFIG ${TPL_ZLIB_INSTALL_DIR}/share/pkgconfig/zlib.pc)
    endif()

    mark_as_advanced(FORCE TPL_ZLIB_INCLUDE_DIR)
    mark_as_advanced(FORCE TPL_ZLIB_LIBRARY_DIR)

    vc_message("TPL_ZLIB_PKGCONFIG   = ${TPL_ZLIB_PKGCONFIG}")
    vc_message("TPL_INSTALL_DIR      = ${TPL_INSTALL_DIR}")
    vc_message("TPL_LIBRARY_DIR      = ${TPL_LIBRARY_DIR}")
    vc_message("TPL_ZLIB_INCLUDE_DIR = ${TPL_ZLIB_INCLUDE_DIR}")
    vc_message("TPL_ZLIB_LIBRARIES   = ${TPL_ZLIB_LIBRARIES}")

endif()

vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
