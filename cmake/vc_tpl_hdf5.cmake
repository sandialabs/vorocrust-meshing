# External Project for HDF5 TPL
#
# URL: https://github.com/gsjaardema/seacas
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

vc_message_var(VOROCRUST_TPL_BUILD_HDF5)
if(VOROCRUST_TPL_BUILD_HDF5)

    # Copy the HDF5 zip file over from the TPLs subdir
    # This is useful when the SNL proxy issues come up
    if(VOROCRUST_TPL_NO_DOWNLOAD)
        configure_file(${CMAKE_SOURCE_DIR}/tpls/hdf5-1.12.0.tar.gz
                       ${TPL_BUILD_ROOT_DIR}/hdf5/src/hdf5-1.12.0.tar.gz
                       COPYONLY)
    endif()

    include(ExternalProject)

    #if(APPLE)
    #    set(HDF5_CONFIGURE_EXTRA_ARGS -DCMAKE_MACOSX_RPATH:STRING=1)
    #    vc_message("HDF5_CONFIGURE_EXTRA_ARGS = ${HDF5_CONFIGURE_EXTRA_ARGS}")
    #endif(APPLE)
    #set(HDF5_CONFIGURE_WARNING_ARGS "-Wno-undef -Wno-cast-qual -Wno-discarded-qualifiers")

    # 'installation' directory for SEACAS & its TPLs within the build directory
    set(TPL_HDF5_INSTALL_DIR ${TPL_INSTALL_ROOT_DIR})
    vc_message("TPL_HDF5_INSTALL_DIR  = ${TPL_HDF5_INSTALL_DIR}")

    # HDF5 Linker Flags
    set(VC_HDF5_EXE_LINKER_FLAGS "")
    if(NOT WIN32)
        set(VC_HDF5_EXE_LINKER_FLAGS "${VC_LINKFLAG_LIBDL}")
    endif(NOT WIN32)

    # Toggle testing of HDF5 TPL
    set(TPL_HDF5_ENABLE_TESTS ${VOROCRUST_ENABLE_TPL_TESTS})

    set(TPL_HDF5_EXTRA_FLAGS "")

    ExternalProject_Add(hdf5
        URL     https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz
        URL_MD5 9e22217d22eb568e09f0cc15fb641d7c

        PREFIX  ${TPL_BUILD_ROOT_DIR}/hdf5

        CMAKE_CACHE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=${TPL_HDF5_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_C_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_EXE_LINKER_FLAGS:STRING=${VC_HDF5_EXE_LINKER_FLAGS}
            -DBUILD_SHARED_LIBS:BOOL=OFF
            -DBUILD_STATIC_LIBS:BOOL=ON
            -DBUILD_TESTING:BOOL=${TPL_HDF5_ENABLE_TESTS}
            -DHDF5_ENABLE_PARALLEL:BOOL=OFF
            -DHDF5_BUILD_FORTRAN:BOOL=OFF
            -DHDF5_BUILD_HL_LIB:BOOL=ON
            -DHDF5_BUILD_EXAMPLES:BOOL=OFF
            -DHDF5_BUILD_TOOLS:BOOL=ON
            -DHDF5_DISABLE_COMPILER_WARNINGS:BOOL=ON
            -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
            -DZLIB_INCLUDE_DIR:PATH=${TPL_ZLIB_INCLUDE_DIR}
            -DZLIB_LIBRARY_RELEASE:FILEPATH=${TPL_ZLIB_LIBRARIES}
            ${TPL_HDF5_EXTRA_FLAGS}

        TEST_BEFORE_INSTALL ${TPL_HDF5_ENABLE_TESTS}
        #TEST_COMMAND "${CMAKE_CTEST_COMMAND} -j ${NUM_CORES}"

        # If logging is enabled then the output goes to file only.
        LOG_DOWNLOAD  ${VOROCRUST_ENABLE_TPL_LOG_DOWNLOAD}
        LOG_CONFIGURE ${VOROCRUST_ENABLE_TPL_LOG_CONF}
        LOG_BUILD     ${VOROCRUST_ENABLE_TPL_LOG_BUILD}
        LOG_TEST      ${VOROCRUST_ENABLE_TPL_LOG_TEST}
        LOG_INSTALL   ${VOROCRUST_ENABLE_TPL_LOG_INST}
    )

    # Add dependency
    ExternalProject_Add_StepDependencies(hdf5 build zlib)


    set(TPL_HDF5_INCLUDE_DIR ${TPL_HDF5_INSTALL_DIR}/include CACHE FILEPATH "Path to HDF5 headers")
    set(TPL_HDF5_LIBRARY_DIR ${TPL_HDF5_INSTALL_DIR}/lib     CACHE FILEPATH "Path to HDF5 library")

    set(TPL_HDF5_CMAKE_CONFIG_DIR ${TPL_HDF5_INSTALL_DIR}/share/cmake/hdf5 CACHE FILEPATH "Path to HDF5 cmake config")

    include_directories(${TPL_HDF5_INCLUDE_DIR})

    # configure link libraries
    set(TPL_HDF5_SUFFIX "")

    if(WIN32)
        set(TPL_HDF5_SUFFIX "_D")
        set(TPL_HDF5_LIBHDF5       ${TPL_HDF5_LIBRARY_DIR}/libhdf5${TPL_HDF5_SUFFIX}.lib)
        set(TPL_HDF5_LIBHDF5_CPP   ${TPL_HDF5_LIBRARY_DIR}/libhdf5_cpp${TPL_HDF5_SUFFIX}.lib)
        set(TPL_HDF5_LIBHDF5_HL    ${TPL_HDF5_LIBRARY_DIR}/libhdf5_hl${TPL_HDF5_SUFFIX}.lib)
        set(TPL_HDF5_LIBHDF5_TOOLS ${TPL_HDF5_LIBRARY_DIR}/libhdf5_tools${TPL_HDF5_SUFFIX}.lib)
        set(TPL_HDF5_CONFIG_DIR    ${TPL_HDF5_INSTALL_DIR}/cmake/hdf5)
    elseif(APPLE)
        if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            set(TPL_HDF5_SUFFIX "_debug")
        endif()
        set(TPL_HDF5_LIBHDF5       ${TPL_HDF5_LIBRARY_DIR}/libhdf5${TPL_HDF5_SUFFIX}.a)
        set(TPL_HDF5_LIBHDF5_CPP   ${TPL_HDF5_LIBRARY_DIR}/libhdf5_cpp${TPL_HDF5_SUFFIX}.a)
        set(TPL_HDF5_LIBHDF5_HL    ${TPL_HDF5_LIBRARY_DIR}/libhdf5_hl${TPL_HDF5_SUFFIX}.a)
        set(TPL_HDF5_LIBHDF5_TOOLS ${TPL_HDF5_LIBRARY_DIR}/libhdf5_tools${TPL_HDF5_SUFFIX}.a)
        set(TPL_HDF5_CONFIG_DIR    ${TPL_HDF5_INSTALL_DIR}/share/cmake/hdf5)
    else()
        if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            set(TPL_HDF5_SUFFIX "_debug")
        endif()
        set(TPL_HDF5_LIBHDF5       ${TPL_HDF5_LIBRARY_DIR}/libhdf5${TPL_HDF5_SUFFIX}.a)
        set(TPL_HDF5_LIBHDF5_CPP   ${TPL_HDF5_LIBRARY_DIR}/libhdf5_cpp${TPL_HDF5_SUFFIX}.a)
        set(TPL_HDF5_LIBHDF5_HL    ${TPL_HDF5_LIBRARY_DIR}/libhdf5_hl${TPL_HDF5_SUFFIX}.a)
        set(TPL_HDF5_LIBHDF5_TOOLS ${TPL_HDF5_LIBRARY_DIR}/libhdf5_tools${TPL_HDF5_SUFFIX}.a)
        set(TPL_HDF5_CONFIG_DIR    ${TPL_HDF5_INSTALL_DIR}/share/cmake/hdf5)
    endif()

    # wcmclen: Notionally this should work but it didn't in my testing so something
    # is probably missing.
    #add_custom_command(TARGET hdf5
    #    PRE_BUILD
    #    COMMAND ${CMAKE_COMMAND} -E copy
    #            ${CMAKE_SOURCE_DIR}/tpls/hdf5-1.12.0.zip
    #            ${TPL_BUILD_ROOT_DIR}/hdf5/src/hdf5-1.12.0.zip
    #    )


    set(TPL_HDF5_LIBRARIES
        ${TPL_HDF5_LIBHDF5_TOOLS}
        ${TPL_HDF5_LIBHDF5_HL}
        ${TPL_HDF5_LIBHDF5_CPP}
        ${TPL_HDF5_LIBHDF5}
    )

    mark_as_advanced(FORCE TPL_HDF5_INCLUDE_DIR)
    mark_as_advanced(FORCE TPL_HDF5_LIBRARY_DIR)
    mark_as_advanced(FORCE TPL_HDF5_CMAKE_CONFIG_DIR)

    vc_message("TPL_HDF5_CONFIG_DIR   = ${TPL_HDF5_CONFIG_DIR}")
    vc_message("TPL_HDF5_INSTALL_DIR  = ${TPL_HDF5_INSTALL_DIR}")
    vc_message("TPL_HDF5_INCLUDE_DIR  = ${TPL_HDF5_INCLUDE_DIR}")
    vc_message("TPL_HDF5_LIBRARY_DIR  = ${TPL_HDF5_LIBRARY_DIR}")
    vc_message("TPL_HDF5_LIBRARIES    = ${TPL_HDF5_LIBRARIES}")

endif()


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
