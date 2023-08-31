# External Project for NetCDF TPL
#
# URL: https://github.com/gsjaardema/seacas
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

if(VOROCRUST_TPL_BUILD_NETCDF)
    vc_message("VOROCRUST_TPL_BUILD_NETCDF = ON")

    include(ExternalProject)

    # Copy the TPL zip file over from the TPLs subdir
    # This is useful when the SNL proxy issues come up
    if(VOROCRUST_TPL_NO_DOWNLOAD)
        configure_file(${CMAKE_SOURCE_DIR}/tpls/netcdf-c-4.8.1.zip
                       ${TPL_BUILD_ROOT_DIR}/netcdf/src/netcdf-c-4.8.1.zip
                       COPYONLY)
    endif()

    # 'installation' directory for SEACAS & its TPLs within the build directory
    set(TPL_NETCDF_INSTALL_DIR ${TPL_INSTALL_ROOT_DIR})

    # TODO: We may need to pull in libdl for windows builds?
    set(TPL_NETCDF_NC_EXTRA_DEPS "")
    if(NOT WIN32)
        set(TPL_NETCDF_NC_EXTRA_DEPS "${VC_LINKFLAG_LIBDL}")
    endif(NOT WIN32)

    set(TPL_NETCDF_EXTRA_FLAGS "")

    # Toggle testing of NetCDF TPL
    set(TPL_NETCDF_ENABLE_TESTS ${VOROCRUST_ENABLE_TPL_TESTS})

    ExternalProject_Add(netcdf
        URL     https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.8.1.zip
        URL_MD5 a8d6933f2b13d00bb8c1c7eb557cb904

        DOWNLOAD_NAME netcdf-c-4.8.1.zip

        PREFIX ${TPL_BUILD_ROOT_DIR}/netcdf

        CMAKE_CACHE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=${TPL_NETCDF_INSTALL_DIR}
            -DCMAKE_INSTALL_LIBDIR:PATH=${TPL_NETCDF_INSTALL_DIR}/lib
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS:BOOL=OFF
            -DBUILD_UTILITIES:BOOL=ON
            -DENABLE_EXAMPLES:BOOL=OFF
            -DENABLE_TESTS:BOOL=${TPL_NETCDF_ENABLE_TESTS}
            -DBUILD_TESTING:BOOL=${TPL_NETCDF_ENABLE_TESTS}
            -DENABLE_UNIT_TESTS:BOOL=${TPL_NETCDF_ENABLE_TESTS}
            -DCMAKE_C_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DENABLE_MMAP:BOOL=ON
            -DENABLE_DAP:BOOL=OFF
            -DENABLE_V2_API:BOOL=OFF
            -DENABLE_CDF5:BOOL=ON
            -DENABLE_NETCDF_4:BOOL=ON
            -DENABLE_CONVERSION_WARNINGS:BOOL=OFF
            # We should just need HDF_DIR, but if HDF5 is on the system it
            # may not be sufficient.
            -DHDF5_DIR:PATH=${TPL_HDF5_CONFIG_DIR}
            # If HDF5_DIR isn't sufficient, we'll need these:
            -DHDF5_INCLUDE_DIR:PATH=${TPL_HDF5_INCLUDE_DIR}
            -DHDF5_C_INCLUDE_DIR:PATH=${TPL_HDF5_INCLUDE_DIR}
            -DHDF5_C_LIBRARY:FILEPATH=${TPL_HDF5_LIBHDF5}
            -DHDF5_C_COMPILER_EXECUTABLE:FILEPATH=${TPL_HDF5_INSTALL_DIR}/bin/h5cc
            -DHDF5_DIFF_EXECUTABLE:FILEPATH=${TPL_HDF5_INSTALL_DIR}/bin/h5diff
            -DHDF5_LIBRARIES:FILEPATH=${TPL_HDF5_LIBRARIES}
            -DHDF5_hdf5_LIBRARY_RELEASE:FILEPATH=${TPL_HDF5_LIBHDF5}
            -DHDF5_hdf5_hl_LIBRARY_RELEASE:FILEPATH=${TPL_HDF5_LIBHDF5_HL}
            -DZLIB_LIBRARY:FILEPATH=${TPL_ZLIB_LIBRARIES}
            -DZLIB_INCLUDE_DIR:PATH=${TPL_ZLIB_INCLUDE_DIR}
            -DNC_EXTRA_DEPS:STRING=${TPL_NETCDF_NC_EXTRA_DEPS}
            ${TPL_NETCDF_EXTRA_FLAGS}


        TEST_BEFORE_INSTALL ${TPL_NETCDF_ENABLE_TESTS}

        # If logging is enabled then the output goes to file only.
        LOG_DOWNLOAD  ${VOROCRUST_ENABLE_TPL_LOG_DOWNLOAD}
        LOG_CONFIGURE ${VOROCRUST_ENABLE_TPL_LOG_CONF}
        LOG_BUILD     ${VOROCRUST_ENABLE_TPL_LOG_BUILD}
        LOG_TEST      ${VOROCRUST_ENABLE_TPL_LOG_TEST}
        LOG_INSTALL   ${VOROCRUST_ENABLE_TPL_LOG_INST}
    )

    ExternalProject_Add_StepDependencies(netcdf build hdf5)

    set(TPL_NETCDF_INCLUDE_DIR ${TPL_NETCDF_INSTALL_DIR}/include CACHE FILEPATH "Path to NetCDF headers")
    set(TPL_NETCDF_LIBRARY_DIR ${TPL_NETCDF_INSTALL_DIR}/lib     CACHE FILEPATH "Path to NetCDF library")

    include_directories(${TPL_NETCDF_INCLUDE_DIR})

    set(TPL_NETCDF_CONFIG_DIR "")
    set(TPL_NETCDF_PKGCONFIG  "")

    if(WIN32)
        set(TPL_NETCDF_LIBRARIES  ${TPL_NETCDF_LIBRARY_DIR}/netcdf.lib)
    elseif(APPLE)
       set(TPL_NETCDF_LIBRARIES   ${TPL_NETCDF_LIBRARY_DIR}/libnetcdf.a)
        set(TPL_NETCDF_CONFIG_DIR ${TPL_NETCDF_INSTALL_DIR}/lib/cmake/netCDF)
        set(TPL_NETCDF_PKGCONFIG  ${TPL_NETCDF_INSTALL_DIR}/lib/pkgconfig/netcdf.pc)

    else()
        set(TPL_NETCDF_LIBRARIES  ${TPL_NETCDF_LIBRARY_DIR}/libnetcdf.a)
        set(TPL_NETCDF_CONFIG_DIR ${TPL_NETCDF_INSTALL_DIR}/lib/cmake/netCDF)
        set(TPL_NETCDF_PKGCONFIG  ${TPL_NETCDF_INSTALL_DIR}/lib/pkgconfig/netcdf.pc)
    endif()

    mark_as_advanced(FORCE TPL_NETCDF_INCLUDE_DIR)
    mark_as_advanced(FORCE TPL_NETCDF_LIBRARY_DIR)

    vc_message("TPL_NETCDF_INSTALL_DIR = ${TPL_NETCDF_INSTALL_DIR}")
    vc_message("TPL_NETCDF_LIBRARY_DIR = ${TPL_NETCDF_LIBRARY_DIR}")
    vc_message("TPL_NETCDF_INCLUDE_DIR = ${TPL_NETCDF_INCLUDE_DIR}")
    vc_message("TPL_NETCDF_LIBRARIES   = ${TPL_NETCDF_LIBRARIES}")

endif()


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
