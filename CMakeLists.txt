cmake_minimum_required(VERSION 3.12.0)
# 3.12 is required due to Kokkos dependency.
# We may need higher though because 3.12 has issues with the installation targets
# but we have a guard around that currently.
project(VoroCrust-Meshing)


# Set C++11 Standard as the minimum required.
set(CMAKE_CXX_STANDARD "14" CACHE STRING "Set by VoroCrust.")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Set by VoroCrust")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Set by VoroCrust")


#========================================
# Import cmake common helper file(s)
#========================================
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include(vc_common)
include(vc_functions)



# Set Build Type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE
        STRING "Choose build type, options are: None Debug Release RelWithDebInfo MinSizeRel"
        FORCE)
endif()



#========================================
# VoroCrust Project Options
#========================================
option(VOROCRUST_VERBOSE_CMAKE "Enable verbosity in CMake output" OFF)

option(VOROCRUST_BUILD_MESH            "Build the VC_Mesh app"             ON)

option(VOROCRUST_ENABLE_MPI    "Enable MPI Compilation"                      OFF)
option(VOROCRUST_ENABLE_OPENMP "Enable OpenMP for VoroCrust"                 ON )
option(VOROCRUST_ENABLE_EXODUS     "Enable Exodus file writing in VoroCrust" OFF)

option(VOROCRUST_ENABLE_LONG_TESTS "Enable long duration tests" OFF)

option(VOROCRUST_TPL_ENABLE_KOKKOS "Enable Kokkos library support"                OFF)
option(VOROCRUST_TPL_ENABLE_KOKKOSKERNELS "Enable Kokkos Kernels library support" OFF)
set(VOROCRUST_TPL_KOKKOS_VERSION "3.5.00" CACHE STRING "Kokkos version. This can be a tag, branch, or SHA1")
mark_as_advanced(VOROCRUST_TPL_KOKKOS_VERSION)
# Note: Kokkos and Kokkos-Kernels versions are synchronized

option(VOROCRUST_TPL_USE_LAPACK     "Use LAPACK"                            OFF)
option(VOROCRUST_TPL_BUILD_OPENBLAS "Download and build our own OpenBLAS"   OFF)

option(VOROCRUST_ENABLE_TPL_LOG_DOWNLOAD "Make ExternalPackage TPLs send download output to logs"  ON)
option(VOROCRUST_ENABLE_TPL_LOG_CONF     "Make ExternalPackage TPLs send configure output to logs" OFF)
option(VOROCRUST_ENABLE_TPL_LOG_BUILD    "Make ExternalPackage TPLs send build output to logs"     OFF)
option(VOROCRUST_ENABLE_TPL_LOG_TEST     "Make ExternalPackage TPLs send test output to logs"      OFF)
option(VOROCRUST_ENABLE_TPL_LOG_INST     "Make ExternalPackage TPLs send install output to logs"   ON)
mark_as_advanced(VOROCRUST_ENABLE_TPL_LOG_DOWNLOAD)
mark_as_advanced(VOROCRUST_ENABLE_TPL_LOG_CONF)
mark_as_advanced(VOROCRUST_ENABLE_TPL_LOG_BUILD)
mark_as_advanced(VOROCRUST_ENABLE_TPL_LOG_TEST)
mark_as_advanced(VOROCRUST_ENABLE_TPL_LOG_INST)

option(VOROCRUST_TPL_BUILD_SEACAS "Enable TPL SEACAS via ExternalProject" OFF)
option(VOROCRUST_TPL_BUILD_NETCDF "Enable TPL NetCDF via ExternalProject" OFF)
option(VOROCRUST_TPL_BUILD_HDF5   "Enable TPL HDF5 via ExternalProject"   OFF)
option(VOROCRUST_TPL_BUILD_ZLIB   "Enable TPL ZLIB via ExternalProject"   OFF)
mark_as_advanced(FORCE VOROCRUST_TPL_BUILD_SEACAS)
mark_as_advanced(FORCE VOROCRUST_TPL_BUILD_NETCDF)
mark_as_advanced(FORCE VOROCRUST_TPL_BUILD_HDF5)
mark_as_advanced(FORCE VOROCRUST_TPL_BUILD_ZLIB)

option(VOROCRUST_ENABLE_DEVTEST "Enable DEVTEST mode (for developers)" OFF)
mark_as_advanced(VOROCRUST_ENABLE_DEVTEST)

option(VOROCRUST_VERSION_USES_GIT_SHA1 "Set VERSION values by checking the Git SHA1" OFF)
mark_as_advanced(VOROCRUST_VERSION_USES_GIT_SHA1)

option(VOROCRUST_ENABLE_TPL_TESTS  "Enable tests for TPLs (warning: LONG!)"  OFF)
mark_as_advanced(VOROCRUST_ENABLE_TPL_TESTS)

option(VOROCRUST_ENABLE_WARNINGS   "Enable a standard set of warnings flags" OFF)
mark_as_advanced(VOROCRUST_ENABLE_WARNINGS)

# This is a work in progress -- we may not be able to use TuckerMPI -- (EXPERIMENTAL)
# wcmclen - it's looking like TuckerMPI isn't really supported anymore... most of the
#           people associated with it have retired and/or are gone.
option(VOROCRUST_ENABLE_TUCKERMPI "Enable TuckerMPI via ExternalProject"   OFF)
mark_as_advanced(VOROCRUST_ENABLE_TUCKERMPI)

# The Sandia firewall / proxy causes issues from time to time, so
# we keep a copy of our TPLs in the tpls/ dir of the repository.
# Set this to OFF to force VoroCrust to download a copy from the
# web, otherwise we'll use the copy we have in tpls/
option(VOROCRUST_TPL_NO_DOWNLOAD "Use locally stored copies of TPL archives if available. Disable to ALWAYS download" OFF)
mark_as_advanced(VOROCRUST_TPL_NO_DOWNLOAD)


#================================================
# Handle OPTION Dependencies
#================================================
if(VOROCRUST_ENABLE_EXODUS)
    vc_message("VOROCRUST_TPL_BUILD_SEACAS enabled due to VOROCRUST_ENABLE_EXODUS")
    set(VOROCRUST_TPL_BUILD_SEACAS ON CACHE BOOL "Enabled for Exodus" FORCE)
endif()
if(VOROCRUST_TPL_BUILD_SEACAS)
    vc_message("VOROCRUST_TPL_BUILD_NETCDF enabled due to VOROCRUST_ENABLE_SEACAS")
    set(VOROCRUST_TPL_BUILD_NETCDF ON CACHE BOOL "Enabled for SEACAS" FORCE)
endif()
if(VOROCRUST_TPL_BUILD_NETCDF)
    vc_message("VOROCRUST_TPL_BUILD_HDF5   enabled due to VOROCRUST_ENABLE_NETCDF")
    set(VOROCRUST_TPL_BUILD_HDF5   ON CACHE BOOL "Enabled for NetCDF" FORCE)
endif()
if(VOROCRUST_TPL_BUILD_HDF5)
    vc_message("VOROCRUST_TPL_BUILD_ZLIB   enabled due to VOROCRUST_ENABLE_HDF5")
    set(VOROCRUST_TPL_BUILD_ZLIB   ON CACHE BOOL "Enabled for HDF5" FORCE)
endif()

# If user wants to build our own OPENBLAS it will override LAPACK
#if(VOROCRUST_TPL_BUILD_OPENBLAS AND VOROCRUST_TPL_USE_LAPACK)
#    vc_message("${Yellow}WARNING${ColorReset}: VOROCRUST_TPL_BUILD_OPENBLAS overrides VOROCRUST_TPL_USE_LAPACK")
#    vc_message("-- ${Red}Setting${ColorReset} ${Magenta}VOROCRUST_TPL_USE_LAPACK=OFF${ColorReset}")
#    set(VOROCRUST_TPL_USE_LAPACK OFF CACHE BOOL "Disabled by VOROCRUST_TPL_BUILD_OPENBLAS" FORCE)
#endif()

# OpenBLAS requires POSIX and will compile with it (I'm not sure how to turn it off optionally)
# so disabling OpenMP will cause linker errors for missing POSIX symbols on CentOS (*nix) systems.
# Let's throw an error if someone disables OpenMP but enables LAPACK.
if((NOT VOROCRUST_ENABLE_OPENMP) AND (VOROCRUST_TPL_USE_LAPACK OR VOROCRUST_TPL_BUILD_OPENBLAS))
    message(FATAL_ERROR "${TagPREFIXERROR} VOROCRUST_ENABLE_OPENMP can not be disabled if LAPACK or OpenBLAS are enabled")
endif()

# If Kokkos Kernels is enabled, Kokkos should be too
if(VOROCRUST_TPL_ENABLE_KOKKOSKERNELS AND (NOT VOROCRUST_TPL_ENABLE_KOKKOS))
    vc_message("VOROCRUST_TPL_ENABLE_KOKKOSKERNELS requires VOROCRUST_TPL_ENABLE_KOKKOS")
    set(VOROCRUST_TPL_ENABLE_KOKKOS ON CACHE BOOL "Enabled for Kokkos-Kernels" FORCE)
endif()


# If TuckerMPI is enabled then we also need to enable other dependent tpls/flags
# TODO: Fix the logic here to include OpenBLAS conditions
if(VOROCRUST_ENABLE_TUCKERMPI AND (NOT VOROCRUST_ENABLE_MPI))
    vc_message("${Yellow}WARNING${ColorReset}: VOROCRUST_ENABLE_TUCKERMPI requires MPI")
    vc_message("-- ${Red}Setting${ColorReset} ${Magenta}VOROCRUST_ENABLE_MPI=ON${ColorReset}")
    set(VOROCRUST_ENABLE_MPI ON CACHE BOOL "Enabled for TuckerMPI" FORCE)
endif()
if(VOROCRUST_ENABLE_TUCKERMPI AND (NOT VOROCRUST_TPL_USE_LAPACK))
    vc_message("${Yellow}WARNING${ColorReset}: VOROCRUST_ENABLE_TUCKERMPI requires LAPACK")
    vc_message("-- ${Red}Setting${ColorReset} ${Magenta}VOROCRUST_TPL_USE_LAPACK=ON${ColorReset}")
    set(VOROCRUST_TPL_USE_LAPACK ON CACHE BOOL "Enabled for TuckerMPI" FORCE)
endif()


#========================================
# Initialize GLOBAL scope vars
#========================================
set(CMAKE_COLOR_CMAKE ON CACHE BOOL "Enable colorized CMake on non-windows platforms")
mark_as_advanced(CMAKE_COLOR_CMAKE)

#========================================
# Enable Testing Globally
#========================================
enable_testing()


#========================================
# Import configuration files
#========================================

include(vc_windows)
include(vc_unix)
include(vc_mpi)
include(vc_openmp)
include(vc_kokkos)
include(vc_kokkoskernels)

vc_message_var(VOROCRUST_TPL_BUILD_OPENBLAS)
vc_message_var(VOROCRUST_TPL_USE_LAPACK)

include(vc_tpl_lapack)

include(vc_tpl_tuckermpi)

if(VOROCRUST_VERSION_USES_GIT_SHA1)
    include(vc_git_tools)
else()
    vc_message("Skip setting version string via git sha1 lookup")
    # Set Default values for GIT_COMMIT_* vars
    set(GIT_COMMIT_SHA1 "UNKNOWN")
    set(GIT_COMMIT_DATE "UNKNOWN")
    vc_message_var(GIT_COMMIT_SHA1)
    vc_message_var(GIT_COMMIT_DATE)
endif()

# SEACAS (Exodus) TPLS
include(vc_tpl_zlib)
include(vc_tpl_hdf5)
include(vc_tpl_netcdf)
include(vc_tpl_seacas)
include(vc_exodus)


# Include the src/ directory
add_subdirectory(src)

# Add src/ directory to the -I<INCLUDE> paths for compilation
include_directories(src)

# Optionally prints out some useful CMake vars
# if VOROCRUST_VERBOSE_CMAKE is enabled.
print_cmake_vars()

# Verify that the compiler is compatible with VoroCrust.
vorocrust_check_compiler_version()

# Include the Global tests subdirectory
add_subdirectory(tests)



#========================================
# Installation
#========================================
vc_print_banner("BEGIN: Configure Installation & Packaging Properties")
if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0" )

    get_property(VOROCRUST_BUILD_TARGETS GLOBAL PROPERTY "VOROCRUST_BUILD_TARGETS")

    set(VOROCRUST_COPYRIGHT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/COPYRIGHT.md")
    set(VOROCRUST_LICENSE_FILE   "")
    set(VOROCRUST_README_FILE    "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

    install(TARGETS ${VOROCRUST_BUILD_TARGETS}
            RUNTIME DESTINATION "bin"
            LIBRARY DESTINATION "lib"
            ARCHIVE DESTINATION "lib"
            COMPONENT applications
            )

    # Add examples to the distribution package
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/tests/datasets/README.md
            DESTINATION share/${PROJECT_NAME}/examples
            )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/datasets/Bird
                      ${CMAKE_CURRENT_SOURCE_DIR}/tests/datasets/Bracket
                      ${CMAKE_CURRENT_SOURCE_DIR}/tests/datasets/Non_Manifold_Box
                      ${CMAKE_CURRENT_SOURCE_DIR}/tests/datasets/OneBox
                      ${CMAKE_CURRENT_SOURCE_DIR}/tests/datasets/Torus
            DESTINATION share/${PROJECT_NAME}/examples
            )
    install(FILES ${VOROCRUST_COPYRIGHT_FILE}
                  ${VOROCRUST_README_FILE}
            DESTINATION share/${PROJECT_NAME}
            )

    # Add project files
    install(FILES ${VOROCRUST_COPYRIGHT_FILE}
                  ${VOROCRUST_README_FILE}
            DESTINATION share/${PROJECT_NAME}
            )

    #========================================
    # BEGIN Packaging
    #========================================
    # See: https://gitlab.kitware.com/cmake/community/-/wikis/doc/cpack/Packaging-With-CPack#using-cmake-variables-to-configure-cpack
    #      for information useful for packaging on windows.
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "VoroCrust mesh generator")
    set(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
    #set(CPACK_RESOURCE_FILE_LICENSE ${VOROCRUST_LICENSE_FILE})
    set(CPACK_RESOURCE_FILE_README  ${VOROCRUST_README_FILE})
    set(CPACK_VERBATIM_VARIABLES true)

    # Default generator is a ZIP file.
    # This can be changed via -G <GENERATOR> or -DCPACK_GENERATOR="ZIP;TGZ"
    set(CPACK_GENERATOR "ZIP")

    # Control filename for output package
    # See: https://cmake.org/cmake/help/v3.8/module/CPack.html#variable:CPACK_PACKAGE_FILE_NAME
    # Default configuration is "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}"
    string(TIMESTAMP TODAY "%Y%m%d")

    set(CPACK_PACKAGE_VERSION ${TODAY})

    set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")

    if(APPLE)
        cmake_host_system_information(RESULT OS_RELEASE QUERY OS_RELEASE)
        #vc_message("OS_RELEASE = ${OS_RELEASE}")
        set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${OS_RELEASE}-${CMAKE_SYSTEM_PROCESSOR}")
    endif(APPLE)

    vc_message_var(CMAKE_SYSTEM_NAME)

    if(CMAKE_SYSTEM_NAME MATCHES "Linux")
        vc_message_var(CMAKE_SYSTEM_NAME)
        get_linux_lsb_release_information()

        vc_message_var(LSB_RELEASE_ID_SHORT)
        vc_message_var(LSB_RELEASE_VERSION_SHORT)
        set(CPACK_SYSTEM_NAME "${LSB_RELEASE_ID_SHORT}-${LSB_RELEASE_VERSION_SHORT}-${CMAKE_SYSTEM_PROCESSOR}")
    endif()

    vc_message_var(CMAKE_SYSTEM_NAME)
    vc_message_var(CMAKE_SYSTEM_PROCESSOR)
    vc_message_var(CPACK_SYSTEM_NAME)

    set(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/resources/icons/VoroCrust.png")

    if (UNIX)
        set(CPACK_STRIP_FILES ${PROJECT_NAME})
        set(CPACK_SOURCE_STRIP_FILES "")
    endif()

    set(CPACK_PACKAGE_EXECUTABLES ${PROJECT_NAME} "VoroCrust Executable")
    include(CPack)
    #========================================
    # END Packaging
    #========================================

else()

    vc_message("Installation and Packaging disabled. To enable please upgrade CMake to 3.16 or greater.")

endif()

vc_print_banner("COMPLETE: Installation & Packaging Properties")

