include_guard(GLOBAL)
#
# External Project for LAPACK
#
include(vc_functions)

vc_print_banner("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

# LAPACK State Variables
set(VOROCRUST_USE_SYSTEM_LAPACK OFF)
set(VOROCRUST_USE_BUILT_LAPACK  OFF)
set(VOROCRUST_HAVE_LAPACK_VAR   "NO_LAPACK")

vc_message_var(VOROCRUST_TPL_USE_LAPACK)

if(VOROCRUST_TPL_USE_LAPACK)
    if(NOT VOROCRUST_TPL_BUILD_OPENBLAS)
        set(BLA_STATIC ON)
        find_package(BLAS REQUIRED)
        vc_message_var(BLAS_FOUND)
        vc_message_var(BLAS_LINKER_FLAGS)
        vc_message_var(BLAS_LIBRARIES)

        find_package(LAPACK REQUIRED)
        vc_message_var(LAPACK_FOUND)
        vc_message_var(LAPACK_LINKER_FLAGS)
        vc_message_var(LAPACK_LIBRARIES)
        vc_message("")

        if(LAPACK_FOUND)
            vc_message("-- SYSTEM LAPACK ${Green}FOUND${ColorReset}")
            set(VOROCRUST_USE_SYSTEM_LAPACK ON)
        else()
            vc_message("-- SYSTEM LAPACK ${Red}NOT FOUND${ColorReset}")
            set(VOROCRUST_USE_BUILT_LAPACK ON)
        endif()
    endif()
endif()


if(VOROCRUST_USE_BUILT_LAPACK OR VOROCRUST_TPL_BUILD_OPENBLAS)
    include(vc_tpl_openblas)
    set(VOROCRUST_USE_BUILT_LAPACK ON)
endif()

vc_message_var(VOROCRUST_USE_SYSTEM_LAPACK)
vc_message_var(VOROCRUST_USE_BUILT_LAPACK)

macro(VOROCRUST_TPLADD_LAPACK TARGET)
    vc_message("vorocrust_tpladd_lapack(${TARGET})")
    if(VOROCRUST_USE_BUILT_LAPACK)
        vc_message("Adding Built BLAS/LAPACK (OpenBLAS) to ${TARGET}")
        add_dependencies(${TARGET} openblas)
        target_include_DIRECTORIES(${TARGET} PUBLIC ${TPL_OPENBLAS_INCLUDE_DIR})
        target_link_libraries(${TARGET} PUBLIC ${TPL_OPENBLAS_LIBRARIES})
    elseif(VOROCRUST_USE_SYSTEM_LAPACK)
        vc_message("Adding System BLAS/LAPACK to ${TARGET}")
        target_link_options(${TARGET}   PUBLIC ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
        target_link_libraries(${TARGET} PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()
endmacro()


vc_print_banner("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")

