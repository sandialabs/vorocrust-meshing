include_guard(GLOBAL)
# include(vc_functions)

#
# Git tools
#
message("BEGIN: ${CMAKE_CURRENT_LIST_FILE}")

# Function: FIND_GIT_ROOT
# Locates the git root of the repository by searching up the directory
# tree starting at the location of CMAKE_SOURCE_DIR until it
# either finds a .git directory or can't go anymore. If not found, then
# returns UNKNOWN.
#
# Usage:
#   find_git_root(GIT_ROOT_DIR)
#
function(FIND_GIT_ROOT git_root_dir_var)
    set(GIT_ROOT_DIR "${CMAKE_SOURCE_DIR}")
    # Maybe also CMAKE_SOURCE_DIR / CMAKE_CURRENT_SOURCE_DIR

    while( NOT EXISTS "${GIT_ROOT_DIR}/.git")
        get_filename_component(NEW_GIT_ROOT_DIR ${GIT_ROOT_DIR} DIRECTORY)
        if( GIT_ROOT_DIR STREQUAL NEW_GIT_ROOT_DIR )
            set(${git_root_dir_var} "UNKNOWN" PARENT_SCOPE)
            return()
        endif()
        set(GIT_ROOT_DIR "${NEW_GIT_ROOT_DIR}")
    endwhile()
    set(${git_root_dir_var} "${GIT_ROOT_DIR}" PARENT_SCOPE)

    if(VOROCRUST_VERBOSE_CMAKE)
        message(STATUS "${TagVERBOSE} GIT_ROOT_DIR: ${GIT_ROOT_DIR}")
    endif()

endfunction()



# Function: GET_GIT_HEAD_REF
# Gets the reference from the <git root>/.git/HEAD file.
# This is formatted as "ref: refs/heads/<current branch name>"
#
# Usage:
#   get_git_head_ref( git_ref_file git_root_dir)
#
function(GET_GIT_HEAD_REF git_head_ref_var git_root_dir_var)
    # HEAD will be formatted as "ref: refs/heads/<current-branch-name>"
    # - strip "ref: " from the string to get the path component
    set(GIT_REF_FILE "UNKNOWN")
    if(EXISTS "${git_root_dir_var}/.git/HEAD")
        file(STRINGS "${git_root_dir_var}/.git/HEAD" GIT_HEAD LIMIT_COUNT 1)
        string(REPLACE "ref: " "" GIT_REF_FILE ${GIT_HEAD})
    endif()
    set(${git_head_ref_var} "${GIT_REF_FILE}" PARENT_SCOPE)

    if(VOROCRUST_VERBOSE_CMAKE)
        message(STATUS "${TagVERBOSE} GIT_REF_FILE: ${GIT_REF_FILE}")
    endif()

endfunction()



# Function : GET_GIT_SHA1
# Get the current SHA1 of the repository (if a git repo) without
# using a git command to do so. This can be useful if the source
# code is downloaded as a .zip file and the building machien does
# not have git installed.
#
# Usage:
#     GET_GIT_SHA1( git_sha1_variable )
#
function(GET_GIT_SHA1 git_sha1_variable)
    set(GIT_SHA1 "UNKNOWN")
    #set(GIT_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    set(GIT_ROOT_DIR "${CMAKE_SOURCE_DIR}")

    find_git_root(GIT_ROOT_DIR)

    if(NOT GIT_ROOT_DIR STREQUAL "UNKNOWN")
        get_git_head_ref(GIT_REF_FILE "${GIT_ROOT_DIR}")

        if(EXISTS "${GIT_ROOT_DIR}/.git/${GIT_REF_FILE}")
            file(STRINGS "${GIT_ROOT_DIR}/.git/${GIT_REF_FILE}" GIT_SHA1 LIMIT_COUNT 1)
        endif()
    endif()

    set(${git_sha1_variable} "${GIT_SHA1}" PARENT_SCOPE)

    if(VOROCRUST_VERBOSE_CMAKE)
        message(STATUS "${TagVERBOSE} GIT_SHA1: ${GIT_SHA1}")
    endif()

endfunction()



# Function : GET_LAST_LINE_FROM_TEXTFILE
# Extract the last line from a text file
#
# Usage:
#     GET_LAST_LINE_FROM_TEXTFILE(last_line "${filename}")
#
function(GET_LAST_LINE_FROM_TEXTFILE output filename)
    set(_output "UNKNOWN")
    if(EXISTS "${filename}")
        file(STRINGS "${filename}" data)
        list(LENGTH data num_lines)
        math(EXPR last_line_idx "${num_lines} - 1")
        list(GET data ${last_line_idx} _output)
    else()
        message(WARNING "get_last_line_from_textfile(): file not found: ${filename}")
    endif()
    set(${output} "${_output}" PARENT_SCOPE)
endfunction()



# Function: GET_GIT_DATETIME
# Gets the date-time stamp of the last commit in a git repo
#
# Usage:
#    GET_GIT_DATETIME(git_datetime)
#
# Todo: This won't work if in headless mode in git. We should make it more robust for when we find git but the
#       refs file doesn't exist.
function(GET_GIT_DATETIME git_datetime_var)
    set(timestamp_unix "UNKNOWN")
    #set(GIT_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    set(GIT_ROOT_DIR "${CMAKE_SOURCE_DIR}")

    find_git_root(GIT_ROOT_DIR)

    if(NOT GIT_ROOT_DIR STREQUAL "UNKNOWN")
        get_git_head_ref(GIT_REF_FILE "${GIT_ROOT_DIR}")

        if(EXISTS "${GIT_ROOT_DIR}/.git/logs/${GIT_REF_FILE}")
            # Get the LAST line from file
            get_last_line_from_textfile(last_log_entry "${GIT_ROOT_DIR}/.git/logs/${GIT_REF_FILE}")

            # Use a REGEX to extract the unix time value out from the Git commit log entry.
            # Note: https://regex101.com/ is a really useful REGEX tester.
            #string(REGEX MATCH "(([0-9]+)( -[0-9]+))[ \t]+((commit.*)|(pull.*)|(branch.*)|(rebase.*)):.*$" args "${last_log_entry}")
            string(REGEX MATCH "(([0-9]+)( [-+][0-9]+))[ \t]+([a-z]+.*):.*$" args "${last_log_entry}")
            set(timestamp_unix "${CMAKE_MATCH_2}")  # 1=unix_time+offset, 2=unix_time, 3=offset, 4=commit (merge):
        endif()
    endif()

    # Error out if the output string is empty
    if("${timestamp_unix}" STREQUAL "")
        message(FATAL_ERROR " get_git_datetime() was unable to generate a timestamp string.")
    endif()

    set(${git_datetime_var} "${timestamp_unix}" PARENT_SCOPE)

    if(VOROCRUST_VERBOSE_CMAKE)
        message(STATUS "${TagVERBOSE} timestamp_unix: ${timestamp_unix}")
    endif()


endfunction()


#==============================================================================
# Execute commands to get the SHA1 and latest COMMIT_DATE
#==============================================================================

message("[VC] -- Looking up git sha1")

# Get the SHA1 of the current branch without using git.
get_git_sha1(GIT_COMMIT_SHA1)
get_git_datetime(GIT_COMMIT_DATE)

message("[VC] -- GIT_COMMIT_SHA1 = ${GIT_COMMIT_SHA1}")
message("[VC] -- GIT_COMMIT_DATE = ${GIT_COMMIT_DATE}")

message ("COMPLETE: ${CMAKE_CURRENT_LIST_FILE}")
