#cmake_minimum_required (VERSION 2.6)
cmake_minimum_required(VERSION 3.3)
set(CMAKE_VERBOSE_MAKEFILE ON)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

function(FIND_LIBRARIES)
    include(CMakeParseArguments)
    CMAKE_PARSE_ARGUMENTS(PARSED_ARGS "" "VARNAME" "FILES;DIRECTORIES" ${ARGN})
    set(libs "")
    foreach(lib ${PARSED_ARGS_FILES})
        if (MSVC)
            find_library(libfile       ${lib} HINTS ${PARSED_ARGS_DIRECTORIES}/Release NO_DEFAULT_PATH)
            find_library(libfile_debug ${lib} HINTS ${PARSED_ARGS_DIRECTORIES}/Debug   NO_DEFAULT_PATH)

            message(Looking for ${lib})
            message(Found: ${libfile})
            message(${PARSED_ARGS_DIRECTORIES}/Release)

            if(NOT libfile_debug MATCHES ".*-NOTFOUND")
                list(APPEND libs debug)
                list(APPEND libs ${libfile_debug})
            endif()
            if(NOT libfile MATCHES ".*-NOTFOUND")
                list(APPEND libs optimized)
                list(APPEND libs ${libfile})
            endif()
            unset(libfile CACHE)
            unset(libfile_debug CACHE)
        else()
            find_library(libfile ${lib} HINTS ${PARSED_ARGS_DIRECTORIES} NO_DEFAULT_PATH)
            if(NOT libfile MATCHES ".*-NOTFOUND")
                list(APPEND libs ${libfile})
            endif()
            unset(libfile CACHE)
        endif()
    endforeach()
    set("${PARSED_ARGS_VARNAME}" "${libs}" PARENT_SCOPE)
    unset(libs)
endfunction()

#===========================================#
#             PROJECT SETTINGS              #
#===========================================#

project (srtm2ply)

if(MSVC)
    # set(CMAKE_CXX_FLAGS "/openmp")
elseif(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}    -fopenmp")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

set(HEADERS
        ./src/Architecture.h
        ./src/Geometry.h
        ./src/Mesh.h
        ./src/Parsing.h
        ./src/SRTM.h
        ./src/Tensor.h
        ./src/Utilities.h
    )

set(SOURCES
        ./src/main.cpp

        ./src/Geometry.cpp
        ./src/Parsing.cpp
        ./src/Mesh.cpp
        ./src/Utilities.cpp
    )

add_executable(srtm2ply ${HEADERS} ${SOURCES})

#===========================================#
#       EXTERNAL PROJECT SETTINGS           #
#===========================================#

add_subdirectory(external)
get_directory_property(EIGEN_INCLUDE_DIR             DIRECTORY external EIGEN_INCLUDE_DIR           )

get_directory_property(CLIPPER_INCLUDE_DIR           DIRECTORY external CLIPPER_INCLUDE_DIR         )
get_directory_property(CLIPPER_SOURCE_DIR            DIRECTORY external CLIPPER_SOURCE_DIR          )
get_directory_property(CLIPPER_BINARY_DIR            DIRECTORY external CLIPPER_BINARY_DIR          )
get_directory_property(CLIPPER_LIBRARY_DIR           DIRECTORY external CLIPPER_LIBRARY_DIR         )
get_directory_property(CLIPPER_LIBRARY_FILES         DIRECTORY external CLIPPER_LIBRARY_FILES       )

get_directory_property(DOCOPT_INCLUDE_DIR            DIRECTORY external DOCOPT_INCLUDE_DIR          )
get_directory_property(DOCOPT_SOURCE_DIR             DIRECTORY external DOCOPT_SOURCE_DIR           )
get_directory_property(DOCOPT_BINARY_DIR             DIRECTORY external DOCOPT_BINARY_DIR           )
get_directory_property(DOCOPT_LIBRARY_DIR            DIRECTORY external DOCOPT_LIBRARY_DIR          )
get_directory_property(DOCOPT_LIBRARY_FILES          DIRECTORY external DOCOPT_LIBRARY_FILES        )

get_directory_property(GEOGRAPHICLIB_INCLUDE_DIR     DIRECTORY external GEOGRAPHICLIB_INCLUDE_DIR   )
get_directory_property(GEOGRAPHICLIB_SOURCE_DIR      DIRECTORY external GEOGRAPHICLIB_SOURCE_DIR    )
get_directory_property(GEOGRAPHICLIB_BINARY_DIR      DIRECTORY external GEOGRAPHICLIB_BINARY_DIR    )
get_directory_property(GEOGRAPHICLIB_LIBRARY_DIR     DIRECTORY external GEOGRAPHICLIB_LIBRARY_DIR   )
get_directory_property(GEOGRAPHICLIB_LIBRARY_FILES   DIRECTORY external GEOGRAPHICLIB_LIBRARY_FILES )

get_directory_property(DIRENT_INCLUDE_DIR            DIRECTORY external DIRENT_INCLUDE_DIR          )
get_directory_property(DIRENT_SOURCE_DIR             DIRECTORY external DIRENT_SOURCE_DIR           )

FIND_LIBRARIES(VARNAME libClipper        FILES ${CLIPPER_LIBRARY_FILES}        DIRECTORIES ${CLIPPER_LIBRARY_DIR}       )
FIND_LIBRARIES(VARNAME libDocopt         FILES ${DOCOPT_LIBRARY_FILES}         DIRECTORIES ${DOCOPT_LIBRARY_DIR}        )
FIND_LIBRARIES(VARNAME libGeographiclib  FILES ${GEOGRAPHICLIB_LIBRARY_FILES}  DIRECTORIES ${GEOGRAPHICLIB_LIBRARY_DIR} )

message(libClipper:        ${libClipper}       )
message(libDocopt:         ${libDocopt}        )
message(libGeographiclib:  ${libGeographiclib} )

#===========================================#
#            COMPILER SETTINGS              #
#===========================================#

# Bump up warning levels appropriately for clang, gcc & msvc
# Also set debug/optimization flags depending on the build type. IDE users choose this when
# selecting the build mode in their IDE
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -std=c++11")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -g")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O2")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
    if (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
endif()

if (#NOT libClipper        OR
    NOT libDocopt         OR
    NOT libGeographiclib  #OR
    )
    message("One or multiple external libraries are not yet compiled.")
    message("CMake will automatically be re-run during your next build.")
    #rerun cmake in initial build
    #will update cmakecache/project files on first build
    #so you may have to reload project after first build
    add_custom_target(Rescan ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}
                      DEPENDS
                        external_eigen
#                       external_clipper
                        external_docopt
                        external_geographiclib
                        external_dirent
                      )
else()
    #Rescan becomes a dummy target after first build
    #this prevents cmake from rebuilding cache/projects on subsequent builds
    add_custom_target(Rescan)

    set(ADDITIONAL_LIBRARIES)
    set(ADDITIONAL_INCLUDE_DIRS)

    message(CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME})
    message(CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID})

    if     ("${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
        #Handle Windows systems

    elseif ("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
        #Handle generic desktop Linux systems

        set(ADDITIONAL_INCLUDE_DIRS "/usr/include")
        set(ADDITIONAL_LIBRARIES)
    endif()

    include_directories(${ADDITIONAL_INCLUDE_DIRS}
                        ${EIGEN_INCLUDE_DIR}
                        ${CLIPPER_INCLUDE_DIR}
                        ${GEOGRAPHICLIB_INCLUDE_DIR}
                        ${DOCOPT_INCLUDE_DIR}
                        ${DIRENT_INCLUDE_DIR}
                       )

    target_link_libraries (srtm2ply
                           ${libClipper}
                           ${libDocopt}
                           ${libGeographiclib}
                           ${libBoost}

                           ${ADDITIONAL_LIBRARIES}
                           )

    if (MSVC)
        target_link_libraries(srtm2ply shlwapi.lib)
    endif()
endif()

add_dependencies(srtm2ply Rescan)
