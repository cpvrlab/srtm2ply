
include(ExternalProject)
set_property(DIRECTORY ${CMAKE_CURRENT_LIST_DIR} PROPERTY EP_PREFIX)

#Eigen
file(DOWNLOAD http://bitbucket.org/eigen/eigen/get/3.2.5.tar.gz
              ${CMAKE_CURRENT_SOURCE_DIR}/eigen.tar.gz
     EXPECTED_MD5 8cc513ac6ec687117acadddfcacf551b
    )
set(EIGEN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
set(EIGEN_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen)
set_directory_properties(PROPERTIES
                         EIGEN_SOURCE_DIR    ${EIGEN_SOURCE_DIR}
                         EIGEN_BINARY_DIR    ${EIGEN_BINARY_DIR}
                         EIGEN_INCLUDE_DIR   ${EIGEN_SOURCE_DIR}
                        )
ExternalProject_Add(
   external_eigen
   URL file://${CMAKE_CURRENT_SOURCE_DIR}/eigen.tar.gz
   CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
        -DEIGEN_SPLIT_LARGE_TESTS:BOOL=OFF
        -DEIGEN_LEAVE_TEST_IN_ALL_TARGET:BOOL=OFF
        -DEIGEN_BUILD_PKGCONFIG:BOOL=OFF
        -DEIGEN_BUILD_BTL:BOOL=OFF
   SOURCE_DIR ${EIGEN_SOURCE_DIR}
   BINARY_DIR ${EIGEN_BINARY_DIR}
)
#add_subdirectory(eigen)

#Clipper
# file(DOWNLOAD http://netassist.dl.sourceforge.net/project/polyclipping/clipper_ver6.2.1.zip
#               ${CMAKE_CURRENT_SOURCE_DIR}/clipper.zip
#      EXPECTED_MD5 040821e66ec529f3d78f8ff7c4e256b2
#     )
# set(CLIPPER_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/clipper)
# set(CLIPPER_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/clipper)
# set_directory_properties(PROPERTIES
#                          CLIPPER_SOURCE_DIR     ${CLIPPER_SOURCE_DIR}/cpp
#                          CLIPPER_BINARY_DIR     ${CLIPPER_BINARY_DIR}
#                          CLIPPER_INCLUDE_DIR    ${CLIPPER_SOURCE_DIR}/cpp
#                          CLIPPER_LIBRARY_DIR    ${CLIPPER_BINARY_DIR}
#                          CLIPPER_LIBRARY_FILES  polyclipping
#                          )
# ExternalProject_Add(
#    external_clipper
#    URL file://${CMAKE_CURRENT_SOURCE_DIR}/clipper.zip
#    SOURCE_DIR ${CLIPPER_SOURCE_DIR}
#    BINARY_DIR ${CLIPPER_BINARY_DIR}
#    CONFIGURE_COMMAND ""
#    BUILD_COMMAND ${CMAKE_COMMAND}
#         -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
#         -DCMAKE_BINARY_DIR:STRING=<BINARY_DIR>
#         --build <SOURCE_DIR>/cpp
# )
# #add_subdirectory(clipper/cpp)

#Docopt
if (NOT MSVC OR NOT MSVC_VERSION LESS 1900)

file(DOWNLOAD https://codeload.github.com/docopt/docopt.cpp/zip/master
              ${CMAKE_CURRENT_SOURCE_DIR}/docopt.zip
     #Docopt.cpp only provides a master branch, so there are no stable releases.
     #EXPECTED_MD5 e3bf762664425fea1da2d0715db9b7e3
    )
set(DOCOPT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/docopt)
set(DOCOPT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/docopt)
set_directory_properties(PROPERTIES
                         DOCOPT_SOURCE_DIR     ${DOCOPT_SOURCE_DIR}
                         DOCOPT_BINARY_DIR     ${DOCOPT_BINARY_DIR}
                         DOCOPT_INCLUDE_DIR    ${DOCOPT_SOURCE_DIR}
                         DOCOPT_LIBRARY_DIR    ${DOCOPT_BINARY_DIR}
                         DOCOPT_LIBRARY_FILES  docopt
                         )
if (MSVC)
  set_directory_properties(PROPERTIES
                           DOCOPT_LIBRARY_FILES  "docopt;docopt_s"
                           )
endif()
ExternalProject_Add(
   external_docopt
   URL file://${CMAKE_CURRENT_SOURCE_DIR}/docopt.zip
   SOURCE_DIR ${DOCOPT_SOURCE_DIR}
   BINARY_DIR ${DOCOPT_BINARY_DIR}
   CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
)
#add_subdirectory(docopt)

endif()

#Geographiclib
file(DOWNLOAD http://heanet.dl.sourceforge.net/project/geographiclib/distrib/GeographicLib-1.44.tar.gz
              ${CMAKE_CURRENT_SOURCE_DIR}/geographiclib.tar.gz
     EXPECTED_MD5 212b1bf07b9b672aeae03f845f01e76b
    )
set(GEOGRAPHICLIB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/geographiclib)
set(GEOGRAPHICLIB_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/geographiclib)
set_directory_properties(PROPERTIES
                         GEOGRAPHICLIB_SOURCE_DIR     ${GEOGRAPHICLIB_SOURCE_DIR}/src
                         GEOGRAPHICLIB_BINARY_DIR     ${GEOGRAPHICLIB_BINARY_DIR}
                         GEOGRAPHICLIB_INCLUDE_DIR    ${GEOGRAPHICLIB_SOURCE_DIR}/include
                         GEOGRAPHICLIB_LIBRARY_DIR    ${GEOGRAPHICLIB_BINARY_DIR}/src
                         GEOGRAPHICLIB_LIBRARY_FILES  Geographic
                         )
if (MSVC)
  set_directory_properties(PROPERTIES
                           GEOGRAPHICLIB_LIBRARY_DIR    ${GEOGRAPHICLIB_BINARY_DIR}/lib
                           GEOGRAPHICLIB_LIBRARY_FILES  "Geographic;Geographic_d"
                           )
endif()
ExternalProject_Add(
   external_geographiclib
   URL file://${CMAKE_CURRENT_SOURCE_DIR}/geographiclib.tar.gz
   SOURCE_DIR ${GEOGRAPHICLIB_SOURCE_DIR}
   BINARY_DIR ${GEOGRAPHICLIB_BINARY_DIR}
   CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
)
#add_subdirectory(geographiclib)

#Dirent on Windows
if (MSVC)
  file(DOWNLOAD https://github.com/tronkko/dirent/archive/master.zip
                ${CMAKE_CURRENT_SOURCE_DIR}/dirent.zip
       EXPECTED_MD5 7727093d6e8b8ac50c4789810ceb9dd3
      )
  set(DIRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dirent)
  set(DIRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/dirent)
  set_directory_properties(PROPERTIES
                           DIRENT_SOURCE_DIR     ${DIRENT_SOURCE_DIR}
                           DIRENT_INCLUDE_DIR    ${DIRENT_SOURCE_DIR}/include
                           )
  ExternalProject_Add(
     external_dirent
     URL file://${CMAKE_CURRENT_SOURCE_DIR}/dirent.zip
     SOURCE_DIR ${DIRENT_SOURCE_DIR}
     BINARY_DIR ${DIRENT_BINARY_DIR}
     CONFIGURE_COMMAND ""
     UPDATE_COMMAND    ""
     BUILD_COMMAND     ""
     PATCH_COMMAND     ""
     INSTALL_COMMAND   ""
  )
  #add_subdirectory(dirent)
else()
  add_custom_target(external_dirent)
endif(MSVC)