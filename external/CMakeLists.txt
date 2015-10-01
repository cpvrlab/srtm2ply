
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
file(DOWNLOAD http://netassist.dl.sourceforge.net/project/polyclipping/clipper_ver6.2.1.zip
              ${CMAKE_CURRENT_SOURCE_DIR}/clipper.zip
     EXPECTED_MD5 040821e66ec529f3d78f8ff7c4e256b2
    )
set(CLIPPER_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/clipper)
set(CLIPPER_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/clipper)
set_directory_properties(PROPERTIES
                         CLIPPER_SOURCE_DIR     ${CLIPPER_SOURCE_DIR}/cpp
                         CLIPPER_BINARY_DIR     ${CLIPPER_BINARY_DIR}
                         CLIPPER_INCLUDE_DIR    ${CLIPPER_SOURCE_DIR}/cpp
                         CLIPPER_LIBRARY_DIR    ${CLIPPER_BINARY_DIR}
                         CLIPPER_LIBRARY_FILES  polyclipping
                         )
ExternalProject_Add(
   external_clipper
   URL file://${CMAKE_CURRENT_SOURCE_DIR}/clipper.zip
   SOURCE_DIR ${CLIPPER_SOURCE_DIR}
   BINARY_DIR ${CLIPPER_BINARY_DIR}
   CONFIGURE_COMMAND ""
   BUILD_COMMAND ${CMAKE_COMMAND}
        -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
        -DCMAKE_BINARY_DIR:STRING=<BINARY_DIR>
        --build <SOURCE_DIR>/cpp
)
#add_subdirectory(clipper/cpp)

#Docopt
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
ExternalProject_Add(
   external_docopt
   URL file://${CMAKE_CURRENT_SOURCE_DIR}/docopt.zip
   SOURCE_DIR ${DOCOPT_SOURCE_DIR}
   BINARY_DIR ${DOCOPT_BINARY_DIR}
   CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
)
#add_subdirectory(docopt)

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
       EXPECTED_MD5 212b1bf07b9b672aeae03f845f01e76b
      )
  set(DIRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dirent)
  set(DIRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/dirent)
  set_directory_properties(PROPERTIES
                           DIRENT_SOURCE_DIR     ${DIRENT_SOURCE_DIR}
                           DIRENT_INCLUDE_DIR    ${DIRENT_SOURCE_DIR}/include
                           )
  ExternalProject_Add(
     external_dirent
     URL file://${CMAKE_CURRENT_SOURCE_DIR}/dirent.tar.gz
     SOURCE_DIR ${DIRENT_SOURCE_DIR}
     BINARY_DIR ${DIRENT_BINARY_DIR}
     UPDATE_COMMAND ""
     BUILD_COMMAND  ""
     INSTALL_COMMAND ""
  )
  #add_subdirectory(dirent)
else()
  add_custom_target(external_dirent)
endif(MSVC)

# #OGLplus
# file(DOWNLOAD https://codeload.github.com/matus-chochlik/oglplus/zip/master
#               ${CMAKE_CURRENT_SOURCE_DIR}/oglplus.zip
#      #EXPECTED_MD5 8cc513ac6ec687117acadddfcacf551b
#     )
# set(OGLPLUS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/oglplus)
# set(OGLPLUS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/oglplus)
# set_directory_properties(PROPERTIES
#                          OGLPLUS_SOURCE_DIR    ${OGLPLUS_SOURCE_DIR}
#                          OGLPLUS_BINARY_DIR    ${OGLPLUS_BINARY_DIR}
#                          OGLPLUS_INCLUDE_DIR   ${OGLPLUS_SOURCE_DIR}
#                         )
# ExternalProject_Add(
#    external_oglplus
#    URL file://${CMAKE_CURRENT_SOURCE_DIR}/oglplus.zip
#    CMAKE_ARGS
#         -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
#         -DOGLPLUS_NO_EXAMPLES:BOOL=ON
#         -DOGLPLUS_NO_DOCS:BOOL=ON
#    SOURCE_DIR ${OGLPLUS_SOURCE_DIR}
#    BINARY_DIR ${OGLPLUS_BINARY_DIR}
# )

# # #GLM
# # file(DOWNLOAD https://codeload.github.com/g-truc/glm/zip/0.9.6
# #               ${CMAKE_CURRENT_SOURCE_DIR}/glm.zip
# #      #EXPECTED_MD5 19caabdac1274824e9e9c667e22c4e94
# #     )
# # set(GLM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glm)
# # set(GLM_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/glm)
# # set_directory_properties(PROPERTIES
# #                          GLM_SOURCE_DIR    ${GLM_SOURCE_DIR}
# #                          GLM_BINARY_DIR    ${GLM_BINARY_DIR}
# #                          GLM_INCLUDE_DIR   ${GLM_SOURCE_DIR}
# #                         )
# # ExternalProject_Add(
# #    external_glm
# #    URL file://${CMAKE_CURRENT_SOURCE_DIR}/glm.zip
# #    CMAKE_ARGS
# #         -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
# #    SOURCE_DIR ${GLM_SOURCE_DIR}
# #    BINARY_DIR ${GLM_BINARY_DIR}
# # )

# #GLFW
# file(DOWNLOAD https://codeload.github.com/glfw/glfw/zip/3.1.1
#               ${CMAKE_CURRENT_SOURCE_DIR}/glfw.zip
#      #EXPECTED_MD5 19caabdac1274824e9e9c667e22c4e94
#     )
# set(GLFW_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glfw)
# set(GLFW_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/glfw)
# set_directory_properties(PROPERTIES
#                          GLFW_SOURCE_DIR    ${GLFW_SOURCE_DIR}
#                          GLFW_BINARY_DIR    ${GLFW_BINARY_DIR}
#                          GLFW_INCLUDE_DIR   ${GLFW_SOURCE_DIR}/include
#                          GLFW_LIBRARY_DIR   ${GLFW_BINARY_DIR}/src
#                          GLFW_LIBRARY_FILES glfw3
#                         )
# ExternalProject_Add(
#    external_glfw
#    URL file://${CMAKE_CURRENT_SOURCE_DIR}/glfw.zip
#    CMAKE_ARGS
#         -DCMAKE_BUILD_TYPE:STRING=Release
#         -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
#         -DBUILD_SHARED_LIBS:BOOL=OFF
#         -DGLFW_BUILD_EXAMPLES:BOOL=OFF
#         -DGLFW_BUILD_TESTS:BOOL=OFF
#         -DGLFW_BUILD_DOCS:BOOL=OFF
#         -DGLFW_DOCUMENT_INTERNALS:BOOL=OFF
#         -DGLFW_INSTALL:BOOL=ON
#    SOURCE_DIR ${GLFW_SOURCE_DIR}
#    BINARY_DIR ${GLFW_BINARY_DIR}
# )

# # #SFML
# # file(DOWNLOAD http://mirror0.sfml-dev.org/files/SFML-2.3.1-sources.zip
# #               ${CMAKE_CURRENT_SOURCE_DIR}/sfml.zip
# #      EXPECTED_MD5 19caabdac1274824e9e9c667e22c4e94
# #     )
# # set(SFML_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sfml)
# # set(SFML_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/sfml)
# # set_directory_properties(PROPERTIES
# #                          CMAKE_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
# #                          CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}
# #                         )
# # ExternalProject_Add(
# #    external_sfml
# #    URL file://${CMAKE_CURRENT_SOURCE_DIR}/sfml.zip
# #    CMAKE_ARGS
# #         -DCMAKE_BUILD_TYPE:STRING=Release
# #         -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
# #         -DBUILD_SHARED_LIBS:BOOL=OFF
# #    SOURCE_DIR ${SFML_SOURCE_DIR}
# #    BINARY_DIR ${SFML_BINARY_DIR}
# # )

# #FLTK
# file(DOWNLOAD http://fltk.org/pub/fltk/1.3.3/fltk-1.3.3-source.tar.gz
#               ${CMAKE_CURRENT_SOURCE_DIR}/fltk.tar.gz
#      EXPECTED_MD5 9ccdb0d19dc104b87179bd9fd10822e3
#     )
# set(FLTK_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/fltk)
# set(FLTK_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/fltk)
# set_directory_properties(PROPERTIES
#                          FLTK_SOURCE_DIR    ${FLTK_SOURCE_DIR}
#                          FLTK_BINARY_DIR    ${FLTK_BINARY_DIR}
#                          FLTK_INCLUDE_DIR   ${FLTK_SOURCE_DIR}
#                          FLTK_LIBRARY_DIR   ${FLTK_BINARY_DIR}/lib
#                          FLTK_LIBRARY_FILES "fltk;fltk_forms;fltk_gl"
#                         )
# ExternalProject_Add(
#    external_fltk
#    URL file://${CMAKE_CURRENT_SOURCE_DIR}/fltk.tar.gz
#    CMAKE_ARGS
#         -DCMAKE_BUILD_TYPE:STRING=Release
#         -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
#         -DBUILD_SHARED_LIBS:BOOL=OFF
#    SOURCE_DIR ${FLTK_SOURCE_DIR}
#    BINARY_DIR ${FLTK_BINARY_DIR}
# )
