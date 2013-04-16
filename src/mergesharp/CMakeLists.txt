PROJECT(MERGESHARP)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(MERGESHARP_DIR "src/mergesharp")
SET(NRRD_LIBDIR "${IJK_DIR}/lib")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("../eigen")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(expat NrrdIO z)

SET(CMAKE_CXX_FLAGS "-std=c++0x")

SET(MERGESHARP_SUB_LIST mergesharpIO.cxx mergesharp.cxx 
                        mergesharp_datastruct.cxx mergesharp_isovert.cxx
                        mergesharp_extract.cxx mergesharp_position.cxx 
                        mergesharp_merge.cxx
                        ijkdualtable.cxx ijkdualtable_ambig.cxx 
                        ijktable_poly.cxx
                        ijktable_ambig.cxx mergesharp_ambig.cxx
                        sharpiso_feature.cxx sharpiso_get_gradients.cxx
                        sharpiso_intersect.cxx
                        sharpiso_svd.cxx
                        sharpiso_closest.cxx)

ADD_EXECUTABLE(mergesharp mergesharp_main.cxx ${MERGESHARP_SUB_LIST} )

SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
INSTALL(TARGETS mergesharp DESTINATION "bin/$ENV{OSTYPE}")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${MERGESHARP_DIR}/isodual3D.tar ${MERGESHARP_DIR}/*.cxx ${MERGESHARP_DIR}/*.h ${MERGESHARP_DIR}/*.txx ${MERGESHARP_DIR}/CMakeLists.txt)
