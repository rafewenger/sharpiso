PROJECT(MERGESHARP)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT DEFINED ${SHARPISO_DIR})
  GET_FILENAME_COMPONENT(SHARPISO_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(SHARPISO_DIR ${SHARPISO_ABSOLUTE_PATH} CACHE PATH "sharpiso directory")
ENDIF (NOT DEFINED ${SHARPISO_DIR})

SET(CMAKE_INSTALL_PREFIX "${SHARPISO_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${SHARPISO_DIR}/lib CACHE PATH "Library directory")
SET(EIGEN_DIR "${SHARPISO_DIR}/src/eigen" CACHE PATH "Eigen directory")
SET(NRRD_LIBDIR "${SHARPISO_DIR}/lib" CACHE PATH "Nrrd library directory")
SET(MERGESHARP_DIR "src/mergesharp")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${EIGEN_DIR}")
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

SET(CMAKE_INSTALL_PREFIX ${SHARPISO_DIR})
INSTALL(TARGETS mergesharp DESTINATION "bin/$ENV{OSTYPE}")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${MERGESHARP_DIR}/mergesharp.tar ${MERGESHARP_DIR}/*.cxx ${MERGESHARP_DIR}/*.h ${MERGESHARP_DIR}/*.txx ${MERGESHARP_DIR}/CMakeLists.txt ${MERGESHARP_DIR}/INSTALL ${MERGESHARP_DIR}/RELEASE_NOTES)

