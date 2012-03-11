PROJECT(ISODUAL3D)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(ISODUAL3D_DIR "src/isodual3D")
SET(NRRD_LIBDIR "${IJK_DIR}/lib")
SET(IJK_ISOTABLE_DIR "${IJK_DIR}/isotable" CACHE PATH "Isotable directory")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")
INCLUDE_DIRECTORIES("${IJK_DIR}/src/sharpiso")
INCLUDE_DIRECTORIES("../eigen")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(expat NrrdIO z)
ADD_DEFINITIONS(-DIJK_ISOTABLE_DIR=\"${IJK_ISOTABLE_DIR}\")

ADD_EXECUTABLE(isodual3D isodual3D_main.cxx isodual3DIO.cxx isodual3D.cxx isodual3D_datastruct.cxx
isodual3D_extract.cxx isodual3D_position.cxx sharpiso_feature.cxx sharpiso_feature.cxx sharpiso_svd.cxx
sh_point_find.cxx sh_point_datastruct.cxx sh_point_svdcal.cxx
sharpiso_findIntersect.cxx)
SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
INSTALL(TARGETS isodual3D DESTINATION "bin/$ENV{OSTYPE}")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${ISODUAL3D_DIR}/isodual3D.tar ${ISODUAL3D_DIR}/*.cxx ${ISODUAL3D_DIR}/*.h ${ISODUAL3D_DIR}/*.txx ${ISODUAL3D_DIR}/CMakeLists.txt)
