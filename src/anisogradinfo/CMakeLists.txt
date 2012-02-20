PROJECT(anisogradinfo)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT DEFINED ${SHARP_DIR})
  GET_FILENAME_COMPONENT(SHARP_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(SHARP_DIR ${SHARP_ABSOLUTE_PATH} CACHE PATH "SHARP directory")
ENDIF (NOT DEFINED ${SHARP_DIR})

SET(CMAKE_INSTALL_PREFIX "${SHARP_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${SHARP_DIR}/lib CACHE PATH "Library directory")
SET(anisograd_DIR "anisograd")
SET(NRRD_LIBDIR "${SHARP_DIR}/lib")
SET(SHARP_ISOTABLE_DIR "${SHARP_DIR}/isotable" CACHE PATH "Isotable directory")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${SHARP_DIR}/include")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(expat NrrdIO z)
ADD_DEFINITIONS(-DSHARP_ISOTABLE_DIR=\"${SHARP_ISOTABLE_DIR}\")

ADD_EXECUTABLE(anisogradinfo anisograd_main.cxx anisograd_operators.cxx  anisograd.cxx)

SET(CMAKE_INSTALL_PREFIX ${SHARP_DIR})
INSTALL(TARGETS anisogradinfo DESTINATION "bin/$ENV{OSTYPE}")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${anisograd_DIR}/isodual3D.tar ${anisograd_DIR}/*.cxx ${anisograd_DIR}/*.h ${anisograd_DIR}/*.txx ${anisograd_DIR}/CMakeLists.txt)

