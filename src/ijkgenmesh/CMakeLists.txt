PROJECT(IJKGENSCALAR)

#---------------------------------------------------------

SET(SRC_DIRNAME "src/ijkgenmesh")
SET(TARFILE ijkgenmesh.tar)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT DEFINED ${IJK_DIR})
  SET(IJK_DIR "../.." CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(NRRD_LIBDIR "${IJK_DIR}/lib")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")

ADD_EXECUTABLE(ijkgenmesh ijkgenmesh_main.cxx ijkgenmesh.cxx)

SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
INSTALL(TARGETS ijkgenmesh DESTINATION "bin/$ENV{OSTYPE}")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${SRC_DIRNAME}/${TARFILE}  ${SRC_DIRNAME}/*.cxx ${SRC_DIRNAME}/*.h ${SRC_DIRNAME}/*.txx ${SRC_DIRNAME}/CMakeLists.txt )

