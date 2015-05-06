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
SET(SHARPISO_SRC_DIR "${SHARPISO_DIR}/src/sharpiso/" CACHE PATH "sharpiso src directory")
SET(NRRD_LIBDIR "${SHARPISO_DIR}/lib" CACHE PATH "Nrrd library directory")
SET(MERGESHARP_DIR "src/mergesharp")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

#Find Expat
include(FindEXPAT)
find_package(EXPAT REQUIRED)
include_directories(${EXPAT_INCLUDE_DIRS})

find_library (ITKZLIB_LIBRARY ITKZLIB PATHS "${SHARPISO_DIR}/lib")
find_library (ZLIB_FOUND ZLIB PATHS "${SHARPISO_DIR}/lib")

#Find itkzlib or z library
if ("${ITKZLIB_LIBRARY}" STREQUAL "ITKZLIB_LIBRARY-NOTFOUND")
  if ("${ZLIB_LIBRARY}" STREQUAL "ZLIB_LIBRARY-NOTFOUND")
    message ("ITKZLIB and ZLIB not found.  Exiting.")
    return(10)
  else()
    set (LIB_ZLIB "z")
  endif()
else()
  message ("ITKZLIB FOUND")
  set (LIB_ZLIB "ITKZLIB")
endif()


INCLUDE_DIRECTORIES("${EIGEN_DIR}")
INCLUDE_DIRECTORIES("${SHARPISO_SRC_DIR}")
INCLUDE_DIRECTORIES("${SHARPISO_DIR}/include")
LINK_DIRECTORIES("${NRRD_LIBDIR}")

IF(MSVC AND MSVC_VERSION GREATER 1600)
# Nothing special needed for MSVC
# Just enable all warnings for debug and set NDEBUG for release
set(CMAKE_CXX_FLAGS_DEBUG "/W4")
set(CMAKE_CXX_FLAGS_RELEASE "/D NDEBUG")
  message(
  "Your compiler probably does support C++11. This project requires C++11")
ELSEIF(CMAKE_COMPILER_IS_GNUCXX)
SET(CMAKE_CXX_FLAGS "-std=c++0x")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
SET(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")
else()
  message(FATAL_ERROR
  "Your compiler probably does not support C++11. This project requires C++11")
ENDIF()

SET(MERGESHARP_SUB_LIST mergesharpIO.cxx mergesharp.cxx 
                        mergesharp_datastruct.cxx 
                        mergesharp_isovert.cxx mergesharp_select.cxx
                        mergesharp_extract.cxx mergesharp_position.cxx 
                        mergesharp_merge.cxx mergesharp_check_map.cxx
                        ijkdualtable.cxx ijkdualtable_ambig.cxx 
                        ijktable_poly.cxx
                        ijktable_ambig.cxx mergesharp_ambig.cxx
                        ${SHARPISO_SRC_DIR}/sharpiso_feature.cxx
                        ${SHARPISO_SRC_DIR}/sharpiso_get_gradients.cxx
                        ${SHARPISO_SRC_DIR}/sharpiso_intersect.cxx
                        ${SHARPISO_SRC_DIR}/sharpiso_svd.cxx
                        ${SHARPISO_SRC_DIR}/sharpiso_closest.cxx)

ADD_EXECUTABLE(mergesharp mergesharp_main.cxx  ${MERGESHARP_SUB_LIST} )
target_link_libraries(mergesharp ${EXPAT_LIBRARIES} NrrdIO ${LIB_ZLIB})

SET(CMAKE_INSTALL_PREFIX ${SHARPISO_DIR})
INSTALL(TARGETS mergesharp DESTINATION "bin/$ENV{OSTYPE}")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${MERGESHARP_DIR}/mergesharp.tar ${MERGESHARP_DIR}/*.cxx ${MERGESHARP_DIR}/*.h ${MERGESHARP_DIR}/CMakeLists.txt ${MERGESHARP_DIR}/INSTALL ${MERGESHARP_DIR}/RELEASE_NOTES)

