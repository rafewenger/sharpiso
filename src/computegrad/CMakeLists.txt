PROJECT(COMPUTEGRADIENT)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT DEFINED ${SHARP_DIR})
  GET_FILENAME_COMPONENT(SHARP_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(SHARP_DIR ${SHARP_ABSOLUTE_PATH} CACHE PATH "SHARP directory")
ENDIF (NOT DEFINED ${SHARP_DIR})

SET(CMAKE_INSTALL_PREFIX "${SHARP_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${SHARP_DIR}/lib CACHE PATH "Library directory")
SET(CGRADIENT_DIR "cgradient")
SET(NRRD_LIBDIR "${SHARP_DIR}/lib")
SET(SHARP_ISOTABLE_DIR "${SHARP_DIR}/isotable" CACHE PATH "Isotable directory")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

IF (CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "-std=c++0x")
ENDIF (CMAKE_COMPILER_IS_GNUCXX)

INCLUDE_DIRECTORIES("${SHARP_DIR}/include")
INCLUDE_DIRECTORIES("${SHARP_DIR}/src/sharpiso")
LINK_DIRECTORIES("${NRRD_LIBDIR}")

include(FindEXPAT)
find_package(EXPAT REQUIRED)
include_directories(${EXPAT_INCLUDE_DIRS})

find_library (ITKZLIB_LIBRARY ITKZLIB PATHS "${SHARPISO_DIR}/libs")
find_library (ZLIB_FOUND ZLIB PATHS "${SHARPISO_DIR}/libs")

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

LINK_LIBRARIES(NrrdIO ${LIB_ZLIB})
ADD_DEFINITIONS(-DSHARP_ISOTABLE_DIR=\"${SHARP_ISOTABLE_DIR}\")


ADD_EXECUTABLE(religrad religrad_main.cxx religrad_computations.cxx)
target_link_libraries(religrad ${EXPAT_LIBRARIES} NrrdIO ${LIB_ZLIB})

SET(CMAKE_INSTALL_PREFIX ${SHARP_DIR})
INSTALL(TARGETS religrad DESTINATION "/usr/local/bin")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${CGRADIENT_DIR}/isodual3D.tar ${CGRADIENT_DIR}/*.cxx ${CGRADIENT_DIR}/*.h ${CGRADIENT_DIR}/*.txx ${CGRADIENT_DIR}/CMakeLists.txt)

