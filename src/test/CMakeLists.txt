PROJECT(TEST_SHARPINFO)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT DEFINED ${SHARP_DIR})
  GET_FILENAME_COMPONENT(SHARP_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(SHARP_DIR ${SHARP_ABSOLUTE_PATH} CACHE PATH "SHARP directory")
ENDIF (NOT DEFINED ${SHARP_DIR})

SET(CMAKE_INSTALL_PREFIX "${SHARP_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${SHARP_DIR}/lib CACHE PATH "Library directory")
SET(SHARPINFO_DIR "src/sharpinfo")
SET(NRRD_LIBDIR "${SHARP_DIR}/lib")
SET(SHARP_ISOTABLE_DIR "${SHARP_DIR}/isotable" CACHE PATH "Isotable directory")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${SHARP_DIR}/include")
INCLUDE_DIRECTORIES("../eigen")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(ITKNrrdIO z)

ADD_EXECUTABLE(testscalar testscalar.cxx)

ADD_EXECUTABLE(testfeature testfeature.cxx sharpiso_feature.cxx sharpiso_intersect.cxx sharpiso_get_gradients.cxx sharpiso_svd.cxx sharpiso_findIntersect.cxx sh_point_find.cxx sh_point_datastruct.cxx sh_point_svdcal.cxx )

ADD_EXECUTABLE(test_scalar_func test_scalar_func.cxx)

ADD_EXECUTABLE(test_contains test_contains.cxx)

ADD_EXECUTABLE(testdist testdist.cxx)

ADD_EXECUTABLE(testintersect testintersect.cxx sharpiso_findIntersect.cxx)

ADD_EXECUTABLE(testambig testambig.cxx ijktable_poly.cxx ijktable_ambig.cxx
                         isodual3D_ambig.cxx 
                         sharpiso_feature.cxx sharpiso_get_gradients.cxx
                         sharpiso_svd.cxx sharpiso_intersect.cxx
                         sharpiso_findIntersect.cxx sh_point_find.cxx
                         sh_point_datastruct.cxx sh_point_svdcal.cxx)

ADD_EXECUTABLE(testbin testbin.cxx)
