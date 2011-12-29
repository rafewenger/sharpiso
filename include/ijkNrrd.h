// ijk functions for interface with Nrrd

#ifndef _IJKNRRD_
#define _IJKNRRD_

#include <iostream>

#include "NrrdIO.h"

namespace {

  template <typename T> void nrrd2scalar( Nrrd *nrrd, T * sdata)
  {
    std::cerr << "Programming error: Illegal instantiation of template nrrd2scalar."
	      << std::endl;
    exit(60);
  }

  template <typename T> 
  void nrrd_check_null(const int numv, T * sdata)
  {
    if (numv > 0 && sdata == NULL) {
      std::cerr << "Programming error detected in nrrd2scalar." << std::endl;
      std::cerr << "  Memory for scalar data not allocated." << std::endl;
      exit(70);
    }
  }

  template <> void nrrd2scalar<float>( Nrrd *nrrd, float * sdata )
  // convert nrrd to <float> data and store in sdata
  // nrrd = nrrd data structure
  // sdata = array of float
  //   Must be preallocated to size at least nrrdElementNumber(nrrd)
  {
    float (*lup) (const void *, size_t iv);

    void * nrrd_data = nrrd->data;
    int numv = nrrdElementNumber(nrrd);
    lup = nrrdFLookup[nrrd->type];

    nrrd_check_null(numv, sdata);

    for (int iv = 0; iv < numv; iv++)
      sdata[iv] = lup(nrrd_data, iv);
  }

  template <> void nrrd2scalar<double>( Nrrd *nrrd, double * sdata )
  // convert nrrd to <double> data and store in sdata
  // nrrd = nrrd data structure
  // sdata = array of float
  //   Must be preallocated to size at least nrrdElementNumber(nrrd)
  {
    double (*lup) (const void *, size_t iv);

    void * nrrd_data = nrrd->data;
    int numv = nrrdElementNumber(nrrd);
    lup = nrrdDLookup[nrrd->type];

    nrrd_check_null(numv, sdata);

    for (int iv = 0; iv < numv; iv++)
      sdata[iv] = lup(nrrd_data, iv);
  }

  template <> void nrrd2scalar<int>( Nrrd *nrrd, int * sdata )
  // convert nrrd to <int> data and store in sdata
  // nrrd = nrrd data structure
  // sdata = array of float
  //   Must be preallocated to size at least nrrdElementNumber(nrrd)
  {
    int (*lup) (const void *, size_t iv);

    void * nrrd_data = nrrd->data;
    int numv = nrrdElementNumber(nrrd);
    lup = nrrdILookup[nrrd->type];

    nrrd_check_null(numv, sdata);

    for (int iv = 0; iv < numv; iv++)
      sdata[iv] = lup(nrrd_data, iv);
  }

}

#endif
