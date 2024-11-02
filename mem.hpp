#ifndef MY_MEM_HPP
#define MY_MEM_HPP
// Header file for memory management related utilities definitions
#include <iostream>
#include <math.h>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <chrono>
#include "vector.hpp"
#include "tensor.hpp"
long double** allocate2d(const unsigned int x, const unsigned int y, long double** array2d);
long double*** allocate3d(const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea);
tensor*** allocate3dt(const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea);
long double**** allocate4d(const unsigned int nc, const unsigned int nx, const unsigned int ny, const unsigned int nz);
vector*** allocate3dv(const unsigned int nc, const unsigned int nx, const unsigned int ny);
#endif
