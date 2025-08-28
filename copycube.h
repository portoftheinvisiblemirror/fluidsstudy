#pragma once
#ifndef FORCEANDTORQUE
#define FORCEANDTORQUE
#include <iostream>
#include <math.h>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include "utilities.hpp"
#include "mem.hpp"
#include "vector.hpp"
#include "tensor.hpp"
#include "vtk.hpp"
double*** compute_spheremesh(const double radius, const unsigned int latitudes, const unsigned int longitudes, double*** spheremesh);
double angleof(double x, double y);
double*** compute_spheremesh(const double radius, const unsigned int latitudes, const unsigned int longitudes);
static double areal(double radius, double xlow, double xup, double ylow, double yup, double zlow, double zup, unsigned int latitudes, unsigned int longitudes, double*** spheremesh, double sx, double sy, double sz, tensor*** st, const double dx, const double dy, const double dz, int x, int y, int z);
static vector forcel(double radius, double xlow, double xup, double ylow, double yup, double zlow, double zup, unsigned int latitudes, unsigned int longitudes, double*** spheremesh, double sx, double sy, double sz, tensor*** st, const double dx, const double dy, const double dz, int x, int y, int z);
static vector torquel(double radius, double xlow, double xup, double ylow, double yup, double zlow, double zup, unsigned int latitudes, unsigned int longitudes, double*** spheremesh, double sx, double sy, double sz, tensor*** st, const double dx, const double dy, const double dz, int x, int y, int z);
tensor gradv(double x, double y, double z, int nx, int ny, int nz, double dx, double dy, double dz, std::vector<std::vector<std::vector<std::vector<double>>>> a);
vector force(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const double nx, const int ny, const int nz, unsigned int longitudes, unsigned int latitudes, double radius, double*** spheremesh, std::vector<std::vector<int>> sphere, tensor*** st, const double R);
vector torque(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const double nx, const int ny, const int nz, unsigned int longitudes, unsigned int latitudes, double radius, double*** spheremesh, std::vector<std::vector<int>> sphere, tensor*** st, const double R);
double area(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const double nx, const int ny, const int nz, unsigned int longitudes, unsigned int latitudes, double radius, double*** spheremesh, std::vector<std::vector<std::vector<bool>>> sphere, tensor*** st, const double R);
vector* forceandtorque(std::vector<std::vector<int>> sphere, std::vector<std::vector<std::vector<double>>> u, std::vector<std::vector<std::vector<double>>> v, std::vector<std::vector<std::vector<double>>> w, std::vector<std::vector<std::vector<double>>> P, int nx, int ny, int nz, double dx, double dy, double dz, double x0, double y0, double z0, double R,double radius);
std::tuple<vector, vector> forceandtorquestag(std::vector<std::vector<int>> sphere, std::vector<std::vector<std::vector<double>>> u, std::vector<std::vector<std::vector<double>>> v, std::vector<std::vector<std::vector<double>>> w, std::vector<std::vector<std::vector<double>>> P, int nx, int ny, int nz, double dx, double dy, double dz, double x0, double y0, double z0, double R, double radius, double *** spheremesh, tensor *** stresses, int latitudes, int longitudes);
#endif






