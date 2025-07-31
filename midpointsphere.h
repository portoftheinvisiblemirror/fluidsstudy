#pragma once
#ifndef MIDPOINTSPHERE
#define MIDPOINTSPHERE
#include <vector>
#include <cmath>
#include <algorithm>
void fillin(std::vector<bool> row, int i1, int i2);
void filledmidpointcircle(const double y0, const double z0, const double dy, const double dz, const int ny, const int nz, const double radius, std::vector<std::vector<bool>> circle, const double R);
std::vector<std::vector<std::vector<bool>>> filledmidpointsphere(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const double nx, const int ny, const int nz, const double radius, const double R);
#endif