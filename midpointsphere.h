#pragma once
#ifndef MIDPOINTSPHERE
#define MIDPOINTSPHERE
#include <vector>
#include <cmath>
#include <algorithm>
void fillin(std::vector<std::vector<int>>& row, int i1, int i2, int X, int Y);
void filledmidpointcircle(const double y0, const double z0, const double dy, const double dz, const int ny, const int nz, const double radius, std::vector<std::vector<int>>& circle, const double R, int X);
std::vector<std::vector<int>> filledmidpointspherex(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz, const double radius, const double R);
void emptiedinside(std::vector<std::vector<int>>& row, int i1, int i2, int ii1, int ii2, int X, int Y);
void emptiedmidpointcircle(const double y0, const double z0, const double dy, const double dz, const int ny, const int nz, const double radius, const double radlow, std::vector<std::vector<int>>& circle, const double R, int X);
std::vector<std::vector<int>> emptiedmidpointspherex(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz, const double radius, const double R);
void fillin(std::vector<bool>& row, int i1, int i2);
void filledmidpointcircle(const double y0, const double z0, const double dy, const double dz, const int ny, const int nz, const double radius, std::vector<std::vector<bool>>& circle, const double R);
std::vector<std::vector<std::vector<bool>>> filledmidpointsphere(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz, const double radius, const double R);
void emptiedinside(std::vector<bool>& row, int i1, int i2, int ii1, int ii2);
void emptiedmidpointcircle(const double y0, const double z0, const double dy, const double dz, const int ny, const int nz, const double radius, const double radlow, std::vector<std::vector<bool>>& circle, const double R);
std::vector<std::vector<std::vector<bool>>> emptiedmidpointsphere(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz, const double radius, const double R);
#endif