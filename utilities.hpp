#ifndef MY_UTILITIES_HPP
#define MY_UTILITIES_HPP
#include "tensor.hpp"
#include "vector.hpp"
tensor shear(double x, double y, double z);
tensor stress(double p, double x, double y, double z);
vector diff(double x, double y, double z, int index);
vector velocity(double x, double y, double z);
double divvadvecv(double **** a, double x, double y, double z,int n, double h);
vector*** divtenall(tensor*** a, const int n, double h);
double partdif(double**** a, int x, int y, int z, int n, double h, int d, int e);
#endif
