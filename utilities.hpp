#ifndef MY_UTILITIES_HPP
#define MY_UTILITIES_HPP
#include "tensor.hpp"
#include "vector.hpp"
tensor shear(double x, double y, double z);
tensor stress(double p, double x, double y, double z);
vector diff(double x, double y, double z, int index);
vector velocity(double x, double y, double z);
vector divvadvecv(double x, double y, double z);
#endif
