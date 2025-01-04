#ifndef VTK_HPP
#define VTK_HPP

#include <vector>
#include <string>
#include <iostream>

// Original 3D version
void writeVTKFile(int timestep, double **** vp4d, int Nx, int Ny, int Nz, 
                 double dx, double dy, double dz);

// New 2D version
void writeVTKFile(const std::string& filename, const std::vector<std::vector<double>>& u,
                 const std::vector<std::vector<double>>& v,
                 const std::vector<std::vector<double>>& p,
                 double dx, double dy);

void writeVTKFile3D(const std::string& filename,
                    const std::vector<std::vector<std::vector<double>>>& u,
                    const std::vector<std::vector<std::vector<double>>>& v,
                    const std::vector<std::vector<std::vector<double>>>& w,
                    const std::vector<std::vector<std::vector<double>>>& p,
                    double dx, double dy, double dz);

#endif
