#include <iostream>
#include <fstream>
#include "vtk.hpp"

void writeVTKFile(int timestep, long double **** vp4d, int Nx, int Ny, int Nz, long double dx, long double dy, long double dz) {
    std::ofstream vtkFile;
    std::string filename = "output_" + std::to_string(timestep) + ".vtk";
    vtkFile.open(filename);
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "SIMPLE Solver Output\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";
    vtkFile << "DIMENSIONS " << Nx << " " << Ny <<  " " << Nz << "\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "SPACING " << dx << " " << dy << " " << dz << "\n";
    
    // Writing Pressure field
    vtkFile << "POINT_DATA " << Nx * Ny * Nz << "\n";
    vtkFile << "SCALARS pressure double\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; ++k) {
			for (int j = 0; j < Ny; ++j) {
					for (int i = 0; i < Nx; ++i) {
							vtkFile << vp4d[3][i][j][k] << " ";
					}
			}
			vtkFile << "\n";
		}

    // Writing Velocity field
    vtkFile << "VECTORS velocity double\n";
    for (int k = 0; k < Nz; ++k) {
			for (int j = 0; j < Ny; ++j) {
					for (int i = 0; i < Nx; ++i) {
							vtkFile << vp4d[0][i][j][k] << " " << vp4d[1][i][j][k] << " " << vp4d[2][i][j][k] << "\n";
					}
			}
		}

    vtkFile.close();
    std::cout << "VTK file " << filename << " written.\n";
}

void writeVTKFile(const std::string& filename, const std::vector<std::vector<double>>& u,
                 const std::vector<std::vector<double>>& v,
                 const std::vector<std::vector<double>>& p,
                 double dx, double dy) {
    std::ofstream vtkFile;
    vtkFile.open(filename);
    
    int Nx = u.size();
    int Ny = u[0].size();
    
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "SIMPLE Solver Output\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";
    vtkFile << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "SPACING " << dx << " " << dy << " 1\n";
    
    // Writing Pressure field
    vtkFile << "POINT_DATA " << Nx * Ny << "\n";
    vtkFile << "SCALARS pressure double\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtkFile << p[i][j] << " ";
        }
        vtkFile << "\n";
    }

    // Writing Velocity field
    vtkFile << "VECTORS velocity double\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtkFile << u[i][j] << " " << v[i][j] << " 0\n";
        }
    }

    vtkFile.close();
    std::cout << "VTK file " << filename << " written.\n";
}

void writeVTKFile3D(const std::string& filename,
                    const std::vector<std::vector<std::vector<double>>>& u,
                    const std::vector<std::vector<std::vector<double>>>& v,
                    const std::vector<std::vector<std::vector<double>>>& w,
                    const std::vector<std::vector<std::vector<double>>>& p,
                    double dx, double dy, double dz) {
    std::ofstream file(filename);
    
    // VTK header
    file << "# vtk DataFile Version 3.0\n";
    file << "Flow field\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    
    int nx = u.size();
    int ny = u[0].size();
    int nz = u[0][0].size();
    
    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    file << "POINTS " << nx*ny*nz << " float\n";
    
    // Write grid points
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                file << i*dx << " " << j*dy << " " << k*dz << "\n";
            }
        }
    }
    
    // Write velocity and pressure data
    file << "POINT_DATA " << nx*ny*nz << "\n";
    
    // Velocity vector
    file << "VECTORS velocity float\n";
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                file << u[i][j][k] << " " << v[i][j][k] << " " << w[i][j][k] << "\n";
            }
        }
    }
    
    // Pressure scalar
    file << "SCALARS pressure float\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                file << p[i][j][k] << "\n";
            }
        }
    }
}
