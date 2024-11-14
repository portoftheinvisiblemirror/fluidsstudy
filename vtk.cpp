#include <iostream>
#include <fstream>
#include "vtk.hpp"

void writeVTKFile(int timestep, long double **** vp4d, int Nx, int Ny, int Nz, long double dx, long double dy, long double dz) {
    std::ofstream vtkFile;
    vtkFile.open("output_" + std::to_string(timestep) + ".vtk");
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
    std::cout << "VTK file output_" << timestep << ".vtk written.\n";
}
