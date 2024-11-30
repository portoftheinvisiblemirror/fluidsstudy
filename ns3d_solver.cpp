#include "ns3d.hpp"
#include "vtk.hpp"
#include <cmath>
#include <iostream>

// Function to check if a point is inside the circular obstacle
bool isInsideObstacle(double x, double y, double z, double centerX, double centerY, double centerZ, double radius) {
    double dx = x - centerX;
    double dy = y - centerY;
    double dz = z - centerZ;
    return (dx*dx + dy*dy + dz*dz) <= radius*radius;
}

int main() {
    // Simulation parameters
    const int nx = 200;           // Grid points in x direction
    const int ny = 80;            // Grid points in y direction
    const int nz = 80;            // Grid points in z direction
    const double Lx = 10.0;       // Domain length
    const double Ly = 4.0;        // Domain height
    const double Lz = 4.0;        // Domain depth
    const double dx = Lx/nx;
    const double dy = Ly/ny;
    const double dz = Lz/nz;      // Grid spacing in z
    const double dt = 0.001;      // Time step
    const double rho = 1.0;       // Density
    const double mu = 0.01;       // Viscosity
    const int nSteps = 1000;      // Number of time steps
    const int saveInterval = 20;  // Save results every N steps

    // Obstacle parameters
    const double obstacleR = 0.2;                  // Radius of obstacle
    const double obstacleX = Lx/4;                 // X position of obstacle
    const double obstacleY = Ly/2;                 // Y position of obstacle
    const double obstacleZ = Lz/2;                // Z position of obstacle

    // Initialize solver
    NavierStokesSolver3D solver(nx, ny, nz, dx, dy, dz, dt, rho, mu);

    // Set initial conditions - inlet velocity profile (parabolic)
    auto& u = const_cast<std::vector<std::vector<std::vector<double>>>&>(solver.getU());
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double y = j*dy;
                double z = k*dz;
                // Parabolic inlet profile
                if (i == 0) {
                    double r2 = pow(y - Ly/2, 2) + pow(z - Lz/2, 2);
                    u[i][j][k] = 1.0 * (1.0 - 4*r2/(Ly*Ly));  // Circular inlet profile
                }
                
                // Set velocity to zero inside obstacle
                double x = i*dx;
                if (isInsideObstacle(x, y, z, obstacleX, obstacleY, obstacleZ, obstacleR)) {
                    u[i][j][k] = 0.0;
                }
            }
        }
    }

    // Time stepping
    for (int step = 0; step < nSteps; step++) {
        // Solve one time step
        solver.step();

        // Save results periodically
        if (step % saveInterval == 0) {
            const auto& u = solver.getU();
            const auto& v = solver.getV();
            const auto& w = solver.getW();    // W velocity component
            const auto& p = solver.getP();
            
            std::string filename = "flow_" + std::to_string(step) + ".vtk";
            writeVTKFile3D(filename, u, v, w, p, dx, dy, dz);
            
            std::cout << "Step " << step << " completed, saved to " << filename << "\n";
        }
    }

    return 0;
} 
