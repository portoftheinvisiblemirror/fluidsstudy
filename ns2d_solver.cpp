#include "ns2d.hpp"
#include "vtk.hpp"
#include <cmath>
#include <iostream>

// Function to check if a point is inside the circular obstacle
bool isInsideObstacle(double x, double y, double centerX, double centerY, double radius) {
    double dx = x - centerX;
    double dy = y - centerY;
    return (dx*dx + dy*dy) <= radius*radius;
}

int main() {
    // Simulation parameters
    const int nx = 200;           // Grid points in x direction
    const int ny = 80;            // Grid points in y direction
    const double Lx = 10.0;       // Domain length
    const double Ly = 4.0;        // Domain height
    const double dx = Lx/nx;
    const double dy = Ly/ny;
    const double dt = 0.001;      // Time step
    const double rho = 1.0;       // Density
    const double mu = 0.01;       // Viscosity
    const int nSteps = 1000;      // Number of time steps
    const int saveInterval = 20;  // Save results every N steps

    // Obstacle parameters
    const double obstacleR = 0.2;                  // Radius of obstacle
    const double obstacleX = Lx/4;                 // X position of obstacle
    const double obstacleY = Ly/2;                 // Y position of obstacle

    // Initialize solver
    NavierStokesSolver2D solver(nx, ny, dx, dy, dt, rho, mu);

    // Set initial conditions - inlet velocity profile (parabolic)
    auto& u = const_cast<std::vector<std::vector<double>>&>(solver.getU());
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double y = j*dy;
            // Parabolic inlet profile
            if (i == 0) {
                u[i][j] = 4.0 * y/Ly * (1.0 - y/Ly);
            }
            
            // Set velocity to zero inside obstacle
            double x = i*dx;
            if (isInsideObstacle(x, y, obstacleX, obstacleY, obstacleR)) {
                u[i][j] = 0.0;
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
            const auto& p = solver.getP();
            
            std::string filename = "flow_" + std::to_string(step) + ".vtk";
            writeVTKFile(filename, u, v, p, dx, dy);
            
            std::cout << "Step " << step << " completed, saved to " << filename << "\n";
        }
    }

    return 0;
} 
