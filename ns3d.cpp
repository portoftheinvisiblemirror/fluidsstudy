#include "ns3d.hpp"
#include <cmath>

NavierStokesSolver3D::NavierStokesSolver3D(int nx_, int ny_, int nz_, double dx_, double dy_, double dz_, 
                                          double dt_, double rho_, double mu_)
    : nx(nx_), ny(ny_), nz(nz_), dx(dx_), dy(dy_), dz(dz_), dt(dt_), rho(rho_), mu(mu_) {
    
    // Initialize velocity and pressure fields
    u = std::vector<std::vector<std::vector<double>>>(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    v = std::vector<std::vector<std::vector<double>>>(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    w = std::vector<std::vector<std::vector<double>>>(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    p = std::vector<std::vector<std::vector<double>>>(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    
    // Initialize temporary arrays
    u_temp = u;
    v_temp = v;
    w_temp = w;
}

void NavierStokesSolver3D::solvePoisson() {
    const int maxIter = 50;
    const double tolerance = 1e-4;
    
    for (int iter = 0; iter < maxIter; iter++) {
        double maxDiff = 0.0;
        
        for (int i = 1; i < nx-1; i++) {
            for (int j = 1; j < ny-1; j++) {
                for (int k = 1; k < nz-1; k++) {
                    double x = i*dx;
                    double y = j*dy;
                    double z = k*dz;
                    
                    // Skip pressure calculation inside obstacle
                    if (isInsideObstacle(x, y, z)) {
                        continue;
                    }

                    double div = (u[i+1][j][k] - u[i-1][j][k])/(2*dx) +
                                (v[i][j+1][k] - v[i][j-1][k])/(2*dy) +
                                (w[i][j][k+1] - w[i][j][k-1])/(2*dz);
                    
                    double p_new = ((p[i+1][j][k] + p[i-1][j][k])/(dx*dx) +
                                  (p[i][j+1][k] + p[i][j-1][k])/(dy*dy) +
                                  (p[i][j][k+1] + p[i][j][k-1])/(dz*dz) -
                                  rho*div/dt) / (2/(dx*dx) + 2/(dy*dy) + 2/(dz*dz));
                    
                    maxDiff = std::max(maxDiff, std::abs(p_new - p[i][j][k]));
                    p[i][j][k] = p_new;
                }
            }
        }
        
        if (maxDiff < tolerance) break;
    }
}

void NavierStokesSolver3D::applyBoundaryConditions() {
    // No-slip boundary conditions for velocities
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            // Inlet (prescribed velocity)
            u[0][j][k] = u[1][j][k];
            v[0][j][k] = 0.0;
            w[0][j][k] = 0.0;
            
            // Outlet (zero gradient)
            u[nx-1][j][k] = u[nx-2][j][k];
            v[nx-1][j][k] = v[nx-2][j][k];
            w[nx-1][j][k] = w[nx-2][j][k];
        }
    }
    
    // Top and bottom walls
    for (int i = 0; i < nx; i++) {
        for (int k = 0; k < nz; k++) {
            u[i][0][k] = 0.0;
            v[i][0][k] = 0.0;
            w[i][0][k] = 0.0;
            
            u[i][ny-1][k] = 0.0;
            v[i][ny-1][k] = 0.0;
            w[i][ny-1][k] = 0.0;
        }
    }
    
    // Front and back walls
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            u[i][j][0] = 0.0;
            v[i][j][0] = 0.0;
            w[i][j][0] = 0.0;
            
            u[i][j][nz-1] = 0.0;
            v[i][j][nz-1] = 0.0;
            w[i][j][nz-1] = 0.0;
        }
    }
}

void NavierStokesSolver3D::step() {
    // Store current velocities
    u_temp = u;
    v_temp = v;
    w_temp = w;
    
    // Solve momentum equations
    for (int i = 1; i < nx-1; i++) {
        for (int j = 1; j < ny-1; j++) {
            for (int k = 1; k < nz-1; k++) {
                double x = i*dx;
                double y = j*dy;
                double z = k*dz;
                
                // Skip calculations inside obstacle
                if (isInsideObstacle(x, y, z)) {
                    u[i][j][k] = 0.0;
                    v[i][j][k] = 0.0;
                    w[i][j][k] = 0.0;
                    continue;
                }

                // Convective terms
                double u_dx = (u_temp[i+1][j][k] - u_temp[i-1][j][k])/(2*dx);
                double u_dy = (u_temp[i][j+1][k] - u_temp[i][j-1][k])/(2*dy);
                double u_dz = (u_temp[i][j][k+1] - u_temp[i][j][k-1])/(2*dz);
                
                double v_dx = (v_temp[i+1][j][k] - v_temp[i-1][j][k])/(2*dx);
                double v_dy = (v_temp[i][j+1][k] - v_temp[i][j-1][k])/(2*dy);
                double v_dz = (v_temp[i][j][k+1] - v_temp[i][j][k-1])/(2*dz);
                
                double w_dx = (w_temp[i+1][j][k] - w_temp[i-1][j][k])/(2*dx);
                double w_dy = (w_temp[i][j+1][k] - w_temp[i][j-1][k])/(2*dy);
                double w_dz = (w_temp[i][j][k+1] - w_temp[i][j][k-1])/(2*dz);
                
                // Viscous terms
                double u_lap = (u_temp[i+1][j][k] - 2*u_temp[i][j][k] + u_temp[i-1][j][k])/(dx*dx) +
                             (u_temp[i][j+1][k] - 2*u_temp[i][j][k] + u_temp[i][j-1][k])/(dy*dy) +
                             (u_temp[i][j][k+1] - 2*u_temp[i][j][k] + u_temp[i][j][k-1])/(dz*dz);
                
                double v_lap = (v_temp[i+1][j][k] - 2*v_temp[i][j][k] + v_temp[i-1][j][k])/(dx*dx) +
                             (v_temp[i][j+1][k] - 2*v_temp[i][j][k] + v_temp[i][j-1][k])/(dy*dy) +
                             (v_temp[i][j][k+1] - 2*v_temp[i][j][k] + v_temp[i][j][k-1])/(dz*dz);
                
                double w_lap = (w_temp[i+1][j][k] - 2*w_temp[i][j][k] + w_temp[i-1][j][k])/(dx*dx) +
                             (w_temp[i][j+1][k] - 2*w_temp[i][j][k] + w_temp[i][j-1][k])/(dy*dy) +
                             (w_temp[i][j][k+1] - 2*w_temp[i][j][k] + w_temp[i][j][k-1])/(dz*dz);
                
                // Update velocities
                u[i][j][k] = u_temp[i][j][k] + dt*(
                    -u_temp[i][j][k]*u_dx - v_temp[i][j][k]*u_dy - w_temp[i][j][k]*u_dz -
                    (1/rho)*(p[i+1][j][k] - p[i-1][j][k])/(2*dx) + 
                    mu*(u_lap)
                );
                
                v[i][j][k] = v_temp[i][j][k] + dt*(
                    -u_temp[i][j][k]*v_dx - v_temp[i][j][k]*v_dy - w_temp[i][j][k]*v_dz -
                    (1/rho)*(p[i][j+1][k] - p[i][j-1][k])/(2*dy) + 
                    mu*(v_lap)
                );
                
                w[i][j][k] = w_temp[i][j][k] + dt*(
                    -u_temp[i][j][k]*w_dx - v_temp[i][j][k]*w_dy - w_temp[i][j][k]*w_dz -
                    (1/rho)*(p[i][j][k+1] - p[i][j][k-1])/(2*dz) + 
                    mu*(w_lap)
                );
            }
        }
    }
    
    // Solve pressure Poisson equation
    solvePoisson();
    
    // Apply boundary conditions
    applyBoundaryConditions();
}

void NavierStokesSolver3D::setObstacle(double x, double y, double z, double r) {
    obstacleX = x;
    obstacleY = y;
    obstacleZ = z;
    obstacleR = r;
    hasObstacle = true;
}

bool NavierStokesSolver3D::isInsideObstacle(double x, double y, double z) const {
    if (!hasObstacle) return false;
    double dx = x - obstacleX;
    double dy = y - obstacleY;
    double dz = z - obstacleZ;
    return (dx*dx + dy*dy + dz*dz) <= obstacleR*obstacleR;
} 