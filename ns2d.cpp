#include "ns2d.hpp"
#include <cmath>
#include <iostream>

void NavierStokesSolver2D::setBoundaryConditions() {
    // Set no-slip boundary conditions
    for (int i = 0; i < nx; i++) {
        u[i][0] = 0.0;  u[i][ny-1] = 0.0;
        v[i][0] = 0.0;  v[i][ny-1] = 0.0;
    }
    for (int j = 0; j < ny; j++) {
        u[0][j] = 0.0;  u[nx-1][j] = 0.0;
        v[0][j] = 0.0;  v[nx-1][j] = 0.0;
    }
}

void NavierStokesSolver2D::solvePoisson() {
    double beta = 1.0;
    int maxIter = 100;
    double tolerance = 1e-5;

    for (int iter = 0; iter < maxIter; iter++) {
        double maxDiff = 0.0;
        for (int i = 1; i < nx-1; i++) {
            for (int j = 1; j < ny-1; j++) {
                double b = rho * (
                    (u[i+1][j] - u[i-1][j])/(2*dx) +
                    (v[i][j+1] - v[i][j-1])/(2*dy)
                ) / dt;
                
                double pNew = (1.0-beta)*p[i][j] + 
                    beta/(2.0/(dx*dx) + 2.0/(dy*dy)) * (
                        (p[i+1][j] + p[i-1][j])/(dx*dx) +
                        (p[i][j+1] + p[i][j-1])/(dy*dy) - b
                    );
                
                maxDiff = std::max(maxDiff, std::abs(p[i][j] - pNew));
                p[i][j] = pNew;
            }
        }
        if (maxDiff < tolerance) break;
    }
}

NavierStokesSolver2D::NavierStokesSolver2D(int nx_, int ny_, double dx_, double dy_, 
                                          double dt_, double rho_, double mu_) 
    : nx(nx_), ny(ny_), dx(dx_), dy(dy_), dt(dt_), rho(rho_), mu(mu_) {
    
    u.resize(nx, std::vector<double>(ny, 0.0));
    v.resize(nx, std::vector<double>(ny, 0.0));
    p.resize(nx, std::vector<double>(ny, 0.0));
    tmp.resize(nx, std::vector<double>(ny, 0.0));
}

void NavierStokesSolver2D::step() {
    // Temporary arrays for velocity components
    std::vector<std::vector<double>> u_tmp = u;
    std::vector<std::vector<double>> v_tmp = v;

    // Solve momentum equations
    for (int i = 1; i < nx-1; i++) {
        for (int j = 1; j < ny-1; j++) {
            // u-momentum
            double du_dx = (u[i+1][j] - u[i-1][j])/(2*dx);
            double du_dy = (u[i][j+1] - u[i][j-1])/(2*dy);
            double d2u_dx2 = (u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dx*dx);
            double d2u_dy2 = (u[i][j+1] - 2*u[i][j] + u[i][j-1])/(dy*dy);
            
            u_tmp[i][j] = u[i][j] + dt*(
                -u[i][j]*du_dx - v[i][j]*du_dy -
                1/rho*(p[i+1][j] - p[i-1][j])/(2*dx) +
                mu/rho*(d2u_dx2 + d2u_dy2)
            );

            // v-momentum
            double dv_dx = (v[i+1][j] - v[i-1][j])/(2*dx);
            double dv_dy = (v[i][j+1] - v[i][j-1])/(2*dy);
            double d2v_dx2 = (v[i+1][j] - 2*v[i][j] + v[i-1][j])/(dx*dx);
            double d2v_dy2 = (v[i][j+1] - 2*v[i][j] + v[i][j-1])/(dy*dy);
            
            v_tmp[i][j] = v[i][j] + dt*(
                -u[i][j]*dv_dx - v[i][j]*dv_dy -
                1/rho*(p[i][j+1] - p[i][j-1])/(2*dy) +
                mu/rho*(d2v_dx2 + d2v_dy2)
            );
        }
    }

    u = u_tmp;
    v = v_tmp;

    // Solve pressure Poisson equation
    solvePoisson();

    // Apply boundary conditions
    setBoundaryConditions();
}

