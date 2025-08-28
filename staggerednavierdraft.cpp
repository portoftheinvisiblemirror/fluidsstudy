#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "midpointsphere.h"
#include "copycube.h"
#include <tuple>
void solve_momentum_equations();
void solve_pressure_correction();
void correct_pressure_velocity();
void check_convergence(double& max_Div, double& max_p_change);
void apply_boundary_conditions();
void output_results(size_t time_step, double time);
void output_final_results();


typedef double real;
//constants
size_t nx = 128, ny = 64, nz = 64;   // Smaller grid for testing
size_t max_iter = 50;             // Maximum time steps
size_t n_correctors = 4;            // Fewer correctors
real dt = 1e-3;          // Much smaller time step
real Re = 50.0;                  // Lower Reynolds number
real nu = 1.0 / Re;                 // Kinematic viscosity
real rho = 1.0;              // Density
real L = 5.0;                  // Shorter tube
real R = 1.0;                    // Tube radius
real dx = L / (nx - 1.);      // Grid spacing in x
real dy = 2.0 * R / (ny - 1.);   // Grid spacing in y
real dz = 2.0 * R / (nz - 1.); // Grid spacing in z
real tolerance = 1.0e-4;          // Relaxed tolerance
real alpha = 0.5;         // Under-relaxation for velocity
real Beta = 0.3;                // Under-relaxation for pressure
real cx = 2, cy = 0, cz = 0, r = 0.3; //center coordinates of a sphere, radius of a sphere
real mass = 1; //mass of the particle
// Arrays
std::vector<std::vector<std::vector<double>>> u, v, w;         // Velocity components
std::vector<std::vector<std::vector<double>>> p;               // Pressure
std::vector<std::vector<std::vector<double>>> u_old, v_old, w_old; // Old velocities
std::vector<std::vector<std::vector<double>>> p_old;           // Old pressure
std::vector<std::vector<std::vector<double>>> p_prime;         // Pressure correction
std::vector<std::vector<std::vector<double>>> Div;                           // Divergence
std::vector<double> x, y, z, xx, yy, zz;                                   // Grid coordinates
std::vector<std::vector<double>> radius, theta_coord;         // Cylindrical coordinates
std::vector<std::vector<std::vector<bool>>> sphere; //boolean ball
std::vector<std::vector<int>> sphere2; //indices of ball boundary
std::vector<std::vector<std::vector<bool>>> sphere3; //boolean ball boundary
void allocate_arrays()
{
    u.resize(nx + 1, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    v.resize(nx, std::vector<std::vector<double>>(ny + 1, std::vector<double>(nz, 0.0)));
    w.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz + 1, 0.0)));
    u_old = u;
    v_old = v;
    w_old = w;
    p.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    p_old = p;
    p_prime = p;
    Div = p;

    x.resize(nx, 0.0);
    y.resize(ny, 0.0);
    z.resize(nz, 0.0);
    xx.resize(nx + 1, 0.0);
    yy.resize(ny + 1, 0.0);
    zz.resize(nz + 1, 0.0);


    radius.resize(ny, std::vector<double>(nz, 0.0));
    theta_coord.resize(ny, std::vector<double>(nz, 0.0));
    sphere = filledmidpointsphere(cx, cy, cz, dx, dy, dz, nx, ny, nz, r, R);
    sphere2 = emptiedmidpointspherex(cx, cy, cz, dx, dy, dz, nx, ny, nz, r, R);
    sphere3 = emptiedmidpointsphere(cx, cy, cz, dx, dy, dz, nx, ny, nz, r, R);
}


void initialize_grid()
{
    // Initialize Cartesian grid coordinates at the center
    for (size_t i = 0; i < nx; ++i)
        x[i] = (i + 0.5) * dx;
    for (size_t i = 0; i < ny; ++i)
        y[i] = -R + (i + 0.5) * dy;
    for (size_t i = 0; i < nz; ++i)
        z[i] = -R + (i + 0.5) * dz;
    // Initialize Cartesian grid coordinates at each face
    for (size_t i = 0; i < nx + 1; ++i)
        xx[i] = (i)*dx;
    for (size_t i = 0; i < ny + 1; ++i)
        yy[i] = -R + (i)*dy;
    for (size_t i = 0; i < nz + 1; ++i)
        zz[i] = -R + (i)*dz;

    // Calculate cylindrical coordinates
    double x_center = 0.0;
    double y_center = 0.0;
    double z_center = 0.0;

    for (size_t j = 0; j < ny; ++j)
    {
        for (size_t k = 0; k < nz; ++k)
        {
            double r_dist = sqrt((y[j] - y_center) * (y[j] - y_center) + (z[k] - z_center) * (z[k] - z_center));
            radius[j][k] = r_dist;
            theta_coord[j][k] = atan2(z[k] - z_center, y[j] - y_center);
        }
    }
    std::cout << "Staggered Grid initialized : nx = " << nx << "ny = " << ny << "nz = " << nz;
    std::cout << "Grid spacing: dx=" << dx << " dy=" << dy << " dz =" << dz;
}

double findvmax(const std::vector< std::vector<std::vector<std::vector<double>>>>& v, const int N)
{
    double Q = 0;
    int x = 0;
    int y = 0;
    int z = 0;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int k = 0; k < N; ++k)
            {
                double mag = sqrt(v[0][i][j][k] * v[0][i][j][k] + v[1][i][j][k] * v[1][i][j][k] + v[2][i][j][k] * v[2][i][j][k]);
                if (mag > Q)
                {
                    x = i, y = j, z = k;
                    Q = mag;
                }
            }
        }
    }
    printf("maxspeed=%f at position (%d,%d,%d)\n", Q, x, y, z);
    return Q;
}
void initialize_flow()
{
    // Initialize velocity and pressure fields
    fill(u.begin(), u.end(), std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    fill(v.begin(), v.end(), std::vector<std::vector<double>>(ny + 1, std::vector<double>(nz, 0.0)));
    fill(w.begin(), w.end(), std::vector<std::vector<double>>(ny, std::vector<double>(nz + 1, 0.0)));
    fill(p.begin(), p.end(), std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));

    // Set inlet velocity profile (parabolic)
    double inlet_velocity = 0.5;  // Reduced inlet velocity
    for (size_t j = 0; j < ny; ++j)
        for (size_t k = 0; k < nz; ++k)
        {
            if (radius[j][k] <= R && !sphere[0][j][k])
            {
                u[0][j][k] = inlet_velocity * (1.0 - (radius[j][k] / R) * (radius[j][k] / R));
            }
            else
            {
                u[0][j][k] = 0;
            }
        }
    u_old = u;
    v_old = v;
    w_old = w;
    p_old = p;

    std::cout << "Flow field initialized with parabolic inlet profile\n";
    std::cout << "Inlet velocity: " << inlet_velocity << " Reynolds number : " << Re << "\n";
}

void solve_momentum_equations()
{
    // Local variables
    double conv_x, conv_y, conv_z;
    double diff_x, diff_y, diff_z;
    double dp_dx, dp_dy, dp_dz;

    // Solve u-momentum equation
    for (size_t i = 1; i < nx; ++i) {
        for (size_t j = 1; j < ny - 1; ++j) {
            for (size_t k = 1; k < nz - 1; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {  // Only solve inside tube
                    // Convection with an upwind
                    conv_x = u[i][j][k] * (u[i][j][k] - u[i - 1][j][k]) / dx;
                    conv_y = v[i][j][k] * (u[i][j + 1][k] - u[i][j - 1][k]) / (2.0 * dy);
                    conv_z = w[i][j][k] * (u[i][j][k + 1] - u[i][j][k - 1]) / (2.0 * dz);

                    // Diffusion terms
                    diff_x = nu * (u[i + 1][j][k] - 2.0 * u[i][j][k] + u[i - 1][j][k]) / (dx * dx);
                    diff_y = nu * (u[i][j + 1][k] - 2.0 * u[i][j][k] + u[i][j - 1][k]) / (dy * dy);
                    diff_z = nu * (u[i][j][k + 1] - 2.0 * u[i][j][k] + u[i][j][k - 1]) / (dz * dz);

                    // Pressure gradient
                    dp_dx = (p[i][j][k] - p[i - 1][j][k]) / (dx * rho);

                    // Update 
                    u[i][j][k] = u_old[i][j][k] + dt * (-conv_x - conv_y - conv_z + diff_x + diff_y + diff_z - dp_dx);
                }
            }
        }
    }

    // Solve v-momentum equation
    for (size_t i = 1; i < nx - 1; ++i) {
        for (size_t j = 1; j < ny; ++j) {
            for (size_t k = 1; k < nz - 1; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {
                    // Convection with an upwind
                    conv_x = u[i][j][k] * (v[i + 1][j][k] - v[i - 1][j][k]) / (2.0 * dx);
                    conv_y = v[i][j][k] * (v[i][j + 1][k] - v[i][j - 1][k]) / (2.0 * dy);
                    conv_z = w[i][j][k] * (v[i][j][k + 1] - v[i][j][k - 1]) / (2.0 * dz);

                    // Diffusion terms
                    diff_x = nu * (v[i + 1][j][k] - 2.0 * v[i][j][k] + v[i - 1][j][k]) / (dx * dx);
                    diff_y = nu * (v[i][j + 1][k] - 2.0 * v[i][j][k] + v[i][j - 1][k]) / (dy * dy);
                    diff_z = nu * (v[i][j][k + 1] - 2.0 * v[i][j][k] + v[i][j][k - 1]) / (dz * dz);

                    // Pressure gradient
                    dp_dy = (p[i][j][k] - p[i][j-1][k]) / (dy * rho); //upwinded?

                    // update 
                    v[i][j][k] = v_old[i][j][k] + dt * (-conv_x -conv_y - conv_z + diff_x + diff_y + diff_z - dp_dy);
                }
            }
        }
    }

    // Solve w-momentum equation
    for (size_t i = 1; i < nx - 1; ++i) {
        for (size_t j = 1; j < ny - 1; ++j) {
            for (size_t k = 1; k < nz; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {
                    // Convection with an upwind
                    conv_x = u[i][j][k] * (w[i + 1][j][k] - w[i - 1][j][k]) / (2.0 * dx);
                    conv_y = v[i][j][k] * (w[i][j + 1][k] - w[i][j - 1][k]) / (2.0 * dy);
                    conv_z = w[i][j][k] * (w[i][j][k + 1] - w[i][j][k - 1]) / (2.0 * dz);
                    // Diffusion terms
                    diff_x = nu * (w[i + 1][j][k] - 2.0 * w[i][j][k] + w[i - 1][j][k]) / (dx * dx);
                    diff_y = nu * (w[i][j + 1][k] - 2.0 * w[i][j][k] + w[i][j - 1][k]) / (dy * dy);
                    diff_z = nu * (w[i][j][k + 1] - 2.0 * w[i][j][k] + w[i][j][k - 1]) / (dz * dz);

                    // Pressure gradient
                    dp_dz = (p[i][j][k] - p[i][j][k-1]) / (dz * rho);

                    // update 
                    w[i][j][k] = w_old[i][j][k] + dt * (-conv_x- conv_y -conv_z + diff_x + diff_y + diff_z - dp_dz);
                }
            }
        }
    }
}
void solve_pressure_correction() {
    // Calculate Divergence
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {
                    Div[i][j][k] = (u[i + 1][j][k] - u[i][j][k]) / (dx) +
                        (v[i][j + 1][k] - v[i][j][k]) / (dy) +
                        (w[i][j][k + 1] - w[i][j][k]) / (dz);
                }
                else {
                    Div[i][j][k] = 0.0;
                }
            }
        }
    }
    double omega = 1.5;
    // Solve pressure correction equation using SOR
    for (size_t sor_iter = 1; sor_iter <= 2; ++sor_iter) {  // Fewer iterations
        for (size_t i = 1; i < nx - 1; ++i) {
            for (size_t j = 1; j < ny - 1; ++j) {
                for (size_t k = 1; k < nz - 1; ++k) {
                    if (radius[j][k] <= R && !sphere[i][j][k]) {
                        double c = 2 * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz *dz));
                        double p_new = (1.0 / c) * ((p_prime[i + 1][j][k] + p_prime[i - 1][j][k])/dx/dx +
                            (p_prime[i][j + 1][k] + p_prime[i][j - 1][k])/dy/dy +
                            (p_prime[i][j][k + 1] + p_prime[i][j][k - 1])/dz/dz -
                            Div[i][j][k]); 

                        p_prime[i][j][k] = p_prime[i][j][k] + omega * (p_new - p_prime[i][j][k]);
                    }
                }
            }
        }
    }
}

void correct_pressure_velocity() {
    // Correct pressure
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {
                    p[i][j][k] = p[i][j][k] + p_prime[i][j][k];
                }
            }
        }
    }

    // Correct velocities
    for (size_t i = 1; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {
                    u[i][j][k] = u[i][j][k] - dt * (p_prime[i][j][k] - p_prime[i - 1][j][k]) / (dx * rho);
                }
            }
        }
    }
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 1; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {
                    v[i][j][k] = v[i][j][k] - dt * (p_prime[i][j][k] - p_prime[i][j - 1][k]) / (dy * rho);
                }
            }
        }
    }
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 1; k < nz; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {
                    w[i][j][k] = w[i][j][k] - dt * (p_prime[i][j][k] - p_prime[i][j][k - 1]) / (dz * rho);
                }
            }
        }
    }
}
void check_convergence(double& max_Div, double& max_p_change) {
    max_Div = 0.0;
    max_p_change = 0.0;

    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {
                    max_Div = std::max(max_Div, std::abs(Div[i][j][k]));
                    max_p_change = std::max(max_p_change, std::abs(p_prime[i][j][k]));
                }
            }
        }
    }
}

void apply_boundary_conditions() {
    // Inlet boundary (parabolic velocity profile)
    double inlet_velocity = 0.5;  // Reduced inlet velocity
    for (size_t j = 0; j < ny; ++j) {
        for (size_t k = 0; k < nz; ++k) {
            if (radius[j][k] <= R && !sphere[0][j][k]) {
                u[0][j][k] = inlet_velocity * (1.0f - std::pow(radius[j][k] / R, 2));
                v[0][j][k] = 0.0f;
                w[0][j][k] = 0.0f;
            }
            else {
                u[0][j][k] = 0.0f;
                v[0][j][k] = 0.0f;
                w[0][j][k] = 0.0f;
            }
        }
    }

    // Outlet boundary (zero gradient)
    for (size_t j = 0; j < ny; ++j) {
        for (size_t k = 0; k < nz; ++k) {
            u[nx][j][k] = u[nx - 1][j][k];
            v[nx - 1][j+1][k] = v[nx - 1][j][k];
            w[nx - 1][j][k+1] = w[nx - 1][j][k];
            p[nx - 1][j][k] = 0.0;
        }
    }

    // Wall boundary (no-slip)
    for (size_t i = 0; i < nx+1; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                if (radius[j][k] > R  || (sphere[std::min(i, nx - 1)][j][k]&& !sphere3[std::min(i, nx - 1)][j][k])) {
                    u[i][j][k] = 0.0;
                }
            }
        }
    }
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny+1; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                if (radius[std::min(j, ny - 1)][k] > R || (sphere[i][std::min(j, ny - 1)][k] && !sphere3[i][std::min(j, ny - 1)][k])) {
                    v[i][j][k] = 0.0;
                }
            }
        }
    }
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz+1; ++k) {
                if (radius[j][std::min(k, nz - 1)] > R ||(sphere[i][j][std::min(k, nz - 1)] && !sphere3[i][j][std::min(k, nz - 1)])) {
                    w[i][j][k] = 0.0;
                }
            }
        }
    }
    // Pressure boundary conditions
    //inlet
    for (size_t j = 0; j < ny; ++j)
    {
        for (size_t k = 0; k < nz; ++k)
        {
            p[0][j][k] = 0.;
            //p[0][j][k] = p[1][j][k];
            //maybe make the inlet ghost cell pressure 0?
        }
    }
}
void output_results(size_t time_step, double time) {

    std::string filename = "stable_piso_velocity_" + std::to_string(time_step) + ".vtk";
    std::ofstream file(filename);

    file << "# vtk DataFile Version 3.0\n";
    file << "3D Navier-Stokes PISO Flow in Circular Tube (Stable)\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";

    file << "POINTS " << nx * ny * nz << " float\n";
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                file << std::fixed << std::setprecision(6) << x[i] << " "
                    << y[j] << " " << z[k] << "\n";
            }
        }
    }

    file << "POINT_DATA " << nx * ny * nz << "\n";
    file << "VECTORS velocity float\n";
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {

                file << std::fixed << std::setprecision(6) << (u[i][j][k] + u[i+1][j][k])/2 << " "
                    << (v[i][j][k] + v[i][j+1][k])/2 << " " << (w[i][j][k] + w[i][j][k+1])/2 << "\n";
            }
        }
    }

    file << "SCALARS pressure float\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                file << std::fixed << std::setprecision(6) << p[i][j][k] << "\n";
            }
        }
    }

    file.close();
}
void output_final_results() {
    // Local variables
    size_t count_points = 0;
    double max_velocity = 0.0, avg_velocity = 0.0, max_pressure = -1.0e10, min_pressure = 1.0e10, total_flow_rate = 0.0, vel_magnitude;

    // Output final results and statistics
    std::ofstream outfile("stable_piso_final_results.txt");

    outfile << "3D Navier-Stokes PISO Flow Solver Results (Stable Version)" << std::endl;
    outfile << "=========================================================" << std::endl;
    outfile << "Grid dimensions: " << nx << " x " << ny << " x " << nz << std::endl;
    outfile << "Reynolds number: " << Re << std::endl;
    outfile << "Time step: " << dt << std::endl;
    outfile << "Tube radius: " << R << std::endl;
    outfile << "Tube length: " << L << std::endl;
    outfile << "PISO correctors: " << n_correctors << std::endl;
    outfile << "Under-relaxation factors: alpha= " << alpha << " Beta= " << Beta << std::endl;

    // Calculate flow statistics
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                if (radius[j][k] <= R && !sphere[i][j][k]) {
                    double u_c = (u[i][j][k] + u[i + 1][j][k]) / 2;
                    double v_c = (v[i][j][k] + v[i][j + 1][k]) / 2;
                    double w_c = (w[i][j][k] + w[i][j][k + 1]) / 2;
                    vel_magnitude = std::sqrt(u_c*u_c +v_c*v_c+w_c*w_c);
                    max_velocity = std::max(max_velocity, vel_magnitude);
                    avg_velocity += vel_magnitude;
                    max_pressure = std::max(max_pressure, p[i][j][k]);
                    min_pressure = std::min(min_pressure, p[i][j][k]);
                    count_points++;

                    if (i == nx / 2) {  // Calculate flow rate at middle
                        total_flow_rate += u_c * dy * dz;
                    }
                }
            }
        }
    }

    if (count_points > 0) {
        avg_velocity /= static_cast<double>(count_points);
    }

    outfile << std::endl;
    outfile << "Flow Statistics:" << std::endl;
    outfile << "Maximum velocity: " << max_velocity << std::endl;
    outfile << "Average velocity: " << avg_velocity << std::endl;
    outfile << "Maximum pressure: " << max_pressure << std::endl;
    outfile << "Minimum pressure: " << min_pressure << std::endl;
    outfile << "Flow rate at middle: " << total_flow_rate << std::endl;

    outfile.close();

    std::cout << "Final results written to stable_piso_final_results.txt" << std::endl;
}
int main()
{
    // Variables
    size_t  i, j, k, iter, time_step, corrector;
    real  time, max_Div, max_p_change;
    real  x_center, y_center, r_dist, inlet_velocity;
    std::string filename;

    std::cout << "Initializing 3D Navier - Stokes PISO Solver(Stable Version)\n";
    std::cout << "Grid: " << nx << " x " << ny << " x " << nz << " dt: " << dt << " Re: " << Re << std::endl;

    // Allocate arrays
    allocate_arrays();

    // Initialize grid
    initialize_grid();

    // Initialize flow field
    initialize_flow();

    vector angv, rcm(cx,cy,cz), ucm;
    for (int time_step = 1; time_step <= max_iter; ++time_step) {
        double time = static_cast<double>(time_step) * dt;

        if (time_step % 10 == 0) {
            std::cout << "Time step: " << time_step << " Time: " << time << std::endl;
        }

        // Save old values
        u_old = u;
        v_old = v;
        w_old = w;
        p_old = p;

        // PISO algorithm
        for (int corrector = 1; corrector <= n_correctors; ++corrector) {
            // Step 1: Solve momentum equations for intermediate velocities
            solve_momentum_equations();

            // Step 2: Solve pressure correction equation
            solve_pressure_correction();
            
            // Step 3: Correct pressure and velocities
            correct_pressure_velocity();

            // Step 4: Check convergence for this corrector
            check_convergence(max_Div, max_p_change);

            if (time_step % 10 == 0) {
                std::cout << "  Corrector: " << corrector << " Max Div: " << max_Div << " Max P change: " << max_p_change << std::endl;
            }

            // If converged, exit corrector loop
            if (max_Div < tolerance && max_p_change < tolerance) {
                break;
            }
        }

        // Apply boundary conditions
        apply_boundary_conditions();
        for (auto& E : sphere2)
        {
            int i = E[0], j = E[1], k = E[2];
            vector position = { i * dx,-R + (j)*dy, -R + (k)*dz };
            vector surfacevelocity = ucm + angv % position;
            u[i][j][k] = surfacevelocity.X(), v[i][j][k] = surfacevelocity.Y(), w[i][j][k] = surfacevelocity.Z();
        }
        // Output results periodically
        if (time_step % 50 == 0) {
            output_results(time_step, time);
        }

        // Check for instability
        if (max_Div > 20.0 || max_p_change > 10.0) {
            std::cout << "WARNING: Simulation becoming unstable" << std::endl;
            std::cout << "Max Div: " << max_Div << " Max P change: " << max_p_change << std::endl;
            std::cout << "Stopping at time step: " << time_step << std::endl;
            output_results(time_step, time);
            break;
        }
        //insert sphere code here
        auto start = std::chrono::high_resolution_clock::now();
        
        vector F, T;
        std::tie(F,T)= forceandtorquestag(sphere2, u, v, w, p, nx, ny, nz, dx, dy, dz, cx, cy, cz, R, r);
        angv = angv + T * dt / (2 * mass * r * r / 5);
        ucm = ucm + F * dt / mass;
        rcm = rcm + ucm * dt + F * dt * dt / (2 * mass);
        cx = rcm.X(), cy = rcm.Y(), cz = rcm.Z();
        sphere = filledmidpointsphere(cx, cy, cz, dx, dy, dz, nx, ny, nz, r, R);
        sphere2 = emptiedmidpointspherex(cx, cy, cz, dx, dy, dz, nx, ny, nz, r, R);
        sphere3 = emptiedmidpointsphere(cx, cy, cz, dx, dy, dz, nx, ny, nz, r, R);
        std::cout << cx << " " << cy << " " << cz << "\n";
        for (auto &E : sphere2)
        {
            int i = E[0], j = E[1], k = E[2];
            vector position = { i * dx,-R + (j)*dy, -R + (k)*dz };
            vector surfacevelocity = ucm + angv % position;
            u[i][j][k] = surfacevelocity.X(), v[i][j][k] = surfacevelocity.Y(), w[i][j][k] = surfacevelocity.Z();
        }

        // Stop the timer
        auto end = std::chrono::high_resolution_clock::now();

        // Calculate the duration
        auto duration = end - start;

        // Convert and display the duration in microseconds
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
        std::cout << "Execution time: " << microseconds << " microseconds" << std::endl;
               
    }

    // Final output
    output_final_results();

    std::cout << "Simulation completed successfully" << std::endl;

    return 0;
}


/* divergence()
* pressure correct
* velocity correct
* converged?
* apply boundary conditions
* output results 
* c
*/
