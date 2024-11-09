#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

const int Nx = 50;  // Number of cells in the x-direction
const int Ny = 50;  // Number of cells in the y-direction
const double Lx = 1.0;  // Domain length in x-direction
const double Ly = 1.0;  // Domain length in y-direction
const double rho = 1.0;  // Density
const double mu = 0.01;  // Viscosity
const double dx = Lx / Nx;
const double dy = Ly / Ny;
const double dt = 0.01;  // Time step
const int maxIter = 1000;
const double tolerance = 1e-5;

// Field variables
std::vector<std::vector<double>> u(Nx+1, std::vector<double>(Ny, 0.0));  // X-velocity
std::vector<std::vector<double>> v(Nx+1, std::vector<double>(Ny, 0.0));  // Y-velocity
std::vector<std::vector<double>> p(Nx, std::vector<double>(Ny, 0.0));    // Pressure
std::vector<std::vector<double>> u_star(Nx+1, std::vector<double>(Ny, 0.0));
std::vector<std::vector<double>> v_star(Nx+1, std::vector<double>(Ny, 0.0));
std::vector<std::vector<double>> p_prime(Nx, std::vector<double>(Ny, 0.0));


void denseMatVec(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& y) ;
double inner_product(const std::vector<double>& a, const std::vector<double>& b) ;

// Solve Ax = b using the Conjugate Gradient method for a symmetric positive-definite dense matrix A.
bool conjugateGradientSolver(const std::vector<std::vector<double>>& A, const std::vector<double>& b,
                             std::vector<double>& x, int maxIter, double tolerance) {
    int n = b.size();
    std::vector<double> r(n), p(n), Ap(n);

    // Initialize x to zero and calculate the initial residual r = b - A * x (since x=0, r=b)
    //std::copy(b.begin(), b.end(), r);
    for (int k = 0 ; k < n; k ++)
      r[k] = b[k];
    p = r;
    double rs_old = inner_product(r, r);

    for (int k = 0; k < maxIter; ++k) {
        // Compute Ap = A * p
        denseMatVec(A, p, Ap);

        // Calculate alpha = rs_old / (p^T * Ap)
        double alpha = rs_old / inner_product(p, Ap);

        // Update x = x + alpha * p
        for (int i = 0; i < n; ++i)
            x[i] += alpha * p[i];

        // Update r = r - alpha * Ap
        for (int i = 0; i < n; ++i)
            r[i] -= alpha * Ap[i];

        // Check for convergence: if ||r|| < tolerance, break
        double rs_new = inner_product(r, r);
        if (std::sqrt(rs_new) < tolerance) {
            std::cout << "Conjugate Gradient converged after " << k + 1 << " iterations.\n";
            return true;
        }

        // Compute beta = rs_new / rs_old
        double beta = rs_new / rs_old;

        // Update p = r + beta * p
        for (int i = 0; i < n; ++i)
            p[i] = r[i] + beta * p[i];

        rs_old = rs_new;
    }

    std::cout << "Conjugate Gradient did not converge within the maximum number of iterations.\n";
    return false;
}

// Computes the matrix-vector product y = A * x for a dense matrix A.
void denseMatVec(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& y) {
    int n = A.size();
    y.assign(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            y[i] += A[i][j] * x[j];
        }
    }
}

// Helper function to calculate the inner product of two vectors
double inner_product(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}


//int main() {
//    // Example 3x3 symmetric positive-definite matrix
//    std::vector<std::vector<double>> A = {
//        {4, -1, 0},
//        {-1, 4, -1},
//        {0, -1, 3}
//    };
//
//    std::vector<double> b = {15, 10, 10};  // Right-hand side vector
//    std::vector<double> x(3, 0.0);         // Solution vector initialized to zero
//    int maxIter = 1000;
//    double tolerance = 1e-6;
//
//    // Solve Ax = b
//    if (conjugateGradientSolver(A, b, x, maxIter, tolerance)) {
//        std::cout << "Solution:\n";
//        for (double xi : x) std::cout << xi << " ";
//        std::cout << std::endl;
//    } else {
//        std::cout << "Solver failed to converge.\n";
//    }
//
//    return 0;
//}


void writeVTKFile(int timestep) {
    std::ofstream vtkFile;
    vtkFile.open("output_" + std::to_string(timestep) + ".vtk");
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
            double u_center = (u[i][j] + u[i+1][j]) / 2.0;
            double v_center = (v[i][j] + v[i][j+1]) / 2.0;
            vtkFile << u_center << " " << v_center << " 0.0\n";
        }
    }

    vtkFile.close();
    std::cout << "VTK file output_" << timestep << ".vtk written.\n";
}

void setBoundaryConditions() {
    // Set velocity boundary conditions (e.g., inlet/outlet)
    for (int j = 0; j < Ny; ++j) {
        u[0][j] = 1.0;  // Inlet velocity
        u[Nx][j] = u[Nx-1][j];  // Outlet boundary
    }
    for (int i = 0; i < Nx+1; ++i) {
        v[i][0] = 0.0;
        v[i][Ny] = 0.0;  // No-slip boundary at top and bottom
    }
}

void solveMomentumEquations() {
    for (int i = 1; i < Nx; ++i) {
        for (int j = 1; j < Ny-1; ++j) {
            // Discretize x-momentum equation
            u_star[i][j] = u[i][j] - dt * (
                (u[i+1][j]*u[i+1][j] - u[i-1][j]*u[i-1][j]) / (2*dx) +
                (u[i][j+1]*v[i][j+1] - u[i][j-1]*v[i][j-1]) / (2*dy) +
                (p[i][j] - p[i-1][j]) / dx / rho -
                mu * ((u[i+1][j] - 2*u[i][j] + u[i-1][j]) / (dx*dx) +
                      (u[i][j+1] - 2*u[i][j] + u[i][j-1]) / (dy*dy))
            );
            // Discretize y-momentum equation
            v_star[i][j] = v[i][j] - dt * (
                (u[i+1][j]*v[i+1][j] - u[i-1][j]*v[i-1][j]) / (2*dx) +
                (v[i][j+1]*v[i][j+1] - v[i][j-1]*v[i][j-1]) / (2*dy) +
                (p[i][j] - p[i][j-1]) / dy / rho -
                mu * ((v[i+1][j] - 2*v[i][j] + v[i-1][j]) / (dx*dx) +
                      (v[i][j+1] - 2*v[i][j] + v[i][j-1]) / (dy*dy))
            );
        }
    }
}

void solvePressureCorrectionEquation() {
    // This part is not correct, the pressure matrix is of size (nx*ny)*(nx*ny) instead of nx*ny
    std::vector<std::vector<double>> residual(Nx, std::vector<double>(Ny, 0.0));
    //std::vector<double> residual(Nx*Ny, 0.0);
    double alpha = 2 (1/(dx*dx) + 1/(dy*dy));
    double beta  = - 1/(dx*dx);
    std::vector<std::vector<double>> M(Nx*Ny, std::vector<double>(Nx*Ny, 0.0));
    std::vector<double> b (Nx*Ny,0);
    for (int i = 1; i < Nx-1; ++i) {
        for (int j = 1; j < Ny-1; ++j) {
            // Pressure correction equation
            double a_p = 1.0 / (dx*dx + dy*dy);
            residual[i][j] =- a_p * (
                (u_star[i][j] - u_star[i-1][j]) / dx +
                (v_star[i][j] - v_star[i][j-1]) / dy
            );
            M[i][j] = alpha;
            M[i+1][j] = beta;
            M[i-1][j] = beta;
            M[i][j+1] = beta;
            M[i][j-1] = beta;
            
//    std::vector<std::vector<double>> A = {
//        {4, -1, 0},
//        {-1, 4, -1},
//        {0, -1, 3}
//    };
//
//    std::vector<double> b = {15, 10, 10};  // Right-hand side vector
//    std::vector<double> x(3, 0.0);         // Solution vector initialized to zero
//    int maxIter = 1000;
//    double tolerance = 1e-6;
//
//    // Solve Ax = b
//    if (conjugateGradientSolver(A, b, x, maxIter, tolerance)) {
//        std::cout << "Solution:\n";
//        for (double xi : x) std::cout << xi << " ";
//        std::cout << std::endl;
//    } else {
//        std::cout << "Solver failed to converge.\n";
//    }
//
        }
    }
}

void correctVelocityAndPressure() {
    for (int i = 1; i < Nx; ++i) {
        for (int j = 1; j < Ny-1; ++j) {
            // Correct velocity and pressure
            u[i][j] = u_star[i][j] - dt * (p_prime[i][j] - p_prime[i-1][j]) / dx;
            v[i][j] = v_star[i][j] - dt * (p_prime[i][j] - p_prime[i][j-1]) / dy;
            p[i][j] += p_prime[i][j];
        }
    }
}

bool checkConvergence() {
    double maxError = 0.0;
    for (int i = 1; i < Nx-1; ++i) {
        for (int j = 1; j < Ny-1; ++j) {
            double continuity = (u[i+1][j] - u[i][j]) / dx + (v[i][j+1] - v[i][j]) / dy;
            maxError = std::max(maxError, std::abs(continuity));
        }
    }
    return maxError < tolerance;
}

void SIMPLESolver() {
		for (auto &row : p) {
			std::fill(row.begin(), row.end(), 0);
    }
    setBoundaryConditions();

    writeVTKFile(0);  // Write the solution to a VTK file every timestep
    for (int iter = 1; iter < maxIter; ++iter) {
        solveMomentumEquations();
        solvePressureCorrectionEquation();
        correctVelocityAndPressure();

        writeVTKFile(iter);  // Write the solution to a VTK file every timestep

        if (checkConvergence()) {
            std::cout << "Converged after " << iter << " iterations.\n";
            break;
        }
    }
}

int main() {
    SIMPLESolver();
    return 0;
}

