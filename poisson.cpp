#include <vector>
#include <iostream>

const double dx = 1.0;  // Grid spacing in x direction
const double dy = 1.0;  // Grid spacing in y direction
const double dt = 0.01; // Time step

// Function to compute the divergence of the velocity field
void computeDivergence(const std::vector<std::vector<double>>& u, const std::vector<std::vector<double>>& v,
                       std::vector<std::vector<double>>& divergence) {
    int nx = u.size();
    int ny = u[0].size();

    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            divergence[i][j] = (u[i+1][j] - u[i][j]) / dx + (v[i][j+1] - v[i][j]) / dy;
        }
    }
}

// Function to build the coefficient matrix and the RHS vector for the Poisson equation
void buildPoissonSystem(const std::vector<std::vector<double>>& divergence,
                        std::vector<std::vector<double>>& A, std::vector<double>& b, int nx, int ny) {
    int n = nx * ny;
    A.resize(n, std::vector<double>(n, 0.0));
    b.resize(n, 0.0);

    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            int idx = i * ny + j;

            A[idx][idx] = -4.0; // Center coefficient

            // Neighbor coefficients
            A[idx][idx + 1] = 1.0;    // Right neighbor
            A[idx][idx - 1] = 1.0;    // Left neighbor
            A[idx][idx + ny] = 1.0;   // Top neighbor
            A[idx][idx - ny] = 1.0;   // Bottom neighbor

            // Right-hand side
            b[idx] = -divergence[i][j] * dx * dx / dt; // Eqn (30), (36), (40)
        }
    }
}

// Example usage
int main() {
    int nx = 5; // Grid points in x direction
    int ny = 5; // Grid points in y direction

    // Initialize velocity fields u and v
    std::vector<std::vector<double>> u(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double>> v(nx, std::vector<double>(ny, 0.0));

    // Initialize divergence field
    std::vector<std::vector<double>> divergence(nx, std::vector<double>(ny, 0.0));

    // Compute divergence of velocity field
    computeDivergence(u, v, divergence);

    // Set up Poisson system (A * p = b)
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    buildPoissonSystem(divergence, A, b, nx, ny);

    // Display A and b for debugging
    std::cout << "Coefficient matrix A:" << std::endl;
    for (const auto& row : A) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Right-hand side vector b:" << std::endl;
    for (double val : b) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}

