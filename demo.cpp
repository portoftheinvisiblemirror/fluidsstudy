#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// Demonstrates how to solve laplace equation numerically for a stick.

const int N = 40;  // Grid size
const double L = 10.0;  // Total length (-5 to 5)
const double dx = L / (N-1);  // Grid spacing
const double V0 = 1.0;  // Central column potential
const double tolerance = 1e-6;  // Convergence tolerance
const int max_iterations = 10000;

// 3D vector to store pressure values and velocity values
std::vector<std::vector<std::vector<double>>> V(N, 
    std::vector<std::vector<double>>(N,
        std::vector<double>(N, 0.0)));
std::vector<std::vector<std::vector<std::vector<double>>>> vel(N, 
    std::vector<std::vector<std::vector<double>>>(N,
    std::vector<std::vector<double>>(N,
        std::vector<double>(N, 0.0))));

void writeVTK(const std::vector<std::vector<std::vector<double>>>& V, const std::string& filename) {
    std::ofstream vtkFile(filename);
    
    // VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "3D Laplace Solution\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << N << " " << N << " " << N << "\n";
    
    // Points
    vtkFile << "POINTS " << N * N * N << " float\n";
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                double x = -L/2.0 + i * dx;
                double y = -L/2.0 + j * dx;
                double z = -L/2.0 + k * dx;
                vtkFile << x << " " << y << " " << z << "\n";
            }
        }
    }
    
    // Point data (potential values)
    vtkFile << "POINT_DATA " << N * N * N << "\n";
    vtkFile << "SCALARS potential float\n";
    vtkFile << "LOOKUP_TABLE default\n";
    
    // Write potential values
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                vtkFile << V[i][j][k] << "\n";
            }
        }
    }
    
    vtkFile.close();
}

// Solve h^2 (partial_i v_j * \partial_j v_i)
double solveB(const std::vector<std::vector<std::vector<std::vector<double>>>>& vel, 
            const int i, const int j, const int k, const double h) {
  double B = 0.0;
	//assert( i > 0);
	//assert( j > 0);
	//assert( k > 0);
	//assert( i < N-1);
	//assert( j < N-1);
	//assert( k < N-1);
  // for cells within the boundary, 2nd order accuracy central difference
  if ( i > 0 && i < N-1 && j > 0 && j < N -1 && k > 0 && k < N - 1){ 
            for (int m = 0; m < 2 ; m ++){ // 0: x, 1: y, 2: z
							for (int n = 0; n < 2 ; n ++){
								if(m == n && m == 0) {
 									B += pow((0.5*(vel[i+1][j][k][m] - vel[i-1][j][k][m])/h), 2); // (partial_x v_x)^2
                  // alternatively, define a pfxpx operator, then say B+=(pfxpx(...))^2
								}
								if(m == n && m == 1) {
 									B += pow((0.5*(vel[i][j+1][k][m] - vel[i][j-1][k][m])/h), 2); // (partial_x v_x)^2
								}
								if(m == n && m == 2) {
 									B += pow((0.5*(vel[i][j][k+1][m] - vel[i][j][k-1][m])/h), 2); // (partial_x v_x)^2
								}
								if(m == 0 && n == 1) {
 									B += (0.5*(vel[i+1][j][k][m] - vel[i-1][j][k][m])/h)*( 0.5*(vel[i][j+1][k][n] - vel[i][j-1][k][n])/h);
								}
								if(m == 0 && n == 2) {
 									B += (0.5*(vel[i+1][j][k][m] - vel[i-1][j][k][m])/h)*( 0.5*(vel[i][j][k+1][n] - vel[i][j][k-1][n])/h);
								}
								if(m == 1 && n == 0) {
 									B += (0.5*(vel[i][j+1][k][m] - vel[i][j-1][k][m])/h)*( 0.5*(vel[i+1][j][k][n] - vel[i-1][j][k][n])/h);
								}
								if(m == 1 && n == 2) {
 									B += (0.5*(vel[i][j+1][k][m] - vel[i][j-1][k][m])/h)*( 0.5*(vel[i][j][k+1][n] - vel[i][j][k-1][n])/h);
								}
								if(m == 2 && n == 0) {
 									B += (0.5*(vel[i][j][k+1][m] - vel[i][j][k-1][m])/h)*( 0.5*(vel[i+1][j][k][n] - vel[i-1][j][k][n])/h);
								}
								if(m == 2 && n == 1) {
 									B += (0.5*(vel[i][j][k+1][m] - vel[i][j][k-1][m])/h)*( 0.5*(vel[i][j+1][k][n] - vel[i][j-1][k][n])/h);
								}
							}
						}
	} // inner cells
	B *= h*h;
	return B;
}

void solveGaussSeidel(std::vector<std::vector<std::vector<double>>>& V) {
    double maxDiff;
    int iter = 0;
    const double h2 = dx * dx;  // Grid spacing squared
    int mid = N/2;

    do {
        maxDiff = 0.0;
        // Interior points - 2nd order central difference
        for(int i = 1; i < N-1; i++) {
            for(int j = 1; j < N-1; j++) {
                for(int k = 1; k < N-1; k++) {
                    // Skip the central charged column
                    if(i == mid && j == mid && (k >= N/4 && k < 3*N/4)) continue;
                    
                    double oldV = V[i][j][k];
                    // Standard 7-point stencil for 3D Laplace equation
                    // This is based on 2nd order central difference
                    //    0 = (V[i+1] - 2V + V[i-1])/dx² +
                    //     (V[j+1] - 2V + V[j-1])/dx² +
                    //     (V[k+1] - 2V + V[k-1])/dx²
                    V[i][j][k] = (V[i+1][j][k] + V[i-1][j][k] + 
                                  V[i][j+1][k] + V[i][j-1][k] +
                                  V[i][j][k+1] + V[i][j][k-1] - solveB(vel, i,j,k, dx) ) / 6.0;
                    
                    maxDiff = std::max(maxDiff, std::abs(V[i][j][k] - oldV));
                }
            }
        }

        // Boundary points - 2nd order one-sided differences
        // x boundaries (i = 0 and i = N-1)
        //    V₁ = V₀ + h(∂V/∂x) + (h²/2)(∂²V/∂x²) + (h³/6)(∂³V/∂x³) + O(h⁴)
        //    V₂ = V₀ + 2h(∂V/∂x) + 2h²(∂²V/∂x²) + (8h³/6)(∂³V/∂x³) + O(h⁴)
        //    V₃ = V₀ + 3h(∂V/∂x) + (9h²/2)(∂²V/∂x²) + (27h³/6)(∂³V/∂x³) + O(h⁴)
        //    From the above difference equations, we solve for V0 at the boundary points
        //    V0 = (-)
        for(int j = 0; j < N; j++) {
            for(int k = 0; k < N; k++) {
                // Left boundary (i = 0) - forward difference
                double oldV = V[0][j][k];
                V[0][j][k] = (5*V[1][j][k] - 4*V[2][j][k] + V[3][j][k])/2.0;  // 2nd order forward difference
                maxDiff = std::max(maxDiff, std::abs(V[0][j][k] - oldV));

                // Right boundary (i = N-1) - backward difference
                oldV = V[N-1][j][k];
                V[N-1][j][k] = (5*V[N-2][j][k] - 4*V[N-3][j][k] + V[N-4][j][k])/2.0;
                maxDiff = std::max(maxDiff, std::abs(V[N-1][j][k] - oldV));
            }
        }

        // y boundaries (j = 0 and j = N-1)
        for(int i = 0; i < N; i++) {
            for(int k = 0; k < N; k++) {
                // Front boundary (j = 0) - forward difference
                double oldV = V[i][0][k];
                V[i][0][k] = (5*V[i][1][k] - 4*V[i][2][k] + V[i][3][k])/2.0;
                maxDiff = std::max(maxDiff, std::abs(V[i][0][k] - oldV));

                // Back boundary (j = N-1) - backward difference
                oldV = V[i][N-1][k];
                V[i][N-1][k] = (5*V[i][N-2][k] - 4*V[i][N-3][k] + V[i][N-4][k])/2.0;
                maxDiff = std::max(maxDiff, std::abs(V[i][N-1][k] - oldV));
            }
        }

        // z boundaries (k = 0 and k = N-1)
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                // Bottom boundary (k = 0) - forward difference
                double oldV = V[i][j][0];
                V[i][j][0] = (5*V[i][j][1] - 4*V[i][j][2] + V[i][j][3])/2.0;
                maxDiff = std::max(maxDiff, std::abs(V[i][j][0] - oldV));

                // Top boundary (k = N-1) - backward difference
                oldV = V[i][j][N-1];
                V[i][j][N-1] = (5*V[i][j][N-2] - 4*V[i][j][N-3] + V[i][j][N-4])/2.0;
                maxDiff = std::max(maxDiff, std::abs(V[i][j][N-1] - oldV));
            }
        }
        // x,y plane boundary intersecting edge x 2
        // y,z plane boundary intersecting edge x 2
        // x,z plane boundary intersecting edge x 2
        // 8 x,y,z plane boundary intersecting corner

        iter++;
        
        if(iter % 100 == 0) {
            std::cout << "Iteration " << iter << ", max difference = " << maxDiff << std::endl;
        }
    } while(maxDiff > tolerance && iter < max_iterations);

    std::cout << "Converged after " << iter << " iterations" << std::endl;
}

double computeLinearDensity(const std::vector<std::vector<std::vector<double>>>& V) {
    double eps0 = 8.85e-12;  // Permittivity of free space
    std::vector<double> rho_linear(N/2, 0.0);
    double total_rho = 0.0;
    int mid = N/2;
    
    // Open output file for linear charge density distribution
    std::ofstream rhoFile("rho_linear.txt");
    rhoFile << "# k-position(m)  rho_linear(C/m)\n";
    
    // Loop over the central column
    for(int k = N/4; k < 3*N/4; k++) {
        // Compute Laplacian using 2nd order central differences
        double d2x = (V[mid+1][mid][k] - 2*V[mid][mid][k] + V[mid-1][mid][k]) / (dx*dx);
        double d2y = (V[mid][mid+1][k] - 2*V[mid][mid][k] + V[mid][mid-1][k]) / (dx*dx);
        double d2z = (V[mid][mid][k+1] - 2*V[mid][mid][k] + V[mid][mid][k-1]) / (dx*dx);
        
        double laplacian = d2x + d2y + d2z;
        double rho = -eps0 * laplacian;  // Charge density at this point
        rho_linear[k-N/4] = rho * dx;  // Integrate over the length element
        total_rho += rho_linear[k-N/4];
        
        // Write position and charge density to file
    }
    double total_q = total_rho * dx;
    for (int k = N/4; k < 3*N/4; k++) {
        double z_position = -L/2.0 + k * dx;  // Physical position along z-axis
        double theoretical_rho = total_q /( M_PI * sqrt(L*L - 4*z_position*z_position));
        //rhoFile << z_position << "\t" << rho_linear[k-N/4] << "\t" << theoretical_rho << "\n";
        rhoFile << z_position << "\t" << rho_linear[k-N/4] << "\n";
    }
    
    rhoFile.close();
    std::cout << "Total linear charge density: " << total_rho << " C/m" << std::endl;
    
    return total_rho;
}

int main() {
    // Initialize V to zero
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            for(int k = 0; k < N; k++) {
                V[i][j][k] = 0.0;
            }
        }
    }

    // Set central column to V0
    int mid = N/2;
    for(int k = N/4; k < 3*N/4; k++) {
        V[mid][mid][k] = V0;
    }

    // Solve using Gauss-Seidel method
    solveGaussSeidel(V);

    // Compute and output the linear charge density
    double total_rho = computeLinearDensity(V);

    // Write output to VTK file
    writeVTK(V, "laplace3d.vtk");

    return 0;
}
