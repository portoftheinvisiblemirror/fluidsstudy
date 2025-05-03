#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <limits>
#include <algorithm>
#include "utilities.hpp"

// Demonstrates how to solve laplace equation numerically for a stick.

const int N = 40;  // Grid size
const double L = 10.0;  // Total length (-5 to 5)
const double dx = L / (N - 1);  // Grid spacing
const double V0 = 1.0;  // Central column potential
const double tolerance = 1e-6;  // Convergence tolerance
const int max_iterations = 1;
double dt;

// 3D vector to store potential values
std::vector<std::vector<std::vector<double>>> V(N,
    std::vector<std::vector<double>>(N,
        std::vector<double>(N, 0.0)));
//velocities 0=x, 1=y, 2=z  (4d)
std::vector<std::vector < std::vector <std::vector<double>>>> vel(3, 
    std::vector < std::vector <std::vector<double>>>(N, 
    std::vector <std::vector<double>>
        (N, std::vector<double>(N, 0))));

std::vector<std::vector < std::vector <std::vector<double>>>> veldiv(3,
    std::vector < std::vector <std::vector<double>>>(N,
        std::vector <std::vector<double>>
        (N, std::vector<double>(N, 0))));

template <typename T>
void statVector(const std::vector<std::vector<std::vector<double>>>& V, T & minVal, T & maxVal, T & avg){
    minVal = std::numeric_limits<T>::max();
    maxVal = std::numeric_limits<T>::min();

		T sum = 0;
		size_t N = 0;
		
    for (const auto& vec2D : V) {
        for (const auto& vec1D : vec2D) {
            for (T val : vec1D) {
                minVal = std::min(minVal, val);
                maxVal = std::max(maxVal, val);
								sum = sum + val;
								N ++;
            }
        }
    }
		if (N > 0) avg = sum/N;
		else avg = 0;
		
		return;
}

template <typename T>
void statVector(const std::vector<std::vector<std::vector<std::vector<double>>>>& V, T & minVal, T & maxVal, T & avg){
    minVal = std::numeric_limits<T>::max();
    maxVal = std::numeric_limits<T>::min();

		T sum = 0;
		size_t N = 0;
		
		for (const auto& vec3D: V){
			for (const auto& vec2D : vec3D) {
					for (const auto& vec1D : vec2D) {
							for (T val : vec1D) {
									minVal = std::min(minVal, val);
									maxVal = std::max(maxVal, val);
									sum = sum + val;
									N ++;
							}
					}
			}
		}
		if (N > 0) avg = sum/N;
		else avg = 0;
		
		return;
}


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
                double x = -L / 2.0 + i * dx;
                double y = -L / 2.0 + j * dx;
                double z = -L / 2.0 + k * dx;
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
std::vector<std::vector<std::vector<double>>> divselfadvect(std::vector< std::vector<std::vector<std::vector<double>>>>& v)
{
    std::vector<std::vector<std::vector<double>>> Q(N,
        std::vector<std::vector<double>>(N,
            std::vector<double>(N, 0.0)));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                //cycle through 3 directions:
                for (int l = 0; l < 3; ++l)
                    for(int m=0; m<3; ++m)
                    Q[i][j][k] += partdif(v, i, j, k, N, dx, l, m) * partdif(v, i, j, k, N, dx, m, l);
                
            }
        }
    }
    return Q;
}
std::vector<double> selfadvect(std::vector< std::vector<std::vector<std::vector<double>>>>& v, int x, int y, int z) //d=0,1,2
{
    std::vector <double> Q(3,0);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Q[i] += v[j][x][y][z] * partdif(v, x, y, z, N, dx, i, j);
        }
    }
    return Q;
}
std::vector<double> gradient(std::vector<std::vector<std::vector<double>>>& V, int x, int y ,int z)
{
    std::vector<double>Q(3, 0);
    for (int i = 0; i < 3; ++i)
    {
        Q[i] = partdif(V, x, y, z, N, dx, i);
    }
    return Q;
}
std::vector<double> laplacian(std::vector< std::vector<std::vector<std::vector<double>>>>& v, int x, int y, int z)
{
    std::vector <double> Q(3, 0);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Q[i] = partdif2(v, x, y, z, N, dx, i, j);
        }
    }
    return Q;
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

void timestepvelocity(std::vector< std::vector<std::vector<std::vector<double>>>>& v, std::vector<std::vector<std::vector<double>>>& P, const double h, const double Re)
{
    double vmax = findvmax(v, N);
    dt = 0.5 * h / vmax;
    std::vector<std::vector < std::vector <std::vector<double>>>> velnew(3,
        std::vector < std::vector <std::vector<double>>>(N,
            std::vector <std::vector<double>>
            (N, std::vector<double>(N, 0))));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                std::vector<double> vnv = selfadvect(v, i, j, k);
                std::vector<double> grad = gradient(P, i, j, k);
                std::vector<double> laplace = laplacian(v, i, j, k);
                for (int l = 0; l < 3; ++l)
                {
                    velnew[l][i][j][k] = vel[l][i][j][k] - dt * (vnv[l] - grad[l] + laplace[l] / Re);
                }
            }
        }
    }
    v = velnew;
}
void boundary(std::vector<std::vector<std::vector<double>>>& V, std::vector< std::vector<std::vector<std::vector<double>>>>& v, double v0) //water flows in through the x direction, walls are the planes parallel to xz and xy planes
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            //set inlet velocity
            v[0][0][i][j] = v0;
            v[1][0][i][j] = 0;
            v[2][0][i][j] = 0;

            //the velocity at a wall is zero
            v[0][i][0][j] = 0;
            v[1][i][0][j] = 0;
            v[2][i][0][j] = 0;
            v[0][i][N - 1][j] = 0;
            v[1][i][N - 1][j] = 0;
            v[2][i][N - 1][j] = 0;
            v[0][i][j][0] = 0;
            v[1][i][j][0] = 0;
            v[2][i][j][0] = 0;
            v[0][i][j][N - 1] = 0;
            v[1][i][j][N - 1] = 0;
            v[2][i][j][N - 1] = 0;

            //the pressure at an outlet is zero
            V[N - 1][i][j] = 0;
            //directional derivative of pressure against wall is zero
        }
    }
}

void solveGaussSeidel(std::vector<std::vector<std::vector<double>>>& V, std::vector< std::vector<std::vector<std::vector<double>>>>& v) {
    double maxDiff;
    int iter = 0;
    const double h2 = dx * dx;  // Grid spacing squared
    int mid = N / 2;
		double minVal, maxVal, avgVal;
    do {
        boundary(V, v,10);
				statVector(v, minVal, maxVal, avgVal);
				std::cout << "velocity minVal = " << minVal << " maxVal = " << maxVal << " avgVal = " << avgVal << std::endl;
        std::vector<std::vector<std::vector<double>>> A = divselfadvect(v);
				statVector(A, minVal, maxVal, avgVal);
				std::cout << "pdv_i v_j * pdv_j v_i minVal = " << minVal << " maxVal = " << maxVal << " avgVal = " << avgVal << std::endl;
        maxDiff = 0.0;
        // Interior points - 2nd order central difference
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                for (int k = 1; k < N - 1; k++) {
                    //// Skip the central charged column
                    //if (i == mid && j == mid && (k >= N / 4 && k < 3 * N / 4)) continue;

                    double oldV = V[i][j][k];
                    // Standard 7-point stencil for 3D Laplace equation
                    // This is based on 2nd order central difference
                    //    0 = (V[i+1] - 2V + V[i-1])/dx² +
                    //     (V[j+1] - 2V + V[j-1])/dx² +
                    //     (V[k+1] - 2V + V[k-1])/dx²
                    V[i][j][k] = (V[i + 1][j][k] + V[i - 1][j][k] +
                        V[i][j + 1][k] + V[i][j - 1][k] +
                        V[i][j][k + 1] + V[i][j][k - 1]+dx*dx*A[i][j][k]) / 6.0;

                    maxDiff = std::max(maxDiff, std::abs(V[i][j][k] - oldV));
                }
            }
        }
        printf("at the end of iteration: %f\n ", maxDiff);

        // Boundary points - 2nd order one-sided differences
        // x boundaries (i = 0 and i = N-1)
        //    V₁ = V₀ + h(∂V/∂x) + (h²/2)(∂²V/∂x²) + (h³/6)(∂³V/∂x³) + O(h⁴)
        //    V₂ = V₀ + 2h(∂V/∂x) + 2h²(∂²V/∂x²) + (8h³/6)(∂³V/∂x³) + O(h⁴)
        //    V₃ = V₀ + 3h(∂V/∂x) + (9h²/2)(∂²V/∂x²) + (27h³/6)(∂³V/∂x³) + O(h⁴)
        //    From the above difference equations, we solve for V0 at the boundary points
        //    V0 = (-)
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                // Left boundary (i = 0) - forward difference
                double oldV = V[0][j][k];
                V[0][j][k] = (5 * V[1][j][k] - 4 * V[2][j][k] + V[3][j][k] + dx * dx * A[0][j][k]) / 2.0;  // 2nd order forward difference
                maxDiff = std::max(maxDiff, std::abs(V[0][j][k] - oldV));

                // Right boundary (i = N-1) - backward difference
                oldV = V[N - 1][j][k];
                V[N - 1][j][k] = (5 * V[N - 2][j][k] - 4 * V[N - 3][j][k] + V[N - 4][j][k] + dx * dx * A[N-1][j][k]) / 2.0;
                maxDiff = std::max(maxDiff, std::abs(V[N - 1][j][k] - oldV));
            }
        }
        printf("at the end of iteration: %f\n ", maxDiff);
        // y boundaries (j = 0 and j = N-1)
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < N; k++) {
                // Front boundary (j = 0) - forward difference
                double oldV = V[i][0][k];
                V[i][0][k] = (5 * V[i][1][k] - 4 * V[i][2][k] + V[i][3][k] + dx * dx * A[i][0][k]) / 2.0;
                maxDiff = std::max(maxDiff, std::abs(V[i][0][k] - oldV));

                // Back boundary (j = N-1) - backward difference
                oldV = V[i][N - 1][k];
                V[i][N - 1][k] = (5 * V[i][N - 2][k] - 4 * V[i][N - 3][k] + V[i][N - 4][k] + dx * dx * A[i][N-1][k]) / 2.0;
                maxDiff = std::max(maxDiff, std::abs(V[i][N - 1][k] - oldV));
            }
        }
        printf("at the end of iteration: %f\n ", maxDiff);
        // z boundaries (k = 0 and k = N-1)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                // Bottom boundary (k = 0) - forward difference
                double oldV = V[i][j][0];
                V[i][j][0] = (5 * V[i][j][1] - 4 * V[i][j][2] + V[i][j][3] + dx * dx * A[i][j][0]) / 2.0;
                maxDiff = std::max(maxDiff, std::abs(V[i][j][0] - oldV));

                // Top boundary (k = N-1) - backward difference
                oldV = V[i][j][N - 1];
                V[i][j][N - 1] = (5 * V[i][j][N - 2] - 4 * V[i][j][N - 3] + V[i][j][N - 4] + dx * dx * A[i][j][N-1]) / 2.0;
                maxDiff = std::max(maxDiff, std::abs(V[i][j][N - 1] - oldV));
            }
        }
        printf("at the end of iteration: %f\n ", maxDiff);
        timestepvelocity(vel, V, dx, 10);
        iter++;

        if (iter % 100 == 0) {
            std::cout << "Iteration " << iter << ", max difference = " << maxDiff << std::endl;
        }
        printf("at the end of iteration: %f\n ", maxDiff);
    } while (maxDiff > tolerance && iter < max_iterations);

    std::cout << "Converged after " << iter << " iterations" << std::endl;
}

int main() {
    double Re = 10;
    // Solve using Gauss-Seidel method
    solveGaussSeidel(V,vel);

    // Write output to VTK file
    writeVTK(V, "laplace3d.vtk");

    return 0;
}
