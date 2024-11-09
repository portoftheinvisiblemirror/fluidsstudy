#include <vector>
#include <cmath>
#include <iostream>

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


int main() {
    // Example 3x3 symmetric positive-definite matrix
    std::vector<std::vector<double>> A = {
        {4, -1, 0},
        {-1, 4, -1},
        {0, -1, 3}
    };

    std::vector<double> b = {15, 10, 10};  // Right-hand side vector
    std::vector<double> x(3, 0.0);         // Solution vector initialized to zero
    int maxIter = 1000;
    double tolerance = 1e-6;

    // Solve Ax = b
    if (conjugateGradientSolver(A, b, x, maxIter, tolerance)) {
        std::cout << "Solution:\n";
        for (double xi : x) std::cout << xi << " ";
        std::cout << std::endl;
    } else {
        std::cout << "Solver failed to converge.\n";
    }

    return 0;
}

