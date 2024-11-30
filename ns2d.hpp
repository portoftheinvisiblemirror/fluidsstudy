#ifndef NS2D_HPP
#define NS2D_HPP

#include <vector>

class NavierStokesSolver2D {
private:
    int nx, ny;                     // Grid dimensions
    double dx, dy, dt;              // Space and time steps
    double rho, mu;                 // Density and viscosity
    std::vector<std::vector<double>> u, v;  // Velocity components
    std::vector<std::vector<double>> p;     // Pressure
    std::vector<std::vector<double>> tmp;   // Temporary array

    void setBoundaryConditions();
    void solvePoisson();

public:
    NavierStokesSolver2D(int nx_, int ny_, double dx_, double dy_, double dt_,
                        double rho_=1.0, double mu_=0.1);

    void step();

    // Getters for visualization
    const std::vector<std::vector<double>>& getU() const { return u; }
    const std::vector<std::vector<double>>& getV() const { return v; }
    const std::vector<std::vector<double>>& getP() const { return p; }
};

#endif 