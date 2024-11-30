#ifndef NS3D_HPP
#define NS3D_HPP

#include <vector>

class NavierStokesSolver3D {
private:
    // Grid parameters
    int nx, ny, nz;
    double dx, dy, dz, dt;
    double rho, mu;
    
    // Flow field variables
    std::vector<std::vector<std::vector<double>>> u; // x-velocity
    std::vector<std::vector<std::vector<double>>> v; // y-velocity
    std::vector<std::vector<std::vector<double>>> w; // z-velocity
    std::vector<std::vector<std::vector<double>>> p; // pressure
    
    // Temporary arrays for calculations
    std::vector<std::vector<std::vector<double>>> u_temp;
    std::vector<std::vector<std::vector<double>>> v_temp;
    std::vector<std::vector<std::vector<double>>> w_temp;

    // Helper methods
    void solvePoisson();
    void applyBoundaryConditions();
    
public:
    NavierStokesSolver3D(int nx, int ny, int nz, double dx, double dy, double dz, 
                        double dt, double rho, double mu);
    
    void step();
    
    // Getters for flow field variables
    const std::vector<std::vector<std::vector<double>>>& getU() const { return u; }
    const std::vector<std::vector<std::vector<double>>>& getV() const { return v; }
    const std::vector<std::vector<std::vector<double>>>& getW() const { return w; }
    const std::vector<std::vector<std::vector<double>>>& getP() const { return p; }
};

#endif 