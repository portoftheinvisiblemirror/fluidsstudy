
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>


void solve_momentum_equations();
void solve_pressure_correction();
void correct_pressure_velocity();
void check_convergence(double& max_div, double& max_p_change);
void apply_boundary_conditions();
void output_results(int time_step, double time);
void output_final_results();
void deallocate_arrays();

using namespace std;
typedef double real;
    //constants
    size_t nx = 32, ny = 32, nz = 64;   // Smaller grid for testing
    size_t max_iter = 6000;             // Maximum time steps
    size_t n_correctors = 2;            // Fewer correctors
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
    real beta = 0.3;                // Under-relaxation for pressure
    
    // Arrays
    std::vector<std::vector<std::vector<double>>> u, v, w;         // Velocity components
    std::vector<std::vector<std::vector<double>>> p;               // Pressure
    std::vector<std::vector<std::vector<double>>> u_old, v_old, w_old; // Old velocities
    std::vector<std::vector<std::vector<double>>> p_old;           // Old pressure
    std::vector<std::vector<std::vector<double>>> u_star, v_star, w_star; // Intermediate velocities
    std::vector<std::vector<std::vector<double>>> p_prime;         // Pressure correction
    std::vector<std::vector<std::vector<double>>> div;                           // Divergence
    std::vector<double> x, y, z;                                   // Grid coordinates
    std::vector<std::vector<double>> radius, theta_coord;         // Cylindrical coordinates
    
int main()
{
    // Variables
    size_t  i, j, k, iter, time_step, corrector;
    real  time, max_div, max_p_change;
    real  x_center, y_center, r_dist, inlet_velocity;
    string filename;
    
    cout << "Initializing 3D Navier - Stokes PISO Solver(Stable Version)\n";
    std::cout << "Grid: " << nx << " x " << ny << " x " << nz << " dt: " << dt << " Re: " << Re << std::endl;
    
    // Allocate arrays
    allocate_arrays();
    
    // Initialize grid
    initialize_grid();
    
    // Initialize flow field
    initialize_flow();
    

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
                double max_div, max_p_change; // Initialize these variables
                check_convergence(max_div, max_p_change);

                if (time_step % 10 == 0) {
                    std::cout << "  Corrector: " << corrector << " Max div: " << max_div << " Max P change: " << max_p_change << std::endl;
                }

                // If converged, exit corrector loop
                if (max_div < tolerance && max_p_change < tolerance) {
                    break;
                }
            }

            // Apply boundary conditions
            apply_boundary_conditions();

            // Output results periodically
            if (time_step % 50 == 0) {
                output_results(time_step, time);
            }

            // Check for instability
            if (max_div > 10.0 || max_p_change > 10.0) {
                std::cout << "WARNING: Simulation becoming unstable" << std::endl;
                std::cout << "Max div: " << max_div << " Max P change: " << max_p_change << std::endl;
                std::cout << "Stopping at time step: " << time_step << std::endl;
                break;
            }
        }

        // Final output
        output_final_results();

        // Deallocate arrays
        deallocate_arrays();

        std::cout << "Simulation completed successfully" << std::endl;

        return 0;
    }


    void allocate_arrays()
    {
        u.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    v.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    w.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    u_old.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    v_old.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    w_old.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    u_star.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    v_star.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    w_star.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    p.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    p_old.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    p_prime.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    div.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));

    x.resize(nx, 0.0);
    y.resize(ny, 0.0);
    z.resize(nz, 0.0);

    radius.resize(ny, std::vector<std::vector<double>>(nz, 0.0));
    theta_coord.resize(ny, std::vector<std::vector<double>>(nz, 0.0));
}
    
    
    void initialize_grid()
    {
        // Initialize Cartesian grid coordinates
        for (size_t i = 1; i <= nx; ++i)
            x[i] = (i - 1) * dx;
        for (size_t i = 1; i <= ny; ++i)
            y[i] = -R+(i - 1) * dy;
        for (size_t i = 1; i <= nz; ++i)
            z[i] = -R+(i - 1) * dz;

                // Calculate cylindrical coordinates
        double x_center = 0.0;
        double y_center = 0.0;
        for (size_t j = 1; j <= ny; ++j)
        {
            for (size_t k = 1; k <= nz; ++k)
            {
                double r_dist = sqrt((y[j] - y_center) * (y[j] - y_center) + (z[k] - y_center) * (z[k] - y_center));
                radius[j][k] = r_dist;
                theta_coord[j][k] = atan2(z[k] - y_center, y[j] - y_center);
            }
        }
        cout << "Grid initialized : nx = " << nx << "ny = " << ny << "nz = " << nz;
        cout << "Grid spacing: dx=" << dx << " dy=" << dy << " dz =" << dz;
    }
    /*
    
    
    subroutine initialize_flow()
        // Initialize velocity and pressure fields
        u = 0.0
        v = 0.0
        w = 0.0
        p = 0.0
        
        // Set inlet velocity profile (parabolic)
        inlet_velocity = 0.5  // Reduced inlet velocity
        do j = 1, ny
            do k = 1, nz
                if (radius(j,k) <= R) then
                    u(1,j,k) = inlet_velocity * (1.0 - (radius(j,k)/R)**2)
                else
                    u(1,j,k) = 0.0
                end if
            end do
        end do
        
        u_old = u
        v_old = v
        w_old = w
        p_old = p
        
        write(*,*) 'Flow field initialized with parabolic inlet profile'
        write(*,*) 'Inlet velocity:', inlet_velocity, 'Reynolds number:', Re
    end subroutine initialize_flow
    
    subroutine solve_momentum_equations()
        // Local variables
        size_t :: i, j, k
        real :: conv_x, conv_y, conv_z
        real :: diff_x, diff_y, diff_z
        real :: dp_dx, dp_dy, dp_dz
        
        // Solve u-momentum equation
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if (radius(j,k) <= R) then  // Only solve inside tube
                        // Convection terms (central differencing for stability)
                        conv_x = u(i,j,k) * (u(i+1,j,k) - u(i-1,j,k)) / (2.0*dx)
                        conv_y = v(i,j,k) * (u(i,j+1,k) - u(i,j-1,k)) / (2.0*dy)
                        conv_z = w(i,j,k) * (u(i,j,k+1) - u(i,j,k-1)) / (2.0*dz)
                        
                        // Diffusion terms
                        diff_x = nu * (u(i+1,j,k) - 2.0*u(i,j,k) + u(i-1,j,k)) / (dx**2)
                        diff_y = nu * (u(i,j+1,k) - 2.0*u(i,j,k) + u(i,j-1,k)) / (dy**2)
                        diff_z = nu * (u(i,j,k+1) - 2.0*u(i,j,k) + u(i,j,k-1)) / (dz**2)
                        
                        // Pressure gradient
                        dp_dx = (p(i+1,j,k) - p(i-1,j,k)) / (2.0*dx*rho)
                        
                        // Update intermediate velocity with under-relaxation
                        u_star(i,j,k) = u_old(i,j,k) + dt * (-conv_x - conv_y - conv_z + diff_x + diff_y + diff_z - dp_dx)
                        u_star(i,j,k) = alpha * u_star(i,j,k) + (1.0 - alpha) * u(i,j,k)
                    end if
                end do
            end do
        end do
        
        // Solve v-momentum equation
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if (radius(j,k) <= R) then
                        // Convection terms (central differencing)
                        conv_x = u(i,j,k) * (v(i+1,j,k) - v(i-1,j,k)) / (2.0*dx)
                        conv_y = v(i,j,k) * (v(i,j+1,k) - v(i,j-1,k)) / (2.0*dy)
                        conv_z = w(i,j,k) * (v(i,j,k+1) - v(i,j,k-1)) / (2.0*dz)
                        
                        // Diffusion terms
                        diff_x = nu * (v(i+1,j,k) - 2.0*v(i,j,k) + v(i-1,j,k)) / (dx**2)
                        diff_y = nu * (v(i,j+1,k) - 2.0*v(i,j,k) + v(i,j-1,k)) / (dy**2)
                        diff_z = nu * (v(i,j,k+1) - 2.0*v(i,j,k) + v(i,j,k-1)) / (dz**2)
                        
                        // Pressure gradient
                        dp_dy = (p(i,j+1,k) - p(i,j-1,k)) / (2.0*dy*rho)
                        
                        // Update intermediate velocity with under-relaxation
                        v_star(i,j,k) = v_old(i,j,k) + dt * (-conv_x - conv_y - conv_z + diff_x + diff_y + diff_z - dp_dy)
                        v_star(i,j,k) = alpha * v_star(i,j,k) + (1.0 - alpha) * v(i,j,k)
                    end if
                end do
            end do
        end do
        
        // Solve w-momentum equation
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if (radius(j,k) <= R) then
                        // Convection terms (central differencing)
                        conv_x = u(i,j,k) * (w(i+1,j,k) - w(i-1,j,k)) / (2.0*dx)
                        conv_y = v(i,j,k) * (w(i,j+1,k) - w(i,j-1,k)) / (2.0*dy)
                        conv_z = w(i,j,k) * (w(i,j,k+1) - w(i,j,k-1)) / (2.0*dz)
                        
                        // Diffusion terms
                        diff_x = nu * (w(i+1,j,k) - 2.0*w(i,j,k) + w(i-1,j,k)) / (dx**2)
                        diff_y = nu * (w(i,j+1,k) - 2.0*w(i,j,k) + w(i,j-1,k)) / (dy**2)
                        diff_z = nu * (w(i,j,k+1) - 2.0*w(i,j,k) + w(i,j,k-1)) / (dz**2)
                        
                        // Pressure gradient
                        dp_dz = (p(i,j,k+1) - p(i,j,k-1)) / (2.0*dz*rho)
                        
                        // Update intermediate velocity with under-relaxation
                        w_star(i,j,k) = w_old(i,j,k) + dt * (-conv_x - conv_y - conv_z + diff_x + diff_y + diff_z - dp_dz)
                        w_star(i,j,k) = alpha * w_star(i,j,k) + (1.0 - alpha) * w(i,j,k)
                    end if
                end do
            end do
        end do
    end subroutine solve_momentum_equations
    
    subroutine solve_pressure_correction()
        // Local variables
        size_t :: i, j, k, sor_iter
        real :: omega, p_new
        
        omega = 1.2  // Reduced SOR relaxation factor
        
        // Calculate divergence of intermediate velocity field
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if (radius(j,k) <= R) then
                        div(i,j,k) = (u_star(i+1,j,k) - u_star(i-1,j,k)) / (2.0*dx) + &
                                    (v_star(i,j+1,k) - v_star(i,j-1,k)) / (2.0*dy) + &
                                    (w_star(i,j,k+1) - w_star(i,j,k-1)) / (2.0*dz)
                    else
                        div(i,j,k) = 0.0
                    end if
                end do
            end do
        end do
        
        // Solve pressure correction equation using SOR
        do sor_iter = 1, 10  // Fewer iterations
            do i = 2, nx-1
                do j = 2, ny-1
                    do k = 2, nz-1
                        if (radius(j,k) <= R) then
                            p_new = (1.0/6.0) * (p_prime(i+1,j,k) + p_prime(i-1,j,k) + &
                                                p_prime(i,j+1,k) + p_prime(i,j-1,k) + &
                                                p_prime(i,j,k+1) + p_prime(i,j,k-1) - &
                                                (dx**2) * div(i,j,k))
                            
                            p_prime(i,j,k) = p_prime(i,j,k) + omega * (p_new - p_prime(i,j,k))
                        end if
                    end do
                end do
            end do
        end do
    end subroutine solve_pressure_correction
    
    subroutine correct_pressure_velocity()
        // Local variables
        size_t :: i, j, k
        
        // Correct pressure with under-relaxation
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if (radius(j,k) <= R) then
                        p(i,j,k) = p(i,j,k) + beta * p_prime(i,j,k)
                    end if
                end do
            end do
        end do
        
        // Correct velocities
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if (radius(j,k) <= R) then
                        u(i,j,k) = u_star(i,j,k) - dt * (p_prime(i+1,j,k) - p_prime(i-1,j,k)) / (2.0*dx*rho)
                        v(i,j,k) = v_star(i,j,k) - dt * (p_prime(i,j+1,k) - p_prime(i,j-1,k)) / (2.0*dy*rho)
                        w(i,j,k) = w_star(i,j,k) - dt * (p_prime(i,j,k+1) - p_prime(i,j,k-1)) / (2.0*dz*rho)
                    end if
                end do
            end do
        end do
    end subroutine correct_pressure_velocity
    
    subroutine check_convergence(max_div, max_p_change)
        real, intent(out) :: max_div, max_p_change
        
        // Local variables
        size_t :: i, j, k
        
        max_div = 0.0
        max_p_change = 0.0
        
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if (radius(j,k) <= R) then
                        max_div = max(max_div, abs(div(i,j,k)))
                        max_p_change = max(max_p_change, abs(p_prime(i,j,k)))
                    end if
                end do
            end do
        end do
    end subroutine check_convergence
    
    subroutine apply_boundary_conditions()
        // Local variables
        size_t :: i, j, k
        
        // Inlet boundary (parabolic velocity profile)
        inlet_velocity = 0.5  // Reduced inlet velocity
        do j = 1, ny
            do k = 1, nz
                if (radius(j,k) <= R) then
                    u(1,j,k) = inlet_velocity * (1.0 - (radius(j,k)/R)**2)
                    v(1,j,k) = 0.0
                    w(1,j,k) = 0.0
                else
                    u(1,j,k) = 0.0
                    v(1,j,k) = 0.0
                    w(1,j,k) = 0.0
                end if
            end do
        end do
        
        // Outlet boundary (zero gradient)
        do j = 1, ny
            do k = 1, nz
                u(nx,j,k) = u(nx-1,j,k)
                v(nx,j,k) = v(nx-1,j,k)
                w(nx,j,k) = w(nx-1,j,k)
                p(nx,j,k) = p(nx-1,j,k)
            end do
        end do
        
        // Wall boundary (no-slip)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    if (radius(j,k) > R) then
                        u(i,j,k) = 0.0
                        v(i,j,k) = 0.0
                        w(i,j,k) = 0.0
                    end if
                end do
            end do
        end do
        
        // Pressure boundary conditions
        p(1,:,:) = p(2,:,:)  // Inlet
        p(nx,:,:) = 0.0      // Outlet (reference pressure)
    end subroutine apply_boundary_conditions
    
    subroutine output_results(time_step, time)
        size_t, intent(in) :: time_step
        real, intent(in) :: time
        
        // Local variables
        size_t :: i, j, k
        
        // Output velocity field
        write(filename, '(A,I4.4,A)') 'stable_piso_velocity_', time_step, '.vtk'
        open(unit=10, file=filename, status='replace')
        
        write(10,'(a)') '# vtk DataFile Version 3.0'
        write(10,'(a)') '3D Navier-Stokes PISO Flow in Circular Tube (Stable)'
        write(10,'(a)') 'ASCII'
        write(10,'(a)') 'DATASET STRUCTURED_GRID'
        write(10,'(a,3(1x,i0))') 'DIMENSIONS', nx, ny, nz
        
        write(10,*) 'POINTS', nx*ny*nz, 'float'
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(10,'(3F12.6)') x(i), y(j), z(k)
                end do
            end do
        end do
        
        write(10,*) 'POINT_DATA', nx*ny*nz
        write(10,*) 'VECTORS velocity float'
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(10,'(3F12.6)') u(i,j,k), v(i,j,k), w(i,j,k)
                end do
            end do
        end do
        
        write(10,*) 'SCALARS pressure float'
        write(10,*) 'LOOKUP_TABLE default'
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(10,'(F12.6)') p(i,j,k)
                end do
            end do
        end do
        
        close(10)
    end subroutine output_results
    
    subroutine output_final_results()
        // Local variables
        size_t :: i, j, k, count_points
        real :: max_velocity, avg_velocity, max_pressure, min_pressure
        real :: total_flow_rate, vel_magnitude
        
        // Output final results and statistics
        open(unit=11, file='stable_piso_final_results.txt', status='replace')
        
        write(11,*) '3D Navier-Stokes PISO Flow Solver Results (Stable Version)'
        write(11,*) '========================================================='
        write(11,*) 'Grid dimensions:', nx, 'x', ny, 'x', nz
        write(11,*) 'Reynolds number:', Re
        write(11,*) 'Time step:', dt
        write(11,*) 'Tube radius:', R
        write(11,*) 'Tube length:', L
        write(11,*) 'PISO correctors:', n_correctors
        write(11,*) 'Under-relaxation factors: alpha=', alpha, 'beta=', beta
        
        // Calculate flow statistics
        max_velocity = 0.0
        avg_velocity = 0.0
        max_pressure = -1.0e10
        min_pressure = 1.0e10
        total_flow_rate = 0.0
        count_points = 0
        
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    if (radius(j,k) <= R) then
                        vel_magnitude = sqrt(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)
                        max_velocity = max(max_velocity, vel_magnitude)
                        avg_velocity = avg_velocity + vel_magnitude
                        max_pressure = max(max_pressure, p(i,j,k))
                        min_pressure = min(min_pressure, p(i,j,k))
                        count_points = count_points + 1
                        
                        if (i == nx/2) then  // Calculate flow rate at middle
                            total_flow_rate = total_flow_rate + u(i,j,k) * dy * dz
                        end if
                    end if
                end do
            end do
        end do
        
        if (count_points > 0) then
            avg_velocity = avg_velocity / real(count_points,8)
        end if
        
        write(11,*) ''
        write(11,*) 'Flow Statistics:'
        write(11,*) 'Maximum velocity:', max_velocity
        write(11,*) 'Average velocity:', avg_velocity
        write(11,*) 'Maximum pressure:', max_pressure
        write(11,*) 'Minimum pressure:', min_pressure
        write(11,*) 'Flow rate at middle:', total_flow_rate
        
        close(11)
        
        write(*,*) 'Final results written to stable_piso_final_results.txt'
    end subroutine output_final_results

end program navier_stokes_3d_piso_stable 
*/