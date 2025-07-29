import numpy as np
import matplotlib.pyplot as plt

def solve_mass_transport_explicit(phi, n_initial, u_g, dx, dt, num_timesteps):
    """
    Solve the explicit mass transport equation (A1) from the image:
    
    (phi_j * n_j^(k+1) - phi_j * n_j^k) / dt + 
    (1/dx) * ([phi*n*u_g]_(j+1/2)^k - [phi*n*u_g]_(j-1/2)^k) = 0
    
    Rearranged to solve for n_j^(k+1):
    n_j^(k+1) = n_j^k - (dt/(phi_j * dx)) * ([phi*n*u_g]_(j+1/2)^k - [phi*n*u_g]_(j-1/2)^k)
    
    Parameters:
    phi: array of porosity values at grid points
    n_initial: initial concentration values at grid points
    u_g: gas velocity (can be array or scalar)
    dx: spatial grid spacing
    dt: time step
    num_timesteps: number of time steps to solve
    
    Returns:
    n_history: array containing concentration at each time step
    """
    
    # Initialize arrays
    num_points = len(n_initial)
    n_current = n_initial.copy()
    n_history = np.zeros((num_timesteps + 1, num_points))
    n_history[0] = n_initial
    
    # Time stepping loop
    for k in range(num_timesteps):
        n_new = np.zeros_like(n_current)
        
        # Calculate flux terms at cell interfaces (j+1/2 and j-1/2)
        for j in range(1, num_points - 1):  # Interior points only
            
            # Calculate [phi*n*u_g] at j+1/2 interface
            if isinstance(u_g, np.ndarray):
                # If u_g is an array, interpolate to interface
                phi_flux_right = 0.5 * (phi[j] * n_current[j] * u_g[j] + 
                                       phi[j+1] * n_current[j+1] * u_g[j+1])
            else:
                # If u_g is scalar, use it directly
                phi_flux_right = 0.5 * u_g * (phi[j] * n_current[j] + 
                                             phi[j+1] * n_current[j+1])
            
            # Calculate [phi*n*u_g] at j-1/2 interface
            if isinstance(u_g, np.ndarray):
                phi_flux_left = 0.5 * (phi[j-1] * n_current[j-1] * u_g[j-1] + 
                                      phi[j] * n_current[j] * u_g[j])
            else:
                phi_flux_left = 0.5 * u_g * (phi[j-1] * n_current[j-1] + 
                                            phi[j] * n_current[j])
            
            # Apply the explicit scheme to calculate n_j^(k+1)
            flux_difference = phi_flux_right - phi_flux_left
            n_new[j] = n_current[j] - (dt / (phi[j] * dx)) * flux_difference
        
        # Apply boundary conditions (assuming zero gradient)
        n_new[0] = n_new[1]
        n_new[-1] = n_new[-2]
        
        # Update for next time step
        n_current = n_new.copy()
        n_history[k + 1] = n_current
    
    return n_history

def check_stability_condition(u_g_max, dx, dt):
    """
    Check the CFL stability condition for the explicit scheme.
    For stability: dt <= dx / |u_g_max|
    """
    cfl_number = abs(u_g_max) * dt / dx
    print(f"CFL number: {cfl_number:.4f}")
    if cfl_number <= 1.0:
        print("Stability condition satisfied (CFL ≤ 1)")
    else:
        print("WARNING: Stability condition violated (CFL > 1)")
    return cfl_number

# Example usage
if __name__ == "__main__":
    # Define grid and parameters
    L = 1.0  # Domain length
    nx = 100  # Number of grid points
    dx = L / (nx - 1)
    x = np.linspace(0, L, nx)
    
    # Time parameters
    dt = 0.001  # Time step
    t_final = 0.5
    num_timesteps = int(t_final / dt)
    
    # Physical parameters
    phi = np.ones(nx) * 0.3  # Constant porosity
    u_g = 1.0  # Constant gas velocity
    
    # Initial condition (Gaussian pulse)
    n_initial = np.exp(-50 * (x - 0.3)**2)
    
    # Check stability
    check_stability_condition(u_g, dx, dt)
    
    # Solve the equation
    print(f"Solving for {num_timesteps} time steps...")
    n_solution = solve_mass_transport_explicit(phi, n_initial, u_g, dx, dt, num_timesteps)
    
    # Plot results
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(x, n_initial, 'b-', label='Initial (t=0)', linewidth=2)
    plt.plot(x, n_solution[num_timesteps//4], 'g--', label=f't={t_final/4:.3f}', linewidth=2)
    plt.plot(x, n_solution[num_timesteps//2], 'r--', label=f't={t_final/2:.3f}', linewidth=2)
    plt.plot(x, n_solution[-1], 'k-', label=f't={t_final:.3f}', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('n (concentration)')
    plt.title('Mass Transport Solution - Spatial Profiles')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    time_points = np.linspace(0, t_final, num_timesteps + 1)
    plt.plot(time_points, n_solution[:, nx//4], 'b-', label=f'x={x[nx//4]:.2f}', linewidth=2)
    plt.plot(time_points, n_solution[:, nx//2], 'r-', label=f'x={x[nx//2]:.2f}', linewidth=2)
    plt.plot(time_points, n_solution[:, 3*nx//4], 'g-', label=f'x={x[3*nx//4]:.2f}', linewidth=2)
    plt.xlabel('Time')
    plt.ylabel('n (concentration)')
    plt.title('Mass Transport Solution - Time Evolution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    print(f"Solution completed. Final time: {t_final}")
    print(f"Grid spacing dx: {dx:.6f}")
    print(f"Time step dt: {dt:.6f}")