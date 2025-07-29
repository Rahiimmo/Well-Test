def calculate_n_next_timestep(phi_j, n_j_current, phi_flux_right, phi_flux_left, dx, dt):
    """
    Calculate n_j^(k+1) using the explicit mass transport equation (A1):
    
    From the image equation (A1):
    (phi_j * n_j^(k+1) - phi_j * n_j^k) / dt + (1/dx) * ([phi*n*u_g]_(j+1/2)^k - [phi*n*u_g]_(j-1/2)^k) = 0
    
    Rearranging to solve for n_j^(k+1):
    n_j^(k+1) = n_j^k - (dt/(phi_j * dx)) * ([phi*n*u_g]_(j+1/2)^k - [phi*n*u_g]_(j-1/2)^k)
    
    Parameters:
    phi_j: porosity at grid point j
    n_j_current: current concentration at grid point j (n_j^k)
    phi_flux_right: [phi*n*u_g] at j+1/2 interface
    phi_flux_left: [phi*n*u_g] at j-1/2 interface
    dx: spatial grid spacing
    dt: time step
    
    Returns:
    n_j_next: concentration at next time step (n_j^(k+1))
    """
    flux_difference = phi_flux_right - phi_flux_left
    n_j_next = n_j_current - (dt / (phi_j * dx)) * flux_difference
    return n_j_next

def calculate_interface_flux(phi_left, n_left, u_g_left, phi_right, n_right, u_g_right):
    """
    Calculate [phi*n*u_g] at an interface using simple averaging.
    This corresponds to the flux terms in equation (A1).
    
    For the interface between points j and j+1, we average:
    [phi*n*u_g]_(j+1/2) = 0.5 * (phi_j * n_j * u_g_j + phi_(j+1) * n_(j+1) * u_g_(j+1))
    
    Returns:
    flux: [phi*n*u_g] at the interface
    """
    flux_left = phi_left * n_left * u_g_left
    flux_right = phi_right * n_right * u_g_right
    flux = 0.5 * (flux_left + flux_right)
    return flux

def solve_single_timestep(phi, n_current, u_g, dx, dt):
    """
    Solve one time step for all interior grid points.
    
    Parameters:
    phi: list of porosity values
    n_current: list of current concentration values
    u_g: gas velocity (scalar or list)
    dx: spatial grid spacing
    dt: time step
    
    Returns:
    n_new: list of concentrations at next time step
    """
    n_points = len(n_current)
    n_new = [0.0] * n_points
    
    # Handle boundary conditions (copy from neighbors)
    n_new[0] = n_current[0]  # Left boundary
    n_new[-1] = n_current[-1]  # Right boundary
    
    # Solve for interior points
    for j in range(1, n_points - 1):
        # Get velocity values
        if isinstance(u_g, list):
            u_g_left = u_g[j-1]
            u_g_center = u_g[j]
            u_g_right = u_g[j+1]
        else:
            u_g_left = u_g_center = u_g_right = u_g
        
        # Calculate fluxes at interfaces
        phi_flux_left = calculate_interface_flux(
            phi[j-1], n_current[j-1], u_g_left,
            phi[j], n_current[j], u_g_center
        )
        
        phi_flux_right = calculate_interface_flux(
            phi[j], n_current[j], u_g_center,
            phi[j+1], n_current[j+1], u_g_right
        )
        
        # Calculate new concentration
        n_new[j] = calculate_n_next_timestep(
            phi[j], n_current[j], phi_flux_right, phi_flux_left, dx, dt
        )
    
    return n_new

# Example usage
if __name__ == "__main__":
    print("Mass Transport Equation Solver")
    print("==============================")
    print("Solving equation (A1) from the image:")
    print("(φⱼnⱼᵏ⁺¹ - φⱼnⱼᵏ)/Δt + (1/Δx)([φnuₘ]ⱼ₊₁/₂ᵏ - [φnuₘ]ⱼ₋₁/₂ᵏ) = 0")
    print()
    
    # Example 1: Single grid point calculation
    print("Example 1: Single Grid Point Calculation")
    print("-" * 40)
    
    phi_j = 0.3
    n_j_current = 1.0
    
    # Neighboring values
    phi_left, n_left, u_g_left = 0.3, 0.8, 1.0
    phi_right, n_right, u_g_right = 0.3, 1.2, 1.0
    
    dx, dt = 0.01, 0.001
    
    # Calculate fluxes
    phi_flux_left = calculate_interface_flux(phi_left, n_left, u_g_left, phi_j, n_j_current, 1.0)
    phi_flux_right = calculate_interface_flux(phi_j, n_j_current, 1.0, phi_right, n_right, u_g_right)
    
    # Calculate next concentration
    n_j_next = calculate_n_next_timestep(phi_j, n_j_current, phi_flux_right, phi_flux_left, dx, dt)
    
    print(f"Current concentration nⱼᵏ = {n_j_current}")
    print(f"Left flux [φnuₘ]ⱼ₋₁/₂ = {phi_flux_left:.6f}")
    print(f"Right flux [φnuₘ]ⱼ₊₁/₂ = {phi_flux_right:.6f}")
    print(f"Flux difference = {phi_flux_right - phi_flux_left:.6f}")
    print(f"Next concentration nⱼᵏ⁺¹ = {n_j_next:.6f}")
    print(f"Change = {n_j_next - n_j_current:.6f}")
    print()
    
    # Example 2: Multiple grid points
    print("Example 2: Multiple Grid Points (5 points)")
    print("-" * 40)
    
    # Setup grid
    phi = [0.3, 0.3, 0.3, 0.3, 0.3]
    n_current = [0.5, 1.0, 1.5, 1.0, 0.5]
    u_g = 1.0
    dx = 0.02
    dt = 0.001
    
    print(f"Initial concentrations: {n_current}")
    
    # Solve one time step
    n_new = solve_single_timestep(phi, n_current, u_g, dx, dt)
    
    print(f"After one time step:   {[f'{x:.6f}' for x in n_new]}")
    
    # Calculate changes
    changes = [n_new[i] - n_current[i] for i in range(len(n_current))]
    print(f"Changes:               {[f'{x:.6f}' for x in changes]}")
    
    # Check stability condition
    cfl_number = abs(u_g) * dt / dx
    print(f"\nCFL number: {cfl_number:.4f} (should be ≤ 1.0 for stability)")
    
    if cfl_number <= 1.0:
        print("✓ Stability condition satisfied")
    else:
        print("⚠ WARNING: Stability condition violated!")