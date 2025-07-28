import numpy as np

def calculate_n_next_timestep(phi_j, n_j_current, phi_flux_right, phi_flux_left, dx, dt):
    """
    Calculate n_j^(k+1) using the explicit mass transport equation (A1):
    
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
    
    Returns:
    flux: [phi*n*u_g] at the interface
    """
    flux = 0.5 * (phi_left * n_left * u_g_left + phi_right * n_right * u_g_right)
    return flux

# Example: Calculate n for a single grid point
if __name__ == "__main__":
    # Example parameters for a single grid point calculation
    phi_j = 0.3  # Porosity at point j
    n_j_current = 1.0  # Current concentration at point j
    
    # Neighboring values for flux calculations
    phi_left = 0.3
    n_left = 0.8
    u_g_left = 1.0
    
    phi_right = 0.3
    n_right = 1.2
    u_g_right = 1.0
    
    # Grid and time parameters
    dx = 0.01
    dt = 0.001
    
    # Calculate fluxes at interfaces
    phi_flux_left = calculate_interface_flux(phi_left, n_left, u_g_left, phi_j, n_j_current, 1.0)
    phi_flux_right = calculate_interface_flux(phi_j, n_j_current, 1.0, phi_right, n_right, u_g_right)
    
    # Calculate concentration at next time step
    n_j_next = calculate_n_next_timestep(phi_j, n_j_current, phi_flux_right, phi_flux_left, dx, dt)
    
    print("Single Grid Point Calculation:")
    print(f"Current concentration n_j^k = {n_j_current}")
    print(f"Left flux [phi*n*u_g]_(j-1/2) = {phi_flux_left:.6f}")
    print(f"Right flux [phi*n*u_g]_(j+1/2) = {phi_flux_right:.6f}")
    print(f"Flux difference = {phi_flux_right - phi_flux_left:.6f}")
    print(f"Next concentration n_j^(k+1) = {n_j_next:.6f}")
    print(f"Change in concentration = {n_j_next - n_j_current:.6f}")