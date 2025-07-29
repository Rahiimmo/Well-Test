# Analysis of your mass transport implementation
# ============================================

def phi(j, dx, L_g, L_w, L_k):
    '''
    Assign porosity to the whole wellbore.
    Cemented region has a porosity of 0.01 and we assume porosity to be 1 in the uncemented region.
    :param j: index of the node
    Returns:
    porosity
    
    ISSUE FIXED: Original had (dx - 1/2) which should be (dx * (j - 1/2))
    '''
    x = j * dx - dx/2  # Position at grid point j (assuming grid starts at j=1)
    if L_g + L_w < x <= L_k:
        return 0.01
    else:
        return 1

def calculate_phi_at_interface(phi_left, phi_right):
    '''
    Calculate phi_(j+1/2)
    Parameters:
    phi_left: porosity at the grid point j (phi_j)
    phi_right: porosity at the grid point j+1 (phi_(j+1))

    Returns:
    phi at the interface (j+1/2)  phi_(j+1/2)
    
    NOTE: Using minimum is consistent with equation (A2) from the image
    '''
    return min(phi_left, phi_right)

def calculate_interface_flux(phi_left, phi_right, n_left, n_right, u_g):
    '''
    Calculate [phi*n*u_g] at the interface using equation (A2).
    This corresponds to the flux terms in equation (A1).

    From equation (A2) in the image:
    [phi*n*u_g]_(j+1/2) = {
        phi_(j+1/2) * n_j * u_g_(j+1/2),     if u_g_(j+1/2) >= 0
        phi_(j+1/2) * n_(j+1) * u_g_(j+1/2), if u_g_(j+1/2) < 0
    }

    Returns:
    flux: [phi*n*u_g] at the interface
    '''
    phi_interface = calculate_phi_at_interface(phi_left, phi_right)
    if u_g >= 0:
        return phi_interface * n_left * u_g   # Upwind: use left value
    else:
        return phi_interface * n_right * u_g  # Upwind: use right value

def calculate_n_next_timestep(n_current, n_right, n_left, phi_current, phi_right, phi_left, 
                             u_g_right, u_g_left, dt, dx):
    '''
    Calculate n_j_(k+1) based on equation (A1) explicitly.
    
    From equation (A1):
    (phi_j * n_j_(k+1) - phi_j * n_j_k) / dt + (1/dx) * ([phi*n*u_g]_(j+1/2)_k - [phi*n*u_g]_(j-1/2)_k) = 0

    Rearranging to solve for n_j_(k+1):
    n_j_(k+1) = n_j_k - (dt/(phi_j * dx)) * ([phi*n*u_g]_(j+1/2)_k - [phi*n*u_g]_(j-1/2)_k)
    
    Parameters:
    n_current: current concentration at grid point j (n_j_k)
    n_right: current concentration at grid point j+1 (n_(j+1)_k)
    n_left: current concentration at grid point j-1 (n_(j-1)_k)
    phi_current: porosity at grid point j (phi_j)
    phi_right: porosity at grid point j+1 (phi_(j+1))
    phi_left: porosity at grid point j-1 (phi_(j-1))
    u_g_right: gas velocity at j+1/2 interface (u_g_(j+1/2))
    u_g_left: gas velocity at j-1/2 interface (u_g_(j-1/2))
    dt: time step
    dx: spatial step
    
    Returns:
    n_j_next: concentration at next time step (n_j_(k+1))
    '''
    # Calculate flux at j-1/2 interface (left side)
    phi_flux_left = calculate_interface_flux(phi_left, phi_current, n_left, n_current, u_g_left)
    
    # Calculate flux at j+1/2 interface (right side)  
    phi_flux_right = calculate_interface_flux(phi_current, phi_right, n_current, n_right, u_g_right)
    
    # Apply explicit scheme
    flux_difference = phi_flux_right - phi_flux_left
    n_next = n_current - (dt / (phi_current * dx)) * flux_difference
    
    return n_next

# IMPROVED MAIN LOOP WITH BETTER STRUCTURE
# ========================================

def solve_transport_timestep(n_list, m_list, P_g_list, P_w_list, S_w_list, S_g_list, 
                           U_g_list, U_w_list, k, N, dt, dx, C_g, C_w, rho_wr,
                           L_g, L_w, L_k):
    '''
    Solve one time step of the transport equations for all spatial points.
    
    IMPROVEMENTS:
    1. Separated the transport calculation from saturation updates
    2. Added proper error handling
    3. Better organization of the calculation steps
    '''
    
    # Store new values temporarily
    n_new = n_list[k + 1, :].copy()  # Initialize with current values
    m_new = m_list[k + 1, :].copy()
    
    # Calculate transport for interior points
    for j in range(1, N-1):  # Avoid boundary points
        
        # Current state variables
        rho_g = P_g_list[k, j] / C_g
        rho_w = P_w_list[k, j] / C_w + rho_wr
        
        # Update mass densities
        m_list[k, j] = S_w_list[k, j] * rho_w
        n_list[k, j] = S_g_list[k, j] * rho_g
        
        # Get porosity values
        phi_j = phi(j, dx, L_g, L_w, L_k)
        phi_j_plus1 = phi(j + 1, dx, L_g, L_w, L_k)
        phi_j_minus1 = phi(j - 1, dx, L_g, L_w, L_k)
        
        # Calculate new gas mass density (n)
        n_new[j] = calculate_n_next_timestep(
            n_list[k, j], n_list[k, j + 1], n_list[k, j - 1],
            phi_j, phi_j_plus1, phi_j_minus1,
            U_g_list[k, j], U_g_list[k, j - 1],  # Note: j-1 for left interface
            dt, dx
        )
        
        # Calculate new water mass density (m)
        m_new[j] = calculate_n_next_timestep(
            m_list[k, j], m_list[k, j + 1], m_list[k, j - 1],
            phi_j, phi_j_plus1, phi_j_minus1,
            U_w_list[k, j], U_w_list[k, j - 1],  # Note: j-1 for left interface
            dt, dx
        )
        
        # Update arrays
        n_list[k + 1, j] = n_new[j]
        m_list[k + 1, j] = m_new[j]
        
        # Calculate intermediate saturations
        s_g_star = n_new[j] / rho_g
        s_w_star = m_new[j] / rho_w
        
        # Normalize saturations
        total_sat = s_w_star + s_g_star
        if total_sat > 0:  # Avoid division by zero
            S_w_list[k + 1, j] = s_w_star / total_sat
            S_g_list[k + 1, j] = s_g_star / total_sat
        else:
            # Handle edge case
            S_w_list[k + 1, j] = S_w_list[k, j]
            S_g_list[k + 1, j] = S_g_list[k, j]
    
    # Check stability condition
    max_velocity = 0
    for j in range(1, N-1):
        max_velocity = max(max_velocity, abs(U_g_list[k, j]), abs(U_w_list[k, j]))
    
    cfl_number = max_velocity * dt / dx
    if cfl_number > 1.0:
        print(f"⚠ WARNING: Stability condition violated at timestep {k}!")
        print(f"CFL number: {cfl_number:.4f} (should be ≤ 1.0)")
        print(f"Consider reducing dt or increasing dx")
        return False  # Indicate instability
    
    return True  # Stable timestep

# EXAMPLE USAGE OF IMPROVED VERSION
# =================================

def run_simulation_example():
    '''
    Example of how to use the improved transport solver
    '''
    
    # Example parameters (you'll need to adjust these to your actual values)
    M = 100  # Number of time steps
    N = 50   # Number of spatial points
    dt = 0.001
    dx = 0.02
    C_g = 1.0
    C_w = 1.0
    rho_wr = 1000.0
    L_g = 0.1
    L_w = 0.2
    L_k = 0.8
    
    # Initialize arrays (you'll have your own initialization)
    import numpy as np
    n_list = np.zeros((M, N))
    m_list = np.zeros((M, N))
    P_g_list = np.ones((M, N)) * 101325  # Example pressure
    P_w_list = np.ones((M, N)) * 101325
    S_w_list = np.ones((M, N)) * 0.5     # Example saturations
    S_g_list = np.ones((M, N)) * 0.5
    U_g_list = np.ones((M, N)) * 0.1     # Example velocities
    U_w_list = np.ones((M, N)) * 0.1
    
    # Main simulation loop
    for k in range(1, M-1):
        stable = solve_transport_timestep(
            n_list, m_list, P_g_list, P_w_list, S_w_list, S_g_list,
            U_g_list, U_w_list, k, N, dt, dx, C_g, C_w, rho_wr,
            L_g, L_w, L_k
        )
        
        if not stable:
            print(f"Simulation stopped due to instability at timestep {k}")
            break
        
        if k % 10 == 0:  # Print progress every 10 steps
            print(f"Completed timestep {k}/{M-1}")

if __name__ == "__main__":
    print("Mass Transport Code Analysis and Improvements")
    print("=" * 50)
    print()
    print("Key Issues Found in Original Code:")
    print("1. phi(j) function had incorrect position calculation")
    print("2. Missing dt, dx parameters in calculate_n_next_timestep")
    print("3. No boundary condition handling")
    print("4. Potential division by zero in saturation normalization")
    print("5. CFL check could be more informative")
    print()
    print("Improvements Made:")
    print("1. Fixed position calculation in phi(j)")
    print("2. Added proper parameter passing")
    print("3. Better error handling and stability checking")
    print("4. Cleaner separation of concerns")
    print("5. More robust saturation normalization")
    print()
    print("Your implementation of the upwind scheme (equation A2) is correct!")
    print("The use of min() for interface porosity is also appropriate.")