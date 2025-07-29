# CORRECTED VERSION OF YOUR CODE
# ===============================

def phi(j, dx, L_g, L_w, L_k):
    '''
    Assign porosity to the whole wellbore.
    Cemented region has a porosity of 0.01 and we assume porosity to be 1 in the uncemented region.
    :param j: index of the node
    Returns:
    porosity
    
    FIXED: Changed (dx - 1/2) to proper position calculation
    '''
    x = j * dx - dx/2  # Position at grid point j (assuming j starts from 1)
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
    '''
    return min(phi_left, phi_right)

def calculate_interface_flux(phi_left, phi_right, n_left, n_right, u_g):
    '''
    Calculate [phi*n*u_g] at the next interface (j+1/2) using equation (A2).
    This corresponds to the flux terms in equation (A1).

    Returns:
    flux: [phi*n*u_g] at the next interface. (j+1/2)
    We will use this function for water phase as well. (using m and u_w instead of n and u_g respectively)
    '''
    phi_average = calculate_phi_at_interface(phi_left, phi_right)
    if u_g >= 0:
        return phi_average * n_left * u_g
    else:
        return phi_average * n_right * u_g

def calculate_n_next_timestep(n_current, n_right, n_left, phi_current, phi_right, phi_left, 
                             u_g_right, u_g_left, dt, dx):
    '''
    We will calculate n_j_(k+1) based on (A1) Equation explicitly.
    From the equation (A1):
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
    dt: time step (ADDED)
    dx: spatial step (ADDED)

    Returns:
    n_j_next: concentration at next time step (n_j_(k+1))
    We will use this function for water phase as well. (using m and u_w instead of n and u_g respectively)
    '''
    phi_flux_left = calculate_interface_flux(phi_left, phi_current, n_left, n_current, u_g_left)
    phi_flux_right = calculate_interface_flux(phi_current, phi_right, n_current, n_right, u_g_right)
    return n_current - (dt / (phi_current * dx)) * (phi_flux_right - phi_flux_left)

# CORRECTED MAIN LOOP
# ===================

def corrected_main_loop():
    '''
    Your corrected main loop with improvements
    '''
    
    # You'll need to define these parameters based on your problem
    # M = number of time steps
    # N = number of spatial points  
    # dt, dx = time and space steps
    # L_g, L_w, L_k = geometry parameters
    # All the arrays: n_list, m_list, P_g_list, etc.
    
    # Example parameters (replace with your actual values)
    M = 100
    N = 50
    dt = 0.001
    dx = 0.02
    L_g = 0.1
    L_w = 0.2
    L_k = 0.8
    C_g = 1.0
    C_w = 1.0
    rho_wr = 1000.0
    
    # Initialize your arrays here (example)
    import numpy as np
    n_list = np.zeros((M, N))
    m_list = np.zeros((M, N))
    P_g_list = np.ones((M, N)) * 101325
    P_w_list = np.ones((M, N)) * 101325
    S_w_list = np.ones((M, N)) * 0.5
    S_g_list = np.ones((M, N)) * 0.5
    U_g_list = np.ones((M, N)) * 0.1
    U_w_list = np.ones((M, N)) * 0.1
    
    print("Starting corrected simulation...")
    
    for k in range(1, M-1):  # Changed range to avoid boundary issues
        
        # Process interior points only
        for j in range(1, N-1):  # Changed to avoid boundary points
            
            # Calculate densities
            rho_g = P_g_list[k, j] / C_g
            rho_w = P_w_list[k, j] / C_w + rho_wr
            
            # Update current mass densities
            m_list[k, j] = S_w_list[k, j] * rho_w
            n_list[k, j] = S_g_list[k, j] * rho_g
            
            # Calculate new gas mass density (n) - CORRECTED FUNCTION CALL
            n_list[k + 1, j] = calculate_n_next_timestep(
                n_list[k, j], n_list[k, j + 1], n_list[k, j - 1], 
                phi(j, dx, L_g, L_w, L_k), phi(j + 1, dx, L_g, L_w, L_k), phi(j - 1, dx, L_g, L_w, L_k), 
                U_g_list[k, j], U_g_list[k, j - 1],
                dt, dx  # ADDED MISSING PARAMETERS
            )
            
            # Calculate new water mass density (m) - CORRECTED FUNCTION CALL
            m_list[k + 1, j] = calculate_n_next_timestep(
                m_list[k, j], m_list[k, j + 1], m_list[k, j - 1], 
                phi(j, dx, L_g, L_w, L_k), phi(j + 1, dx, L_g, L_w, L_k), phi(j - 1, dx, L_g, L_w, L_k), 
                U_w_list[k, j], U_w_list[k, j - 1],
                dt, dx  # ADDED MISSING PARAMETERS
            )
            
            # Calculate intermediate saturations
            s_g_j_interface_star = n_list[k + 1, j] / rho_g
            s_w_j_interface_star = m_list[k + 1, j] / rho_w

            # Normalize the computed saturations - IMPROVED WITH SAFETY CHECK
            total_saturation = s_w_j_interface_star + s_g_j_interface_star
            if total_saturation > 1e-12:  # Avoid division by very small numbers
                s_w_j_interface = s_w_j_interface_star / total_saturation
                s_g_j_interface = s_g_j_interface_star / total_saturation
                
                # Update saturation arrays
                S_w_list[k + 1, j] = s_w_j_interface
                S_g_list[k + 1, j] = s_g_j_interface
            else:
                # Handle edge case - maintain previous values
                S_w_list[k + 1, j] = S_w_list[k, j]
                S_g_list[k + 1, j] = S_g_list[k, j]
                print(f"Warning: Very small total saturation at k={k}, j={j}")

        # Check stability condition - IMPROVED
        max_velocity = 0
        for j in range(1, N-1):
            max_velocity = max(max_velocity, abs(U_g_list[k, j]), abs(U_w_list[k, j]))
        
        cfl_number = max_velocity * dt / dx
        if cfl_number > 1.0:
            print(f"⚠ WARNING: Stability condition violated at timestep {k}!")
            print(f"CFL number: {cfl_number:.4f} (should be ≤ 1.0)")
            print(f"Max velocity: {max_velocity:.4f}")
            print(f"Consider reducing dt from {dt} or increasing dx from {dx}")
            # Optionally break or continue with warning
            # break  # Uncomment to stop simulation on instability
        
        # Progress reporting
        if k % 10 == 0:
            print(f"Completed timestep {k}/{M-1}, CFL = {cfl_number:.4f}")
    
    print("Simulation completed!")
    return n_list, m_list, S_w_list, S_g_list

# SUMMARY OF CORRECTIONS MADE
# ===========================

def print_corrections():
    print("CORRECTIONS MADE TO YOUR ORIGINAL CODE:")
    print("=" * 50)
    print()
    print("1. FIXED phi(j) function:")
    print("   Original: x = j * (dx - 1/2)  # INCORRECT")
    print("   Fixed:    x = j * dx - dx/2   # CORRECT")
    print()
    print("2. ADDED missing parameters to calculate_n_next_timestep:")
    print("   Added: dt, dx parameters")
    print("   Updated function calls to include these parameters")
    print()
    print("3. IMPROVED boundary handling:")
    print("   Changed loop range from range(1, N) to range(1, N-1)")
    print("   This avoids accessing n_list[k, j+1] when j = N-1")
    print()
    print("4. ENHANCED saturation normalization:")
    print("   Added check for total_saturation > 1e-12 to avoid division by zero")
    print("   Added fallback to previous values if normalization fails")
    print()
    print("5. IMPROVED CFL stability check:")
    print("   More informative error messages")
    print("   Shows actual values and suggestions")
    print()
    print("6. ADDED progress reporting:")
    print("   Shows simulation progress every 10 timesteps")
    print()
    print("Your core implementation was correct! These are mainly safety improvements.")

if __name__ == "__main__":
    print_corrections()
    print()
    print("Running example simulation...")
    # Uncomment the next line to run the example
    # corrected_main_loop()