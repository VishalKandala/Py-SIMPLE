import numpy as np

def set_fluid_properties(re=200, pr=0.71, rho=1.225, mu=1.81e-5):
    """Returns a dictionary of fluid properties."""
    return {
        're': re,
        'pr': pr,
        'rho': rho,
        'mu': mu,
        'k': 2.57e-2,  # Thermal conductivity (example for air)
        'cp': 1005,    # Specific heat (example for air)
    }

def apply_boundary_conditions_channel(solver):
    """Applies boundary conditions for a fully-developed channel flow case."""
    grid = solver.grid
    
    # Calculate inlet velocity from Reynolds number
    dh = 2 * grid.width
    u_inlet = (solver.props['re'] * solver.props['mu']) / (solver.props['rho'] * dh)
    
    # --- Velocity BCs ---
    # U-velocity
    solver.u[0, 1:grid.jumax-1] = u_inlet  # Inlet
    solver.u[:, 0] = 0.0                   # Bottom wall (no-slip)
    solver.u[:, grid.jumax-1] = 0.0        # Top wall (no-slip)
    
    # V-velocity
    solver.v[0, 1:grid.jvmax-1] = 0.0      # Inlet
    solver.v[:, 0] = 0.0                   # Bottom wall (no-slip)
    solver.v[:, grid.jvmax-1] = 0.0        # Top wall (no-slip)
    
    # --- Pressure BCs ---
    solver.P[grid.ipmax-1, 1:grid.jpmax-1] = 0.0 # Outlet pressure
    
    # --- Temperature BCs (if applicable) ---
    if hasattr(solver, 'T'):
        T_inlet = 300 # K
        T_wall = 350 # K
        solver.T[0, 1:grid.jtmax-1] = T_inlet
        solver.T[:, 0] = T_wall
        solver.T[:, grid.jtmax-1] = T_wall

def apply_boundary_conditions_cavity(solver):
    """Applies boundary conditions for a lid-driven cavity case."""
    grid = solver.grid
    # Get the lid velocity from the solver object, where it was placed by main_cavity.py
    u_lid = solver.U_lid 
    
    # --- Velocity BCs (all walls are no-slip, except top lid) ---
    # U-velocity
    solver.u[:, 0] = 0.0                  # Bottom wall
    solver.u[:, grid.jumax-1] = u_lid     # Top wall (lid)
    solver.u[0, :] = 0.0                  # Left wall
    solver.u[grid.iumax-1, :] = 0.0       # Right wall

    # V-velocity
    solver.v[:, 0] = 0.0                  # Bottom wall
    solver.v[:, grid.jvmax-1] = 0.0       # Top wall
    solver.v[0, :] = 0.0                  # Left wall
    solver.v[grid.ivmax-1, :] = 0.0       # Right wall

    # --- Pressure BC is setting one reference point ---
    solver.P[1, 1] = 0.0
