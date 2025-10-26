from src.grid import StaggeredGrid
from src.fvm_solver import Solver
from src.case_setup import set_fluid_properties, apply_boundary_conditions_cavity
from src.postprocessing import plot_velocity_contours, plot_pressure_contour, plot_streamlines

def main():
    """
    Main script to set up and run the Lid-Driven Cavity simulation.
    """
    # 1. --- Simulation Setup ---
    # Geometry
    length = 1.0
    width = 1.0
    nx = 32  # Number of cells in x-direction
    ny = 32  # Number of cells in y-direction

    # Fluid properties based on Reynolds number
    # Re = rho * U_lid * L / mu
    U_lid = 1.0
    Re = 400
    rho = 1.0
    mu = (rho * U_lid * length) / Re
    
    fluid_props = set_fluid_properties(re=Re, rho=rho, mu=mu)
    
    # Control parameters for the solver
    controls = {
        'iter_max': 15000,
        'convergence_tol': 1e-6,
        'monitor_interval': 100,
        'alpha_u': 0.7,  # Relaxation factor for u
        'alpha_v': 0.7,  # Relaxation factor for v
        'alpha_p': 0.3,  # Relaxation factor for p
    }

    # 2. --- Create Grid ---
    print(f"Setting up a {nx}x{ny} grid for a lid-driven cavity.")
    grid = StaggeredGrid(length, width, nx, ny)

    # 3. --- Initialize and Run Solver ---
    print(f"Running simulation for Re = {Re}...")
    solver = Solver(grid, fluid_props, controls)
    
    # The `apply_boundary_conditions_cavity` function needs the lid velocity,
    # so we pass it through the solver object. A bit of a hack, but simple.
    # A more advanced setup might pass a dictionary of BC parameters.
    solver.U_lid = U_lid 
    
    solver.solve(apply_boundary_conditions_cavity)

    # 4. --- Post-processing ---
    if solver.is_converged:
        print("Simulation converged. Generating plots...")
        plot_velocity_contours(solver)
        plot_pressure_contour(solver)
        plot_streamlines(solver)
        print("Plots saved in the root directory.")
    else:
        print("Simulation did not converge. Still generating plots for the final state...")
        plot_velocity_contours(solver)
        plot_pressure_contour(solver)
        plot_streamlines(solver)
        print("Plots saved in the root directory.")


if __name__ == "__main__":
    main()
