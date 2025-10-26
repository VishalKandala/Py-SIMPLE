from src.grid import StaggeredGrid
from src.fvm_solver import Solver
from src.case_setup import set_fluid_properties, apply_boundary_conditions_channel
from src.postprocessing import plot_velocity_contours, plot_pressure_contour, plot_streamlines

def main():
    """
    Main script to set up and run the Channel Flow simulation.
    """
    # 1. --- Simulation Setup ---
    length = 2.0
    width = 0.02
    nx = 50
    ny = 25
    
    # Control parameters
    controls = {
        'iter_max': 25000,
        'convergence_tol': 1e-5,
        'monitor_interval': 100,
        'alpha_u': 0.7, # Relaxation factor for u
        'alpha_v': 0.7, # Relaxation factor for v
        'alpha_p': 0.3, # Relaxation factor for p
    }

    # 2. --- Create Grid and Properties ---
    grid = StaggeredGrid(length, width, nx, ny)
    fluid_props = set_fluid_properties(re=200)

    # 3. --- Initialize and Run Solver ---
    solver = Solver(grid, fluid_props, controls)
    solver.solve(apply_boundary_conditions_channel)

    # 4. --- Post-processing ---
    if solver.is_converged:
        print("Generating plots...")
        plot_velocity_contours(solver)
        plot_pressure_contour(solver)
        plot_streamlines(solver)
        print("Plots saved in the root directory.")

if __name__ == "__main__":
    main()
