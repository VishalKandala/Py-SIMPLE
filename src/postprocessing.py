import matplotlib.pyplot as plt
import numpy as np

def plot_velocity_contours(solver):
    """Plots and saves contours for u and v velocities."""
    grid = solver.grid
    
    plt.figure(figsize=(10, 5))
    plt.contourf(grid.xu, grid.yu, solver.u.T, levels=50, cmap='jet')
    plt.colorbar(label='U-velocity (m/s)')
    plt.title('U-Velocity Contour')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.axis('equal')
    plt.savefig('u_velocity_contour.png')
    plt.close()

    plt.figure(figsize=(10, 5))
    plt.contourf(grid.xv, grid.yv, solver.v.T, levels=50, cmap='jet')
    plt.colorbar(label='V-velocity (m/s)')
    plt.title('V-Velocity Contour')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.axis('equal')
    plt.savefig('v_velocity_contour.png')
    plt.close()

def plot_pressure_contour(solver):
    """Plots and saves contours for pressure."""
    grid = solver.grid
    
    plt.figure(figsize=(10, 5))
    plt.contourf(grid.xp, grid.yp, solver.P.T, levels=50, cmap='viridis')
    plt.colorbar(label='Pressure (Pa)')
    plt.title('Pressure Contour')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.axis('equal')
    plt.savefig('pressure_contour.png')
    plt.close()

def plot_streamlines(solver):
    """Plots streamlines of the flow field."""
    grid = solver.grid
    
    # Interpolate velocities to cell centers (P-grid) for streamline plotting
    u_centered = (solver.u[:-1, :] + solver.u[1:, :]) / 2
    v_centered = (solver.v[:, :-1] + solver.v[:, 1:]) / 2
    
    # Ensure dimensions match for plotting
    u_plot = u_centered[:, 1:-1]
    v_plot = v_centered[1:-1, :]
    
    plt.figure(figsize=(10, 5))
    plt.streamplot(grid.xp[1:-1, 1:-1], grid.yp[1:-1, 1:-1], u_plot.T, v_plot.T, density=2)
    plt.title('Streamlines')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.axis('equal')
    plt.savefig('streamlines.png')
    plt.close()
