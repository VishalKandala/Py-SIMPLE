# Py-SIMPLE

A 2D Finite Volume Method (FVM) solver for incompressible, steady-state flows, written in Python with NumPy. This project implements the **SIMPLE (Semi-Implicit Method for Pressure Linked Equations)** algorithm on a staggered grid.

The solver is structured to be modular and can be easily configured to solve different fluid dynamics problems. Two example cases are provided:
1.  **Laminar Channel Flow:** Fully developed flow entering a rectangular channel.
2.  **Lid-Driven Cavity:** A classic CFD benchmark problem.

This project demonstrates the implementation of core CFD algorithms from first principles, without relying on external CFD libraries.

## Features

-   **Finite Volume Discretization:** The solver is based on the integral form of the governing equations.
-   **Staggered Grid:** Velocities are stored on cell faces, and pressure is stored at cell centers to prevent pressure-velocity decoupling.
-   **SIMPLE Algorithm:** An iterative algorithm for coupling the pressure and velocity fields.
-   **Hybrid Differencing Scheme:** Used for discretizing convective terms, providing stability.
-   **Line-by-Line TDMA Solver:** A Tri-Diagonal Matrix Algorithm is used to solve the discretized algebraic equations.
-   **Modular, Object-Oriented Design:** The code is structured into classes for the grid and solver, making it clean and extensible.
-   **Post-processing:** Includes scripts for generating contour plots of velocity and pressure, as well as streamlines.

## Numerical Method

The solver follows the standard SIMPLE algorithm for steady-state flows:
1.  Guess the pressure field, P*.
2.  Solve the momentum equations (for u* and v*) using the guessed pressure.
3.  Solve the pressure-correction equation, which is derived from the continuity equation, to find the pressure correction, P'.
4.  Correct the pressure field: P = P* + α_p * P'.
5.  Correct the face velocities to satisfy continuity: u = u* + u' and v = v* + v'.
6.  Repeat from step 2 until the solution residuals fall below a specified tolerance.

## Project Structure
├── src/ # Source code modules
│ ├── fvm_solver.py # The main Solver class
│ ├── grid.py # StaggeredGrid class
│ ├── case_setup.py # Boundary conditions and properties
│ └── postprocessing.py # Plotting functions
├── main_channel.py # Main script for the Channel Flow case
├── main_cavity.py # Main script for the Lid-Driven Cavity case
└── README.md

## How to Run

### Prerequisites

-   Python 3.x
-   NumPy
-   Matplotlib

You can install the required packages using pip:
```sh
pip install numpy matplotlib
```
Execution
To run a simulation, execute one of the main scripts from the root directory of the project.
To run the Channel Flow simulation:
code
Sh
python main_channel.py
To run the Lid-Driven Cavity simulation:
code
Sh
python main_cavity.py
The script will print the residuals at regular intervals. Once the solution converges, it will automatically generate and save plot images (u_velocity_contour.png, streamlines.png, etc.) in the root directory.
Sample Results (Channel Flow)
Below are sample output plots for a converged channel flow simulation at Re=200.
U-Velocity Contour
![alt text](u_velocity_contour.png)

The plot shows the development of the parabolic velocity profile along the channel.
Streamlines
![alt text](streamlines.png)

The streamlines are parallel, as expected for a fully-developed channel flow.
code
Code
**Note:** You will need to run the code to generate your own `u_velocity_contour.png` and `streamlines.png` files and place them in the root directory for the images to show up in the README.
