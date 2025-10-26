# Python Incompressible Flow Solvers (Finite Volume Method)

![Language](https://img.shields.io/badge/Language-Python-blue.svg)
![Libraries](https://img.shields.io/badge/Libraries-NumPy%20%7C%20Matplotlib-orange.svg)
![Methods](https://img.shields.io/badge/Methods-FVM%20%7C%20SIMPLE%20%7C%20Artificial%20Compressibility-red.svg)

This repository contains custom, from-scratch Python solvers for the 2D incompressible Navier-Stokes equations. The project demonstrates the implementation of fundamental CFD principles and advanced numerical methods for solving classic fluid dynamics benchmark problems, without relying on external CFD libraries.

The primary goal is to showcase a deep, first-principles understanding of solver architecture, numerical schemes, and validation techniques.

## Core Numerical Methods Implemented

This collection of solvers demonstrates proficiency across a range of essential CFD techniques:

*   **Discretization:** Cell-centered **Finite Volume Method (FVM)** on structured grids.
*   **Grid Systems:**
    *   **Staggered Grid:** To prevent pressure-velocity decoupling and avoid spurious checkerboard pressure fields.
    *   **Stretched Curvilinear Grid:** For efficient mesh clustering in regions with high gradients (e.g., near walls and steps).
*   **Pressure-Velocity Coupling Algorithms:**
    *   **SIMPLE Algorithm:** The classic Semi-Implicit Method for Pressure-Linked Equations, an iterative segregation-based approach.
    *   **Artificial Compressibility Method:** A pseudo-time marching method that allows the use of explicit schemes (like Runge-Kutta) for incompressible flow.
*   **Schemes:**
    *   **Power-Law Scheme:** For robust and stable discretization of convection-diffusion terms.
    *   **2nd-Order Central Differencing:** For spatial discretization of fluxes.
*   **Stability & Time-Stepping:**
    *   **4th-Order Artificial Dissipation:** To suppress numerical oscillations and maintain stability in convection-dominated flows.
    *   **4th-Order Low-Storage Runge-Kutta (LSRK4):** For high-accuracy explicit pseudo-time integration with optimized memory usage.
*   **Linear System Solvers:**
    *   **Line-by-Line TDMA (Thomas Algorithm):** For solving the resulting algebraic equations.
    *   **Strongly Implicit Procedure (SIP):** A combination of Successive Under-Relaxation and TDMA.

---

## Showcased Simulations & Validations

### 1. Laminar Channel Flow (with Heat Transfer)

This case simulates the developing flow of viscous oil (Re=200) in a 2D channel, including the conjugate heat transfer problem. It serves as a validation for the **SIMPLE algorithm**.

![`U Contours for Channel Flow`](ch_u_ctr.png)

*Figure 3: U velocity contours at Re=200

![`V Countours for Channel Flow`](ch_v_ctr.png)

*Figure 4: U velocity contours at Re=200

**Results & Validation:**
The solver correctly captures the development of the parabolic velocity profile from a uniform inlet. Crucially, the calculated Fanning friction factor converges to the theoretical value of `24/Re = 0.12` for fully developed flow between parallel plates.

### 2. Lid-Driven Cavity Flow

This classic CFD benchmark simulates the vortex formation within a square cavity where the top lid moves at a constant velocity. It is an excellent test for the solver's robustness and accuracy in handling strong vortical flows.

![`U Contours for Lid-Driven Cavity`](cavity_u_ctr.png)

*Figure 3: U velocity contours at Re=200

![`V Countours for Lid-Driven Cavity`](cavity_v_ctr.png)

*Figure 4: U velocity contours at Re=200


## How to Run

### Prerequisites

*   Python 3.x
*   NumPy
*   Matplotlib

You can install the required packages using pip:
```sh
pip install numpy matplotlib
```

### Execution

To run a simulation, execute one of the main scripts from the root directory.

**To run the Channel Flow simulation:**
```sh
python main_channel.py
```

**To run the Backward-Facing Step simulation:**
```sh
python main_bfs.py
```

The scripts will print convergence residuals to the console. Upon completion, plot images (e.g., `u_contour.png`, `streamlines.png`) will be saved to the root directory.
