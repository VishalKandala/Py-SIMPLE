import numpy as np
import time

class Solver:
    """
    Finite Volume Method solver for 2D incompressible, steady-state flows
    using the SIMPLE algorithm on a staggered grid.
    """
    def __init__(self, grid, props, controls):
        self.grid = grid
        self.props = props
        self.controls = controls

        # --- Fluid Properties ---
        self.rho = self.props['rho']
        self.mu = self.props['mu']

        # --- Solution Arrays ---
        self.u = np.zeros((grid.iumax, grid.jumax))
        self.v = np.zeros((grid.ivmax, grid.jvmax))
        self.P = np.zeros((grid.ipmax, grid.jpmax))
        self.p_corr = np.zeros((grid.ipmax, grid.jpmax))

        # --- Momentum Equation Coefficients ---
        self.aUp = np.ones((grid.iumax, grid.jumax))
        self.aUe = np.zeros((grid.iumax, grid.jumax))
        self.aUw = np.zeros((grid.iumax, grid.jumax))
        self.aUn = np.zeros((grid.iumax, grid.jumax))
        self.aUs = np.zeros((grid.iumax, grid.jumax))
        self.bUp = np.zeros((grid.iumax, grid.jumax))

        self.aVp = np.ones((grid.ivmax, grid.jvmax))
        self.aVe = np.zeros((grid.ivmax, grid.jvmax))
        self.aVw = np.zeros((grid.ivmax, grid.jvmax))
        self.aVn = np.zeros((grid.ivmax, grid.jvmax))
        self.aVs = np.zeros((grid.ivmax, grid.jvmax))
        self.bVp = np.zeros((grid.ivmax, grid.jvmax))

        # --- Pressure Correction Equation Coefficients ---
        self.aPp = np.ones((grid.ipmax, grid.jpmax))
        self.aPe = np.zeros((grid.ipmax, grid.jpmax))
        self.aPw = np.zeros((grid.ipmax, grid.jpmax))
        self.aPn = np.zeros((grid.ipmax, jpmax))
        self.aPs = np.zeros((grid.ipmax, jpmax))
        self.bPp = np.zeros((grid.ipmax, jpmax))

        # --- Solver State ---
        self.iteration = 0
        self.is_converged = False
        self.residuals = {'u': 1e5, 'v': 1e5, 'mass': 1e5}
        self.mass_imbalance = 1.0

    def solve(self, boundary_conditions_func):
        """Main solution loop for the SIMPLE algorithm."""
        boundary_conditions_func(self)
        start_time = time.time()

        while self.iteration < self.controls['iter_max'] and not self.is_converged:
            self.iteration += 1

            self.calculate_u_coefficients()
            self.solve_u_momentum()

            self.calculate_v_coefficients()
            self.solve_v_momentum()
            
            self.calculate_p_correction_coefficients()
            self.solve_pressure_correction()
            
            self.correct_fields()
            self.update_boundary_conditions()
            
            self.check_convergence()
            self.print_residuals()

        end_time = time.time()
        if self.is_converged:
            print(f"Solution Converged in {self.iteration} iterations.")
        else:
            print(f"Solver stopped after reaching max iterations ({self.controls['iter_max']}).")
        print(f"Total time: {end_time - start_time:.2f} seconds.")

    def _hybrid_scheme_ahead(self, Pe):
        """Hybrid differencing scheme for positive flux direction."""
        return max(0, (1 - 0.1 * abs(Pe))**5) + max(-Pe, 0)

    def _hybrid_scheme_behind(self, Pe):
        """Hybrid differencing scheme for negative flux direction."""
        return self._hybrid_scheme_ahead(abs(Pe)) + max(Pe, 0)

    def calculate_u_coefficients(self):
        """Calculates coefficients for the U-momentum equation."""
        g = self.grid
        for iu in range(1, g.iumax - 1):
            for ju in range(1, g.jumax - 1):
                # Geometric params
                dxp = g.dxp
                dyp = g.dyp
                
                # Diffusion coefficients
                De = self.mu * dyp / dxp
                Dw = self.mu * dyp / dxp
                Dn = self.mu * dxp / dyp
                Ds = self.mu * dxp / dyp
                
                # Convective fluxes (mass flow rates)
                Fe = self.rho * 0.5 * (self.u[iu, ju] + self.u[iu + 1, ju]) * dyp
                Fw = self.rho * 0.5 * (self.u[iu - 1, ju] + self.u[iu, ju]) * dyp
                Fn = self.rho * 0.5 * (self.v[iu, ju] + self.v[iu + 1, ju]) * dxp
                Fs = self.rho * 0.5 * (self.v[iu, ju - 1] + self.v[iu + 1, ju - 1]) * dxp
                
                # Peclet numbers
                Pe_e, Pe_w = Fe / De, Fw / Dw
                Pe_n, Pe_s = Fn / Dn, Fs / Ds

                # Coefficients using Hybrid scheme
                self.aUe[iu, ju] = De * self._hybrid_scheme_ahead(Pe_e)
                self.aUw[iu, ju] = Dw * self._hybrid_scheme_behind(Pe_w)
                self.aUn[iu, ju] = Dn * self._hybrid_scheme_ahead(Pe_n)
                self.aUs[iu, ju] = Ds * self._hybrid_scheme_behind(Pe_s)
                
                self.aUp[iu, ju] = self.aUe[iu, ju] + self.aUw[iu, ju] + self.aUn[iu, ju] + self.aUs[iu, ju]
                self.bUp[iu, ju] = (self.P[iu, ju] - self.P[iu + 1, ju]) * dyp

    def solve_u_momentum(self):
        """Solves the U-momentum equation using line-by-line TDMA sweeps."""
        g = self.grid
        alpha_u = self.controls['alpha_u']
        
        # --- Implicit X-sweep (line by line in Y) ---
        for ju in range(1, g.jumax - 1):
            A = np.zeros((g.iumax, 3))
            B = np.zeros((g.iumax, 1))
            
            # Assemble system for this line
            for iu in range(1, g.iumax - 1):
                A[iu, 0] = -self.aUw[iu, ju]
                A[iu, 1] = self.aUp[iu, ju] / alpha_u
                A[iu, 2] = -self.aUe[iu, ju]
                B[iu, 0] = self.aUn[iu, ju] * self.u[iu, ju + 1] + self.aUs[iu, ju] * self.u[iu, ju - 1] + \
                           self.bUp[iu, ju] + (1 - alpha_u) / alpha_u * self.aUp[iu, ju] * self.u[iu, ju]
            
            # Apply BCs to the linear system
            B[1, 0] += self.aUw[1, ju] * self.u[0, ju] # West BC
            A[1, 0] = 0
            B[g.iumax-2, 0] += self.aUe[g.iumax-2, ju] * self.u[g.iumax-1, ju] # East BC
            A[g.iumax-2, 2] = 0
            
            # Solve with TDMA
            x = self.tdma_solver_x(A, B, g.iumax-2)
            self.u[1:g.iumax-1, ju] = x.flatten()
            
    def calculate_v_coefficients(self):
        """Calculates coefficients for the V-momentum equation."""
        g = self.grid
        for iv in range(1, g.ivmax - 1):
            for jv in range(1, g.jvmax - 1):
                dxp, dyp = g.dxp, g.dyp
                
                De = self.mu * dyp / dxp
                Dw = self.mu * dyp / dxp
                Dn = self.mu * dxp / dyp
                Ds = self.mu * dxp / dyp
                
                Fe = self.rho * 0.5 * (self.u[iv, jv] + self.u[iv, jv + 1]) * dyp
                Fw = self.rho * 0.5 * (self.u[iv - 1, jv] + self.u[iv - 1, jv + 1]) * dyp
                Fn = self.rho * 0.5 * (self.v[iv, jv] + self.v[iv, jv + 1]) * dxp
                Fs = self.rho * 0.5 * (self.v[iv, jv - 1] + self.v[iv, jv]) * dxp
                
                Pe_e, Pe_w = Fe / De, Fw / Dw
                Pe_n, Pe_s = Fn / Dn, Fs / Ds

                self.aVe[iv, jv] = De * self._hybrid_scheme_ahead(Pe_e)
                self.aVw[iv, jv] = Dw * self._hybrid_scheme_behind(Pe_w)
                self.aVn[iv, jv] = Dn * self._hybrid_scheme_ahead(Pe_n)
                self.aVs[iv, jv] = Ds * self._hybrid_scheme_behind(Pe_s)
                
                self.aVp[iv, jv] = self.aVe[iv, jv] + self.aVw[iv, jv] + self.aVn[iv, jv] + self.aVs[iv, jv]
                self.bVp[iv, jv] = (self.P[iv, jv] - self.P[iv, jv + 1]) * dxp
                
    def solve_v_momentum(self):
        """Solves the V-momentum equation using line-by-line TDMA sweeps."""
        g = self.grid
        alpha_v = self.controls['alpha_v']
        
        # --- Implicit Y-sweep (line by line in X) ---
        for iv in range(1, g.ivmax - 1):
            A = np.zeros((g.jvmax, 3))
            B = np.zeros((g.jvmax, 1))
            
            for jv in range(1, g.jvmax - 1):
                A[jv, 0] = -self.aVs[iv, jv]
                A[jv, 1] = self.aVp[iv, jv] / alpha_v
                A[jv, 2] = -self.aVn[iv, jv]
                B[jv, 0] = self.aVe[iv, jv] * self.v[iv + 1, jv] + self.aVw[iv, jv] * self.v[iv - 1, jv] + \
                           self.bVp[iv, jv] + (1 - alpha_v) / alpha_v * self.aVp[iv, jv] * self.v[iv, jv]

            B[1, 0] += self.aVs[iv, 1] * self.v[iv, 0] # South BC
            A[1, 0] = 0
            B[g.jvmax-2, 0] += self.aVn[iv, g.jvmax-2] * self.v[iv, g.jvmax-1] # North BC
            A[g.jvmax-2, 2] = 0
            
            x = self.tdma_solver_y(A, B, g.jvmax-2)
            self.v[iv, 1:g.jvmax-1] = x.flatten()
            
    def calculate_p_correction_coefficients(self):
        """Calculates coefficients for the pressure correction equation."""
        g = self.grid
        for ip in range(1, g.ipmax - 1):
            for jp in range(1, g.jpmax - 1):
                de = g.dyp / self.aUp[ip, jp] if self.aUp[ip, jp] != 0 else 0
                dw = g.dyp / self.aUp[ip - 1, jp] if self.aUp[ip-1, jp] != 0 else 0
                dn = g.dxp / self.aVp[ip, jp] if self.aVp[ip, jp] != 0 else 0
                ds = g.dxp / self.aVp[ip, jp - 1] if self.aVp[ip, jp-1] != 0 else 0
                
                self.aPe[ip, jp] = self.rho * de * g.dyp
                self.aPw[ip, jp] = self.rho * dw * g.dyp
                self.aPn[ip, jp] = self.rho * dn * g.dxp
                self.aPs[ip, jp] = self.rho * ds * g.dxp
                
                # Mass imbalance (source term for pressure correction)
                self.bPp[ip, jp] = self.rho * (self.u[ip-1, jp] - self.u[ip, jp]) * g.dyp + \
                                   self.rho * (self.v[ip, jp-1] - self.v[ip, jp]) * g.dxp
        
        # Apply BCs for pressure correction (zero flux)
        self.aPs[:, 1] = 0
        self.aPn[:, g.jpmax-2] = 0
        self.aPw[1, :] = 0
        self.aPe[g.ipmax-2, :] = 0
        
        # Central coefficient
        self.aPp = self.aPe + self.aPw + self.aPn + self.aPs
    
    def solve_pressure_correction(self):
        """Solves the pressure correction equation using TDMA sweeps."""
        g = self.grid
        self.p_corr = np.zeros_like(self.P) # Reset correction field
        
        for _ in range(10): # Iterate a few times to help the pressure field
            # --- Implicit X-sweep (line by line in Y) ---
            for jp in range(1, g.jpmax - 1):
                A = np.zeros((g.ipmax, 3))
                B = np.zeros((g.ipmax, 1))
                for ip in range(1, g.ipmax - 1):
                    A[ip, 0] = -self.aPw[ip, jp]
                    A[ip, 1] = self.aPp[ip, jp]
                    A[ip, 2] = -self.aPe[ip, jp]
                    B[ip, 0] = self.bPp[ip, jp] + self.aPn[ip, jp] * self.p_corr[ip, jp + 1] + self.aPs[ip, jp] * self.p_corr[ip, jp - 1]
                
                x = self.tdma_solver_x(A, B, g.ipmax-2)
                self.p_corr[1:g.ipmax-1, jp] = x.flatten()
    
    def correct_fields(self):
        """Corrects pressure and velocities based on the p_corr field."""
        g = self.grid
        alpha_p = self.controls['alpha_p']
        
        # Correct Pressure
        self.P += alpha_p * self.p_corr
        
        # Correct U-velocity
        for iu in range(1, g.iumax - 1):
            for ju in range(1, g.jumax - 1):
                d_u = g.dyp / self.aUp[iu, ju] if self.aUp[iu, ju] != 0 else 0
                self.u[iu, ju] += alpha_p * d_u * (self.p_corr[iu, ju] - self.p_corr[iu + 1, ju])

        # Correct V-velocity
        for iv in range(1, g.ivmax - 1):
            for jv in range(1, g.jvmax - 1):
                d_v = g.dxp / self.aVp[iv, jv] if self.aVp[iv, jv] != 0 else 0
                self.v[iv, jv] += alpha_p * d_v * (self.p_corr[iv, jv] - self.p_corr[iv, jv + 1])

    def update_boundary_conditions(self):
        """Updates Neumann boundary conditions after each iteration."""
        g = self.grid
        
        # Outlet U-velocity (mass conservation)
        mass_in = np.sum(self.rho * self.u[0, :] * g.dyp)
        mass_out_current = np.sum(self.rho * self.u[g.iumax - 2, :] * g.dyp)
        if mass_out_current != 0:
            self.u[g.iumax - 1, :] = self.u[g.iumax - 2, :] * (mass_in / mass_out_current)
        
        # Wall pressures (zero normal gradient)
        self.P[:, 0] = self.P[:, 1]
        self.P[:, g.jpmax-1] = self.P[:, g.jpmax-2]
        self.P[0, :] = self.P[1, :]
        # Outlet pressure is Dirichlet, so it's not updated here.

    def check_convergence(self):
        """Checks the residuals of the momentum and continuity equations."""
        g = self.grid
        tol = self.controls['convergence_tol']
        
        # --- U-Momentum Residual ---
        num, den = 0, 0
        for iu in range(1, g.iumax-1):
            for ju in range(1, g.jumax-1):
                num += abs(self.aUp[iu,ju]*self.u[iu,ju] - (self.aUe[iu,ju]*self.u[iu+1,ju] + self.aUw[iu,ju]*self.u[iu-1,ju] + self.aUn[iu,ju]*self.u[iu,ju+1] + self.aUs[iu,ju]*self.u[iu,ju-1]) - self.bUp[iu,ju])
                den += abs(self.aUp[iu,ju]*self.u[iu,ju])
        self.residuals['u'] = num / den if den != 0 else 0

        # --- V-Momentum Residual ---
        num, den = 0, 0
        for iv in range(1, g.ivmax-1):
            for jv in range(1, g.jvmax-1):
                num += abs(self.aVp[iv,jv]*self.v[iv,jv] - (self.aVe[iv,jv]*self.v[iv+1,jv] + self.aVw[iv,jv]*self.v[iv-1,jv] + self.aVn[iv,jv]*self.v[iv,jv+1] + self.aVs[iv,jv]*self.v[iv,jv-1]) - self.bVp[iv,jv])
                den += abs(self.aVp[iv,jv]*self.v[iv,jv])
        self.residuals['v'] = num / den if den != 0 else 0

        # --- Mass Residual ---
        self.residuals['mass'] = np.sum(np.abs(self.bPp))
        
        if self.residuals['u'] < tol and self.residuals['v'] < tol and self.residuals['mass'] < tol:
            self.is_converged = True

    def print_residuals(self):
        if self.iteration % self.controls['monitor_interval'] == 0 or self.is_converged:
            res_u = self.residuals['u']
            res_v = self.residuals['v']
            res_m = self.residuals['mass']
            print(f"Iter: {self.iteration:5d}, Res(u): {res_u:.4e}, Res(v): {res_v:.4e}, Res(mass): {res_m:.4e}")

    def tdma_solver_x(self, A, B, n):
        """Line-by-line TDMA solver for X-direction sweeps."""
        T = np.zeros((self.grid.iumax, 1))
        for i in range(1, n + 1):
            m = A[i, 0] / A[i - 1, 1]
            A[i, 1] -= m * A[i - 1, 2]
            B[i, 0] -= m * B[i - 1, 0]
        
        T[n] = B[n, 0] / A[n, 1]
        for i in range(n - 1, 0, -1):
            T[i] = (B[i, 0] - A[i, 2] * T[i + 1]) / A[i, 1]
        return T

    def tdma_solver_y(self, A, B, n):
        """Line-by-line TDMA solver for Y-direction sweeps."""
        T = np.zeros((self.grid.jvmax, 1))
        for i in range(1, n + 1):
            m = A[i, 0] / A[i - 1, 1]
            A[i, 1] -= m * A[i - 1, 2]
            B[i, 0] -= m * B[i - 1, 0]
            
        T[n] = B[n, 0] / A[n, 1]
        for i in range(n - 1, 0, -1):
            T[i] = (B[i, 0] - A[i, 2] * T[i + 1]) / A[i, 1]
        return T
