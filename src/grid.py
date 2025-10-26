import numpy as np

class StaggeredGrid:
    """
    Generates and holds all grid information for a 2D staggered grid.
    """
    def __init__(self, length, width, nx, ny):
        self.length = length
        self.width = width
        self.nx = nx
        self.ny = ny
        
        # Main control volume (P-grid) dimensions
        self.ipmax = self.nx + 2
        self.jpmax = self.ny + 2
        
        # U-velocity grid dimensions
        self.iumax = self.nx + 1
        self.jumax = self.ny + 2
        
        # V-velocity grid dimensions
        self.ivmax = self.nx + 2
        self.jvmax = self.ny + 1
        
        # T-scalar grid (same as P-grid)
        self.itmax = self.nx + 2
        self.jtmax = self.ny + 2
        
        self.generate_grids()

    def generate_grids(self):
        """Generates the coordinate arrays for P, U, V, and T control volumes."""
        # --- P-Grid (and T-Grid) Generation ---
        self.dxp = self.length / self.nx
        self.dyp = self.width / self.ny
        
        xp_nodes = np.zeros(self.ipmax)
        xp_nodes[1] = self.dxp / 2
        for i in range(2, self.ipmax - 1):
            xp_nodes[i] = xp_nodes[i-1] + self.dxp
        xp_nodes[-1] = self.length
        
        yp_nodes = np.zeros(self.jpmax)
        yp_nodes[1] = self.dyp / 2
        for j in range(2, self.jpmax - 1):
            yp_nodes[j] = yp_nodes[j-1] + self.dyp
        yp_nodes[-1] = self.width
        
        # --- U-Grid Generation ---
        xu_nodes = np.linspace(0, self.length, self.iumax)
        yu_nodes = yp_nodes
        
        # --- V-Grid Generation ---
        xv_nodes = xp_nodes
        yv_nodes = np.linspace(0, self.width, self.jvmax)
        
        # --- Create Meshgrids ---
        self.xp, self.yp = np.meshgrid(xp_nodes, yp_nodes, indexing='ij')
        self.xt, self.yt = self.xp, self.yp # T-grid is co-located with P-grid
        self.xu, self.yu = np.meshgrid(xu_nodes, yu_nodes, indexing='ij')
        self.xv, self.yv = np.meshgrid(xv_nodes, yv_nodes, indexing='ij')
