import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import biceps as bp


# Define the expressions for the bottom and surface elevations
def zb_expr(x): return 0.0  # Flat bottom surface
def zs_expr(x): return H + z0*np.cos(np.pi*x/L)  # Undulating surface elevation


# Define external force functions
def force_x(x, z): return 0.0  # No horizontal force
def force_z(x, z): return -1e-3*bp.GRAVITY*bp.ICE_DENSITY  # Gravity-driven force in the z-direction


# Define domain and grid parameters
x0 = 0.0
x1 = 100.0
L = x1 - x0  # Length of the domain
H = 1  # Mean height of the domain
z0 = 0.5  # Amplitude of surface undulation

nx = 50  # Number of elements in x-direction
nz = 5  # Number of elements in z-direction
deg_u = 2  # Polynomial degree for velocity field
deg_p = 1  # Polynomial degree for pressure field
gp = 5  # Number of Gauss points per element
cell_type = bp.CELL_TYPE_2D.TRIANGLE_LEFT  # Type of triangular element

# Create structured meshes for velocity and pressure fields
u_mesh_2d = bp.StructuredMesh(nx, nz, deg_u, cell_type)
p_mesh_2d = bp.StructuredMesh(nx, nz, deg_p, cell_type)

# Extrude the meshes in the x and z directions
u_mesh_2d.extrude_x(x0, x1)
u_mesh_2d.extrude_z(zb_expr, zs_expr)
p_mesh_2d.extrude_x(x0, x1)
p_mesh_2d.extrude_z(zb_expr, zs_expr)

# Define material properties and solver parameters
A = 100.0  # Coefficient for viscosity term
n_i = 3.0  # Flow law exponent (nonlinear)
fssa_version = bp.FSSA_VERSION.FSSA_NONE  # No FSSA stabilization
fssa_param = 0  # No additional parameter for FSSA
eps_reg_2 = 1e-10  # Regularization parameter

max_iter = 100  # Maximum number of iterations for solver
stol = 1e-10  # Convergence tolerance for solver

# Apply boundary conditions
ux_boundary_id = (
    bp.DOMAIN_IDS_2D.NORTH_WEST_ID |
    bp.DOMAIN_IDS_2D.WEST_ID |
    bp.DOMAIN_IDS_2D.BED_ID |
    bp.DOMAIN_IDS_2D.EAST_ID |
    bp.DOMAIN_IDS_2D.NORTH_EAST_ID
)
uz_boundary_id = bp.DOMAIN_IDS_2D.BED_ID

# Assemble and solve the nonlinear pStokes problem using fixed-point iterations
psp = bp.pStokesProblem(u_mesh_2d, p_mesh_2d)
psp.solve_nonlinear_system_picard(
    A, n_i, eps_reg_2, fssa_version, fssa_param, force_x, force_z, ux_boundary_id, uz_boundary_id, max_iter, stol, gp
)

# Extract velocity field solutions
ux_vec = psp.velocity_x().vals
uz_vec = psp.velocity_z().vals

# Extract velocity at mesh vertices
vdofs = u_mesh_2d.extract_vertex_dof_inds(bp.DOMAIN_IDS_2D.DOMAIN_ID)
ux_vec_p1 = ux_vec[vdofs]

# Plot solution
tri = Triangulation(p_mesh_2d.pmat[:, 0], p_mesh_2d.pmat[:, 1], p_mesh_2d.cmat)
fig, ax = plt.subplots()
contour = ax.tricontourf(tri, ux_vec_p1, levels=10)  # Contour plot of velocity field
plt.colorbar(contour)

plt.show()
