using PyCall
using Plots
import Printf.@printf
bp = pyimport("biceps")

GRAVITY = 9.8  # Gravitational acceleration
ICE_DENSITY = 910  # Ice mass density

x0 = 0.0  # Left end
x1 = 100.0  # Right end
L = x1 - x0  # Length of the domain
H = 1  # Mean height of the domain
z0 = 0.1  # Amplitude of surface undulation

py"""
import numpy as np

# Define the expressions for the bottom and surface elevations
def zb_expr(x): return 0.0  # Flat bottom surface
def zs_expr(x): return $H + $z0*np.cos(np.pi*x/$L)  # Undulating surface elevation


# Define external force functions
def force_x(x, z): return 0.0  # No horizontal force
def force_z(x, z): return -1e-3*$GRAVITY*$ICE_DENSITY  # Gravity force in the z-direction
"""

A = 100.0  # Ice softness parameter
n_i = 3.0  # glen exponent
eps_reg_2 = 1e-10  # Regularization parameter
fssa_version = bp.FSSA_VERSION.FSSA_NONE  # No FSSA stabilization
fssa_param = 0  # Stabilization parameter in FSSA

nx = 50  # Number of elements in x-direction
nz = 5  # Number of elements in z-direction
nt = 100  # Number of time steps
dt = 35  # Time step size

deg_u = 2  # Polynomial degree for velocity field
deg_p = 1  # Polynomial degree for pressure field
deg_h = 1  # Polynomial degree for height field

max_iter = 100  # Maximum number of iterations for solver
stol = 1e-6  # Convergence tolerance for solver
cell_type = bp.CELL_TYPE_2D.TRIANGLE_LEFT  # 2D mesh cell type

# Create structured meshes for velocity, pressure, and height fields
u_mesh_2d = bp.StructuredMesh(nx, nz, deg_u, cell_type)
p_mesh_2d = bp.StructuredMesh(nx, nz, deg_p, cell_type)

# Extrude the meshes in the x and z directions
u_mesh_2d.extrude_x(x0, x1)
u_mesh_2d.extrude_z(py"zb_expr", py"zs_expr")
p_mesh_2d.extrude_x(x0, x1)
p_mesh_2d.extrude_z(py"zb_expr", py"zs_expr")

# Extract degrees of freedom for surface nodes
sdofs_u = u_mesh_2d.extract_dof_inds(bp.DOMAIN_IDS_2D.SURFACE_ID)
sdofs_h = u_mesh_2d.extract_vertex_dof_inds(bp.DOMAIN_IDS_2D.SURFACE_ID)
# Compensate for julia being 1-indexed
sdofs_u = 1 .+ sdofs_u
sdofs_h = 1 .+ sdofs_h

# Extract surface coordinates
spmat_u = copy(u_mesh_2d.pmat[sdofs_u, :])
spmat_h = copy(u_mesh_2d.pmat[sdofs_h, :])
xs_vec = copy(spmat_h[:, 1])
zs_vec = copy(spmat_h[:, 2])
# Project mesh to z=0
spmat_u[:, 2] .= 0
spmat_h[:, 2] .= 0

# Create 1D meshes for velocity and height on the surface
u_mesh_1d = bp.IntervalMesh(spmat_u, deg_u)
h_mesh_1d = bp.IntervalMesh(spmat_h, deg_h)

# Initialize FEM functions for height, velocity, and accumulation
ux_func = bp.FEMFunction1D(u_mesh_1d)
uz_func = bp.FEMFunction1D(u_mesh_1d)
h0_func = bp.FEMFunction1D(h_mesh_1d)
ac_func = bp.FEMFunction1D(h_mesh_1d)

# For L2 norm calculation
h0_func.assemble_mass_matrix()

# Initialize the pStokes
psp = bp.pStokesProblem(
    A, n_i, eps_reg_2, py"force_x", py"force_z", u_mesh_2d, p_mesh_2d
)

# Set Dirichlet BC masks
# Horizontal velocity component
psp.ux_dirichlet_bc_mask = (
    bp.DOMAIN_IDS_2D.NORTH_WEST_ID |
    bp.DOMAIN_IDS_2D.WEST_ID |
    bp.DOMAIN_IDS_2D.BED_ID |
    bp.DOMAIN_IDS_2D.EAST_ID |
    bp.DOMAIN_IDS_2D.NORTH_EAST_ID
)
# Vertical velocity component
psp.uz_dirichlet_bc_mask = bp.DOMAIN_IDS_2D.BED_ID

# Initialize the free-surface problem
fsp = bp.FreeSurfaceProblem(h_mesh_1d, u_mesh_1d)

# Plot initial surface profile
p = plot()
plot!(xs_vec, zs_vec, lw=2)
# Time-stepping loop
for k=1:nt
    # Assemble and solve the nonlinear pStokes problem
    psp.solve_nonlinear_system()
    # Clear lhs matrix and rhs vector
    psp.reset_system()

    # Extract velocity field solutions
    ux_vec = psp.velocity_x().vals
    uz_vec = psp.velocity_z().vals
    # Set free surface velocity
    ux_func.vals = ux_vec[sdofs_u]
    uz_func.vals = uz_vec[sdofs_u]
    # Set initial height
    h0_func.vals = copy(zs_vec)

    # Print surface energy and domain area
    @printf("||E|| = %.16f\n", h0_func.L2_norm())
    @printf("A = %.16f\n", u_mesh_2d.area())

    # Solve the free surface problem using explicit time stepping
    fsp.assemble_lhs_explicit()
    fsp.commit_lhs()
    fsp.assemble_rhs_explicit(
        h0_func, ux_func, uz_func, ac_func, dt
    )
    fsp.solve_linear_system()
    # Clear lhs matrix and rhs vector
    fsp.reset_system()

    # Update surface elevation
    global zs_vec = copy(fsp.zs_vec)

    # Update mesh with new surface elevation
    u_mesh_2d.extrude_z(zs_vec)
    p_mesh_2d.extrude_z(zs_vec)
end

# Plot final surface
plot!(xs_vec, zs_vec, lw=2)
gui()
# Pause program to display plot
readline()
