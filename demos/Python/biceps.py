#
# Copyright (C) 2025 André Löfgren
#
# This file is part of Biceps.
#
# Biceps is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Biceps is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Biceps. If not, see <https://www.gnu.org/licenses/>.
#
import numpy as np
import matplotlib.pyplot as plt
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

A = 100.0  # Coefficient for viscosity term
n_i = 3.0  # Flow law exponent
fssa_version = bp.FSSA_VERSION.FSSA_NONE  # No FSSA stabilization
fssa_param = 0  # No additional parameter for FSSA
eps_reg_2 = 1e-10  # Regularization parameter

nx = 50  # Number of elements in x-direction
nz = 5  # Number of elements in z-direction
nt = 10  # Number of time steps
dt = 1  # Time step size

deg_u = 2  # Polynomial degree for velocity field
deg_p = 1  # Polynomial degree for pressure field
deg_h = 1  # Polynomial degree for height field
gp = 5  # Number of Gauss points per element

max_iter = 100  # Maximum number of iterations for solver
stol = 1e-10  # Convergence tolerance for solver

cell_type = bp.CELL_TYPE_2D.TRIANGLE_LEFT  # Type of triangular element

# Create structured meshes for velocity, pressure, and height fields
u_mesh_2d = bp.StructuredMesh(nx, nz, deg_u, cell_type)
p_mesh_2d = bp.StructuredMesh(nx, nz, deg_p, cell_type)

# Extrude the meshes in the x and z directions
u_mesh_2d.extrude_x(x0, x1)
u_mesh_2d.extrude_z(zb_expr, zs_expr)
p_mesh_2d.extrude_x(x0, x1)
p_mesh_2d.extrude_z(zb_expr, zs_expr)

# Extract degrees of freedom for surface nodes
sdofs_u = u_mesh_2d.extract_dof_inds(bp.DOMAIN_IDS_2D.SURFACE_ID)
sdofs_h = u_mesh_2d.extract_vertex_dof_inds(bp.DOMAIN_IDS_2D.SURFACE_ID)

# Extract surface coordinates
spmat_u = u_mesh_2d.pmat[sdofs_u, :].copy()
spmat_h = u_mesh_2d.pmat[sdofs_h, :].copy()
xs_vec = spmat_h[:, 0].copy()
zs_vec = spmat_h[:, 1].copy()
spmat_u[:, 1] = 0
spmat_h[:, 1] = 0

# Create 1D meshes for velocity and height on the surface
u_mesh_1d = bp.IntervalMesh(spmat_u, deg_u)
h_mesh_1d = bp.IntervalMesh(spmat_h, deg_h)

# Apply boundary conditions
ux_boundary_id = (
    bp.DOMAIN_IDS_2D.NORTH_WEST_ID |
    bp.DOMAIN_IDS_2D.WEST_ID |
    bp.DOMAIN_IDS_2D.BED_ID |
    bp.DOMAIN_IDS_2D.EAST_ID |
    bp.DOMAIN_IDS_2D.NORTH_EAST_ID
)
uz_boundary_id = bp.DOMAIN_IDS_2D.BED_ID

# Initialize FEM functions for velocity and height
ux_func = bp.FEMFunction1D(u_mesh_1d)
uz_func = bp.FEMFunction1D(u_mesh_1d)
h0_func = bp.FEMFunction1D(h_mesh_1d)
ac_func = bp.FEMFunction1D(h_mesh_1d)

# Initialize Stokes and Free Surface Problems
psp = bp.pStokesProblem(u_mesh_2d, p_mesh_2d)
fsp = bp.FreeSurfaceProblem(h_mesh_1d, u_mesh_1d)

# Plot initial surface profile
plt.plot(xs_vec, zs_vec)

# Time-stepping loop
for k in range(nt):
    # Assemble and solve the nonlinear pStokes problem
    psp.solve_nonlinear_system_picard(
        A, n_i, eps_reg_2, fssa_version, fssa_param, force_x, force_z,
        ux_boundary_id, uz_boundary_id, max_iter, stol, gp
    )

    # Extract velocity field solutions
    ux_vec = psp.velocity_x().vals
    uz_vec = psp.velocity_z().vals
    ux_func.vals = ux_vec[sdofs_u]
    uz_func.vals = uz_vec[sdofs_u]
    h0_func.vals = zs_vec  # Update initial height values

    # Print surface energy and domain area
    print(f"||E|| = {h0_func.L2_norm()}")
    print(f"A = {u_mesh_2d.area()}")

    # Solve the free surface problem using explicit time stepping
    fsp.assemble_lhs_explicit(gp)
    fsp.commit_lhs()
    fsp.assemble_rhs_explicit(h0_func, ux_func, uz_func, ac_func, dt, gp)
    fsp.solve_linear_system()

    # Update surface elevation
    zs_vec = fsp.zs_vec

    # Update mesh with new surface elevation
    u_mesh_2d.extrude_z(zs_vec)
    p_mesh_2d.extrude_z(zs_vec)

    # Plot updated surface profile
    plt.plot(xs_vec, zs_vec)

# Show final plot
plt.show()
