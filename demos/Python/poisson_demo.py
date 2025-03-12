import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import biceps as bp

# Define domain and mesh parameters
nx = 16  # Number of elements in x-direction
nz = 16  # Number of elements in z-direction
deg = 1  # Polynomial degree of the finite element basis functions
cell_type = bp.CELL_TYPE_2D.TRIANGLE_LEFT  # Type of element

# Create a structured triangular mesh
mesh = bp.StructuredMesh(nx, nz, deg, cell_type)

# Set up and solve the Poisson problem
pp = bp.PoissonProblem(mesh)
pp.assemble_stiffness_block(lambda x, z: 1, lambda x, z: 1, deg + 2)
pp.commit_lhs_mat()
pp.assemble_force_rhs(lambda x, z: 1.0, deg + 2)
pp.apply_dirichlet_bc(lambda x, z: 0.0, bp.DOMAIN_IDS_2D.BOUNDARY_ID)
pp.solve_linear_system()
solution_vec = pp.solve().vals

# Extract and plot the solution
tri = Triangulation(mesh.pmat[:, 0], mesh.pmat[:, 1], mesh.cmat)
fig, ax = plt.subplots()
plt.colorbar(ax.tricontourf(tri, solution_vec, levels=10))
plt.show()
