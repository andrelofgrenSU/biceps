#include <enums.hpp>
#include <pstokes_fem.hpp>
#include <free_surface_fem.hpp>

#define GRAVITY 9.8
#define ICE_DENSITY 910

// Define domain and grid parameters
FloatType x0 = 0.0;  // Left end
FloatType x1 = 100.0;  // Right end
FloatType L = x1 - x0;  // Length of the domain
FloatType H = 1.0;  // Mean height of the domain
FloatType z0 = 0.5;  // Amplitude of surface undulation

FloatType A = 100.0  // Ice softness parameter
FloatType n_i = 3.0;  // Glen exponent
FloatType eps_reg_2 = 1e-10;  // Regularization parameter
int fssa_version = FSSA_NONE  // No FSSA stabilization
FloatType fssa_param = 0  // Stabilization parameter in FSSA 

int nx = 50;  // Number of elements in x-direction
int nz = 5;  // Number of elements in z-direction
int nt = 10;  // Number of time steps
FloatType dt = 1.0;  // Time step size
int deg_u = 2;  // Polynomial degree for velocity field
int deg_p = 1;  // Polynomial degree for pressure field
int deg_h = 1;  // Polynomial degree for height field
int gauss_precision = 5;  // Number of Gauss points in each direction per element
int max_iter = 100;  // Maximum number of iterations for solver
FloatType stol = 1e-10;  // Convergence tolerance for solver
int cell_type = MESH2D::TRIANGLE_LEFT; 

FloatType zb_expr(FloatType x)
{
    return 0.0;
}

FloatType zs_expr(FloatType x)
{
    return H + z0*SIN_FUNC(PI_CONST*x/L);
}

FloatType force_x(FloatType x, FloatType z) {
    return 0.0;
}

FloatType force_z(FloatType x, FloatType z) {
    return 1e-3*ICE_DENSITY*GRAVITY;
}

int main(int argc, char *argv[])
{
    // Create structured meshes for velocity, pressure, and height fields
    StructuredMesh u_mesh_2d(nx, nz, deg_u, cell_type);
    StructuredMesh p_mesh_2d(nx, nz, deg_p, cell_type);

    // Extrude the meshes in the x and z directions
    u_mesh_2d.extrude_x(x0, x1);
    u_mesh_2d.extrude_z(zb_expr, zs_expr);
    p_mesh_2d.extrude_x(x0, x1);
    p_mesh_2d.extrude_z(zb_expr, zs_expr);

    // Extract degrees of freedom for surface nodes
    std::vector<int> sdofs_u = u_mesh_2d.extract_dof_inds(MESH2D::SURFACE_ID);
    std::vector<int> sdofs_h = u_mesh_2d.extract_vertex_dof_inds(MESH2D::SURFACE_ID);

    // Extract surface coordinates
    Eigen::MatrixX<FloatType> spmat_u = u_mesh_2d.pmat(sdofs_u, Eigen::all);
    Eigen::MatrixX<FloatType> spmat_h = u_mesh_2d.pmat(sdofs_h, Eigen::all);
    // Project mesh to z=0
    spmat_u(Eigen::all, 1).array() = 0.0;
    spmat_h(Eigen::all, 1).array() = 0.0;

    // Create 1D meshes for velocity and height on the surface
    IntervalMesh u_mesh_1d = IntervalMesh(spmat_u, deg_u);
    IntervalMesh h_mesh_1d = IntervalMesh(spmat_h, deg_h);

    // Define ids for Dirichlet boundary condition
    int ux_boundary_id = (
        MESH2D::NORTH_WEST_ID |
        MESH2D::WEST_ID |
        MESH2D::BED_ID |
        MESH2D::EAST_ID |
        MESH2D::NORTH_EAST_ID
    );
    int uz_boundary_id = MESH2D::BED_ID;

    // Initialize FEM functions for height, velocity, and accumulation
    FEMFunction1D ux_func = FEMFunction1D(u_mesh_1d);
    FEMFunction1D uz_func = FEMFunction1D(u_mesh_1d);
    FEMFunction1D h0_func = FEMFunction1D(h_mesh_1d);
    FEMFunction1D ac_func = FEMFunction1D(h_mesh_1d);

    // Initialize the pStokes and Free Surface Problem
    pStokesProblem psp(u_mesh_2d, p_mesh_2d);
    FreeSurfaceProblem fsp(h_mesh_1d, u_mesh_1d);

    for (int k = 0; k < nt; k++) {
        psp.solve_nonlinear_system_picard(
            A, n_i, eps_reg_2, fssa_version, fssa_param, force_x, force_z,
            ux_boundary_id, uz_boundary_id, max_iter, stol, gauss_precision
        );

        // To be continued...
        // Requires installing matplot++ and gnuplot
    }

    return 0;
}
