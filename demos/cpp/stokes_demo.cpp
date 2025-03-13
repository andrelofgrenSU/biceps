#include <enums.hpp>
#include <pstokes_fem.hpp>

#define LENGTH 100.0
#define HEIGHT 100.0
#define AMPLITUDE 0.1

FloatType zb_expr(FloatType x)
{
    return 0.0;
}

FloatType zs_expr(FloatType x)
{
    return HEIGHT + AMPLITUDE*SIN_FUNC(PI_CONST*x/LENGTH);
}

FloatType force_x(FloatType x, FloatType z) {
    return 0.0;
}

FloatType force_z(FloatType x, FloatType z) {
    return 1e-3*ICE_DENSITY*GRAVITY;
}

int main(int argc, char *argv[])
{
    int nx = 10;
    int nz = 2;
    int deg_u = 2;
    int deg_p = 1;
    int cell_type = MESH2D::TRIANGLE_LEFT;
    int gp = 4;

    FloatType eta = 1e12 * PA_TO_MPA / SEC_PER_YEAR;
    FloatType A = 1/(2*eta);
    FloatType n_i = 1.0;
    FloatType eps_reg_2 = 1e-10;
    StructuredMesh u_mesh_2d(nx, nz, deg_u, cell_type);
    StructuredMesh p_mesh_2d(nx, nz, deg_p, cell_type);
    pStokesProblem psp(u_mesh_2d, p_mesh_2d);
    psp.assemble_stress_block(A, n_i, eps_reg_2, gp);
    psp.assemble_incomp_block(gp);
    psp.assemble_rhs_vec(force_x, force_z, gp);
    psp.commit_lhs_mat();
    int ux_boundary_id = MESH2D::NORTH_WEST_ID | MESH2D::WEST_ID | MESH2D::BED_ID | MESH2D::EAST_ID | MESH2D::NORTH_EAST_ID;
    int uz_boundary_id = MESH2D::BED_ID;
    psp.apply_zero_dirichlet_bc(ux_boundary_id, uz_boundary_id);
    psp.solve_linear_system();
    std::cout << psp.rhs_vec << std::endl;

    return 0;
}
