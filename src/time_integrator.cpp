#include <time_integrator.hpp>

TimeIntegrator::TimeIntegrator(
    pStokesProblem &psp, FreeSurfaceProblem &fsp
) : psp(psp), fsp(fsp)
{
    svinds = psp.u_mesh.extract_vertex_dof_inds(MESH2D::SURFACE_ID);
}

Eigen::VectorX<FloatType> TimeIntegrator::step_explicit(FloatType dt)
{
    psp.solve_nonlinear_system();
    psp.reset_system();

    FEMFunction1D
        h0_fem_func(fsp.h_mesh),
        ux_fem_func(fsp.u_mesh),
        uz_fem_func(fsp.u_mesh),
        ac_fem_func(fsp.h_mesh);

    FEMFunction2D ux_fem_func_p2 = psp.velocity_x();
    FEMFunction2D uz_fem_func_p2 = psp.velocity_z();
    h0_fem_func.vals =
        ux_fem_func_p2.extract_vertex_subvec(MESH2D::SURFACE_ID)(Eigen::all, 1);
    ux_fem_func.vals =
        ux_fem_func_p2.extract_dof_subvec(MESH2D::SURFACE_ID)(Eigen::all, 2);
    uz_fem_func.vals =
        uz_fem_func_p2.extract_dof_subvec(MESH2D::SURFACE_ID)(Eigen::all, 2);

    fsp.assemble_lhs_explicit();
    fsp.commit_lhs();
    fsp.assemble_rhs_explicit(
        h0_fem_func, ux_fem_func, uz_fem_func, ac_fem_func, dt
    );
    fsp.solve_linear_system();
    fsp.reset_system();

    Eigen::VectorX<FloatType> dh_vec = fsp.zs_vec - h0_fem_func.vals;
    return dh_vec;
}

void TimeIntegrator::extrude_mesh_z(const Eigen::VectorX<FloatType> &zs_vec)
{
    psp.u_mesh.extrude_z(zs_vec);
    psp.p_mesh.extrude_z(zs_vec);
}
