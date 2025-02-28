#include <time_integrator.hpp>

TimeIntegrator::TimeIntegrator(
	pStokesProblem &psp, FreeSurfaceProblem &fsp
) : psp(psp), fsp(fsp)
{
	svinds = psp.u_mesh.extract_vertex_dof_inds(SURFACE_ID);
}

Eigen::VectorX<FloatType> TimeIntegrator::step_explicit(FloatType dt, bool reset_pstokes)
{
	psp.solve_picard();

	FEMFunction1D
		h0_fem_func(fsp.h_mesh),
		ux_fem_func(fsp.u_mesh),
		uz_fem_func(fsp.u_mesh);

	FEMFunction2D ux_fem_func_p2 = psp.velocity_x();
	FEMFunction2D uz_fem_func_p2 = psp.velocity_z();
	h0_fem_func.vals =
		ux_fem_func_p2.extract_vertex_subvec(SURFACE_ID)(Eigen::all, 1);
	ux_fem_func.vals =
		ux_fem_func_p2.extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);
	uz_fem_func.vals =
		uz_fem_func_p2.extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);

	fsp.assemble_lhs_explicit();
	fsp.assemble_rhs_explicit(
		h0_fem_func, ux_fem_func, uz_fem_func, dt
	);
	fsp.commit_lhs();
	fsp.solve_linear_system();

	if (reset_pstokes)
		psp.reset_system();

	Eigen::VectorX<FloatType> dh_vec = fsp.zs_vec - h0_fem_func.vals;
	return dh_vec;
}

Eigen::MatrixX<FloatType> TimeIntegrator::transient_jacobian(
	int time_scheme, int jac_stencil, FloatType dt_eps
) {
	int nx = psp.u_mesh.nx;
	// Assemble jacobian
	Eigen::MatrixX<FloatType> dF = Eigen::MatrixX<FloatType>::Zero(nx+1, nx+1);
	Eigen::VectorX<FloatType> e = Eigen::VectorX<FloatType>::Zero(nx+1);
	if (jac_stencil == FORWARD) {
		Eigen::VectorX<FloatType> h0_vec = psp.u_mesh.pmat(svinds, 1);
		Eigen::VectorX<FloatType> un_0_vec;
		if (time_scheme == EXPLICIT)
			un_0_vec = step_explicit(1.0, true);
		else if (time_scheme == SEMI_IMPLICIT)
			un_0_vec = step_semi_implicit(1.0, true);
		else {
			//TODO: throw exception
		}
		for (int i = 0; i < nx+1; i++) {
			e(i) = 1.0;
			// Step forward
			Eigen::VectorX<FloatType> h1_vec = h0_vec + dt_eps*e;
			extrude_mesh_z(h1_vec);
			Eigen::VectorX<FloatType> un_1_vec;
			if (time_scheme == EXPLICIT)
				un_1_vec = step_explicit(1.0, true);
			else if (time_scheme == SEMI_IMPLICIT)
				un_0_vec = step_semi_implicit(1.0, true);
			else {
				//TODO: throw exception
			}
			// Calculate difference
			dF(i, Eigen::all) = (un_1_vec - un_0_vec)/dt_eps;
			e(i) = 0.0;
		}
	} else if (jac_stencil == CENTER) {
		Eigen::VectorX<FloatType> h1_vec = psp.u_mesh.pmat(svinds, 1);
		for (int i = 0; i < nx+1; i++) {
			// Step backward
			e(i) = 1.0;
			Eigen::VectorX<FloatType> h0_vec = h1_vec - dt_eps*e;
			extrude_mesh_z(h0_vec);
			Eigen::VectorX<FloatType> un_0_vec;
			if (time_scheme == EXPLICIT)
				un_0_vec = step_explicit(1.0, true);
			else if (time_scheme == SEMI_IMPLICIT)
				un_0_vec = step_semi_implicit(1.0, true);
			else {
				//TODO: throw exception
			}
			// Step forward
			Eigen::VectorX<FloatType> h2_vec = h1_vec + dt_eps*e;
			extrude_mesh_z(h2_vec);
			Eigen::VectorX<FloatType> un_2_vec;
			if (time_scheme == EXPLICIT)
				un_2_vec = step_explicit(1.0, true);
			else if (time_scheme == SEMI_IMPLICIT)
				un_2_vec = step_semi_implicit(1.0, true);
			else {
				//TODO: throw exception
			}
			// Calculate difference
			dF(i, Eigen::all) = 0.5*(un_2_vec - un_0_vec)/dt_eps;
			e(i) = 0.0;
		}
	} else {
		//TODO: throw exception
	}
	return dF;
}

Eigen::VectorX<FloatType> TimeIntegrator::step_semi_implicit(FloatType dt, bool reset_pstokes)
{
	psp.solve_picard();

	FEMFunction1D
		h0_fem_func(fsp.h_mesh),
		ux_fem_func(fsp.u_mesh),
		uz_fem_func(fsp.u_mesh);

	FEMFunction2D ux_fem_func_p2 = psp.velocity_x();
	FEMFunction2D uz_fem_func_p2 = psp.velocity_z();
	h0_fem_func.vals =
		ux_fem_func_p2.extract_vertex_subvec(SURFACE_ID)(Eigen::all, 1);
	ux_fem_func.vals =
		ux_fem_func_p2.extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);
	uz_fem_func.vals =
		uz_fem_func_p2.extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);

	fsp.assemble_lhs_semi_implicit(ux_fem_func, dt);
	fsp.assemble_rhs_semi_implicit(
		h0_fem_func, uz_fem_func, dt
	);
	fsp.commit_lhs();
	fsp.solve_linear_system();

	if (reset_pstokes)
		psp.reset_system();

	Eigen::VectorX<FloatType> dh_vec = fsp.zs_vec - h0_fem_func.vals;
	return dh_vec;
}

Eigen::VectorX<FloatType> TimeIntegrator::step_rk2(FloatType dt, bool reset_pstokes)
{
	Eigen::VectorX<FloatType> h0_vec = psp.u_mesh.pmat(svinds, 1);
	Eigen::VectorX<FloatType> k1 = step_explicit(dt, true);
	Eigen::VectorX<FloatType> h1_vec = h0_vec + k1;

	extrude_mesh_z(h1_vec);
	Eigen::VectorX<FloatType> k2 = step_explicit(dt, reset_pstokes);

	return 0.5*(k1 + k2);
}

void TimeIntegrator::extrude_mesh_z(const Eigen::VectorX<FloatType> &zs_vec)
{
	psp.u_mesh.extrude_z(zs_vec);
	psp.p_mesh.extrude_z(zs_vec);
}
