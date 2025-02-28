#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <biceps.hpp>

Biceps::Biceps(
	long double A,
	long double n_i,
	long double eps_reg_2,
	std::function<long double(long double, long double)> force_x,
	std::function<long double(long double, long double)> force_z,
	std::function<long double(long double, long double)> ac_expr,
	std::function<long double(long double)> surf_expr,
	std::function<long double(long double)> bed_expr,
	long double x0,
	long double x1, 
	int nx,
	int nz,
	int cell_type,
	int mixed_element,
	std::string output_dirname
) : A(A), n_i(n_i), eps_reg_2(eps_reg_2), force_x(force_x), force_z(force_z), io_handler(output_dirname)
{
	if (mixed_element == P1P1) {
		psp = pStokesProblem(nx, nz, 1, 1, cell_type);
	} else if (mixed_element == P2P1) {
		psp = pStokesProblem(nx, nz, 2, 1, cell_type);
	} else {
		print_msg((boost::format{"Element with id %d is not defined."} %mixed_element).str(), ERROR);
		exit(-1);
	}

	extrude_mesh_x(x0, x1);
	extrude_mesh_z(surf_expr, bed_expr);

	fsp = FreeSurfaceProblem(x0, x1, nx, 2, 1);
	ac_fem_func = FEMFunction1D(fsp.h_mesh);
	ac_fem_func.assign(ac_expr);
	svinds = psp.u_mesh.extract_vertex_dof_inds(SURFACE_ID);
	bvinds = psp.u_mesh.extract_vertex_dof_inds(BED_ID);
}

// TODO: Picard solver should ideally be defined inside the pStokes class
void Biceps::assemble_and_solve_pstokes(
	int fssa_version,
	long double fssa_param,
	long double picard_rtol,
	int picard_max_iter,
	int gp,
	bool reset_system
) {
	Eigen::SparseMatrix<long double> M = FEM2D::assemble_mass_matrix(psp.u_mesh, gp);
	psp.assemble_incomp_block(gp);
	if (fssa_version == FSSA_VERTICAL)
		psp.assemble_fssa_vertical_block(force_x, force_z, fssa_param, gp);
	else if (fssa_version == FSSA_NORMAL)
		psp.assemble_fssa_normal_block(force_z, fssa_param, gp);
	psp.assemble_rhs_vec(force_x, force_z, gp);
	int nof_pushed_elements = psp.nof_pushed_elements;

	long double error_abs;
	FEMFunction2D ux_old = psp.velocity_x();
	FEMFunction2D uz_old = psp.velocity_z();
	FEMFunction2D ux = FEMFunction2D(psp.u_mesh);
	FEMFunction2D uz = FEMFunction2D(psp.u_mesh);
	long double area_sqrt = sqrt(FEM2D::calculate_area(psp.u_mesh, gp));
	// Nonlinear iterations (Picard)
	int ip = 0;
	do {
		psp.lhs_coeffs.resize(nof_pushed_elements);
		psp.nof_pushed_elements = nof_pushed_elements;
		psp.assemble_stress_block(A, n_i, eps_reg_2, gp);
		psp.commit_lhs_mat();

		psp.apply_standard_dirichlet_bc();
		psp.prune_lhs(1e-12, 1);
		psp.solve_linear("umfpack");

		// Estimate the error
		FEMFunction2D ex_func = ux - ux_old;
		FEMFunction2D ez_func = uz - uz_old;
		long double ex = FEM2D::L2_norm(ex_func, M);
		long double ez = FEM2D::L2_norm(ez_func, M);
		error_abs = sqrt(ex*ex + ez*ez)/area_sqrt;
		print_msg((boost::format{"<Picard> iteration %d: e(abs) = %g"} %(ip+1) %error_abs).str(), INFO);

		ux_old.vals = ux.vals;
		uz_old.vals = uz.vals;
		ip++;
	} while (error_abs > picard_rtol && ip < picard_max_iter);

	print_ps_summary(ip, error_abs);

	if (reset_system)
		psp.reset_system();
}

Eigen::MatrixX<long double> Biceps::assemble_evolution_matrix(
	int fssa_version,
	long double fssa_param,
	long double picard_rtol,
	int picard_max_iter,
	int gp,
	bool reset_system
) {
	Eigen::MatrixX<long double> L;
	std::vector<int> bed_dinds = psp.u_mesh.extract_dof_inds(BED_ID);
	std::vector<int> surf_vinds = psp.u_mesh.extract_vertex_dof_inds(SURFACE_ID);

	Eigen::SparseMatrix<long double> M = FEM2D::assemble_mass_matrix(psp.u_mesh, gp);
	psp.assemble_incomp_block(gp);
	if (fssa_version == FSSA_VERTICAL)
		psp.assemble_fssa_vertical_block(force_x, force_z, fssa_param, gp);
	else if (fssa_version == FSSA_NORMAL)
		psp.assemble_fssa_normal_block(force_z, fssa_param, gp);
	psp.assemble_rhs_vec(force_x, force_z, gp);
	int nof_pushed_elements = psp.nof_pushed_elements;

	long double error_rel;
	int ip = 0;
	FEMFunction2D u_old(psp.u_mesh);
	// Nonlinear iterations (Picard)
	do {
		psp.lhs_coeffs.resize(nof_pushed_elements);
		psp.nof_pushed_elements = nof_pushed_elements;
		psp.assemble_stress_block(A, n_i, eps_reg_2, gp);
		psp.commit_lhs_mat();

		psp.apply_dirichlet_bc(
			 NORTH_WEST_ID | WEST_ID | NORTH_EAST_ID | EAST_ID, HORIZONTAL, 0.0
		);
		psp.apply_noslip_bc(BED_ID);
		psp.prune_lhs(1e-12, 1);

		psp.solve_linear("umfpack");

		// Estimate the error
		FEMFunction2D ux = psp.velocity_x();
		FEMFunction2D uz = psp.velocity_z();
		FEMFunction2D u(psp.u_mesh);
		u.vals = (
			ux.vals.cwiseProduct(ux.vals) + uz.vals.cwiseProduct(uz.vals)
		).cwiseSqrt();
		FEMFunction2D eu = u - u_old;
		error_rel = FEM2D::L2_norm(eu, M)/FEM2D::L2_norm(u, M);
		u_old = u;
		ip++;
	} while (error_rel > picard_rtol && ip < picard_max_iter);

	print_ps_summary(ip, error_rel);

	// Assemble the propagation matrix
	Eigen::MatrixX<long double> K = psp.lhs_mat.toDense();
	Eigen::MatrixX<long double> K_inv = K.inverse();
	Eigen::MatrixX<long double> Gzz = K_inv(psp.uz_v2d, psp.uz_v2d);
	Eigen::MatrixX<long double> Gzz_p1_s = Gzz(surf_vinds, Eigen::all);
	Eigen::MatrixX<long double> P_plus = FEM2D::assemble_expansion_matrix(psp.u_mesh, force_z, gp);
	P_plus(bed_dinds, Eigen::all).setZero();
	L = Gzz_p1_s*P_plus;
	if (reset_system)
		psp.reset_system();
    return L;
}

Eigen::MatrixX<long double> Biceps::assemble_jacobian_matrix(
	long double dh,
	int fssa_version,
	long double fssa_param,
	long double picard_rtol,
	int picard_max_iter,
	int gp_ps,
	int gp_fs
) {
	int nx = psp.u_mesh.nx;
	Eigen::VectorX<long double> h0_vec = psp.u_mesh.pmat(svinds, 1);
	Eigen::MatrixX<long double> dF = Eigen::MatrixX<long double>::Zero(nx+1, nx+1);
	Eigen::VectorX<long double> e = Eigen::VectorX<long double>::Zero(nx+1);
	Eigen::VectorX<long double> un_0_vec = step_rk1(1.0, fssa_version, fssa_param, picard_rtol, picard_max_iter, gp_ps, gp_fs, true);
	for (int i = 0; i < nx+1; i++) {
		e(i) = 1;
		Eigen::VectorX<long double> h1_vec = h0_vec + dh*e;
		extrude_mesh_z(h1_vec);
		Eigen::VectorX<long double> un_1_vec = step_rk1(1.0, fssa_version, fssa_param, picard_rtol, picard_max_iter, gp_ps, gp_fs, true);
		dF(i, Eigen::all) = (un_1_vec - un_0_vec)/dh;
		e(i) = 0;
	}
	// Reset to initial mesh
	extrude_mesh_z(h0_vec);
	return dF;
}

Eigen::MatrixX<long double> Biceps::assemble_error_evolution_matrix(
	long double dh,
	int fssa_version,
	long double fssa_param,
	long double picard_rtol,
	int picard_max_iter,
	int gp,
	bool reset_system 
) {
	int nx = psp.u_mesh.nx;
	Eigen::VectorX<long double> h0_vec = psp.u_mesh.pmat(svinds, 1);
	Eigen::MatrixX<long double> DL = Eigen::MatrixX<long double>::Zero(nx+1, nx+1);
	Eigen::VectorX<long double> e = Eigen::VectorX<long double>::Zero(nx+1);
	Eigen::MatrixX<long double> L0 = assemble_evolution_matrix(fssa_version, fssa_param, 1e-8, 1, 4, true);
	for (int i = 0; i < nx+1; i++) {
		e(i) = dh;
		Eigen::VectorX<long double> h1_vec = h0_vec + e;
		extrude_mesh_z(h1_vec);
		Eigen::MatrixX<long double> L1 = assemble_evolution_matrix(fssa_version, fssa_param, 1e-8, 1, 4, true);
		DL = DL + (L1 - L0)/dh;
		e(i) = 0;
	}
	// Reset to initial mesh
	extrude_mesh_z(h0_vec);
	return L0 + DL;
}

Eigen::VectorX<long double> Biceps::step_rk1(
	long double dt,
	int fssa_version,
	long double fssa_param,
	long double picard_ps_rtol,
	int picard_ps_max_iter,
	int gp_ps,
	int gp_fs,
	bool reset_system
) {
	assemble_and_solve_pstokes(
		fssa_version, fssa_param, picard_ps_rtol, picard_ps_max_iter, gp_ps, reset_system
	);

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

	fsp.assemble_lhs_explicit_euler(gp_fs);
	fsp.assemble_rhs_explicit_euler(
		h0_fem_func, ux_fem_func, uz_fem_func, ac_fem_func, dt, gp_fs
	);
	fsp.commit_lhs();
	fsp.solve("umfpack");

	Eigen::VectorX<long double> dh_vec = fsp.zs_sol - h0_fem_func.vals;
	return dh_vec;
}

Eigen::VectorX<long double> Biceps::step_rk2(
	long double dt,
	long double picard_ps_rtol,
	int picard_ps_max_iter,
	int gp_ps,
	int gp_fs
) {
	Eigen::VectorX<long double> h0_vec = psp.u_mesh.pmat(svinds, 1);
	Eigen::VectorX<long double> k1 = step_rk1(
		dt, FSSA_NONE, 0.0, picard_ps_rtol, picard_ps_max_iter, gp_ps, gp_fs, true
	);
	Eigen::VectorX<long double> h1_vec = h0_vec + k1;

	extrude_mesh_z(h1_vec);
	Eigen::VectorX<long double> k2 = step_rk1(
		dt, FSSA_NONE, 0.0, picard_ps_rtol, picard_ps_max_iter, gp_ps, gp_fs, true
	);

	return 0.5*(k1 + k2);
}

Eigen::VectorX<long double> Biceps::step_rk4(
	long double dt,
	long double picard_ps_rtol,
	int picard_ps_max_iter,
	int gp_ps,
	int gp_fs
) {
	Eigen::VectorX<long double> h0_vec = psp.u_mesh.pmat(svinds, 1);

	Eigen::VectorX<long double> k1 = step_rk1(
		dt, FSSA_NONE, 0.0, picard_ps_rtol, picard_ps_max_iter, gp_ps, gp_fs, true
	);
	Eigen::VectorX<long double> h1_vec = h0_vec + 0.5*k1;
	extrude_mesh_z(h1_vec);

	Eigen::VectorX<long double> k2 = step_rk1(
		dt, FSSA_NONE, 0.0, picard_ps_rtol, picard_ps_max_iter, gp_ps, gp_fs, true
	);
	Eigen::VectorX<long double> h2_vec = h0_vec + 0.5*k2;
	extrude_mesh_z(h2_vec);

	Eigen::VectorX<long double> k3 = step_rk1(
		dt, FSSA_NONE, 0.0, picard_ps_rtol, picard_ps_max_iter, gp_ps, gp_fs, true
	);
	Eigen::VectorX<long double> h3_vec = h0_vec + k3;
	extrude_mesh_z(h3_vec);

	Eigen::VectorX<long double> k4 = step_rk1(
		dt, FSSA_NONE, 0.0, picard_ps_rtol, picard_ps_max_iter, gp_ps, gp_fs, true
	);

	return (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
}

Eigen::VectorX<long double> Biceps::step_semi_implicit(
	long double dt,
	int fssa_version,
	long double fssa_param,
	long double picard_ps_rtol,
	int picard_ps_max_iter,
	int gp_ps,
	int gp_fs,
	bool reset_system
) {
	assemble_and_solve_pstokes(
		fssa_version, fssa_param, picard_ps_rtol, picard_ps_max_iter, gp_ps, reset_system
	);

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

	fsp.assemble_lhs_implicit_euler(ux_fem_func, dt, gp_fs);
	fsp.assemble_rhs_implicit_euler(
		h0_fem_func, uz_fem_func, ac_fem_func, dt, gp_fs
	);
	fsp.commit_lhs();
	fsp.solve("umfpack");

	Eigen::VectorX<long double> dh_vec = fsp.zs_sol - h0_fem_func.vals;
	return dh_vec;
}

Eigen::VectorX<long double> Biceps::step_midpoint(
	long double dt,
	long double picard_ps_rtol,
	long double picard_fs_atol,
	int picard_ps_max_iter,
	int picard_fs_max_iter,
	int gp_ps,
	int gp_fs
) {
	std::vector<int> svinds = psp.u_mesh.extract_vertex_dof_inds(SURFACE_ID);
	Eigen::VectorX<long double> svec_mid_fi = psp.u_mesh.pmat(svinds, 1);
	Eigen::VectorX<long double> z0 = svec_mid_fi;
	int p;
	long double error = 0;
	for (p = 0; p < picard_fs_max_iter; p++) {
		Eigen::VectorX<long double> sm_vec = 0.5*(svec_mid_fi + z0);
		extrude_mesh_z(sm_vec);
		assemble_and_solve_pstokes(
			FSSA_NONE, 0.0, picard_ps_rtol, picard_ps_max_iter, gp_ps, true
		);

		FEMFunction1D
			s0_fem_func(fsp.h_mesh),
			ux_fem_func(fsp.u_mesh),
			uz_fem_func(fsp.u_mesh),
			ac_fem_func(fsp.h_mesh);

		FEMFunction2D ux_fem_func_p2 = psp.velocity_x();
		FEMFunction2D uz_fem_func_p2 = psp.velocity_z();
		s0_fem_func.vals = svec_mid_fi;
		ux_fem_func.vals =
			ux_fem_func_p2.extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);
		uz_fem_func.vals =
			uz_fem_func_p2.extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);

		// Solve the free-surface equation for new height
		fsp.assemble_lhs_implicit_midpoint(ux_fem_func, dt, 3);
		fsp.assemble_rhs_implicit_midpoint(s0_fem_func, ux_fem_func, uz_fem_func, ac_fem_func, dt, 3);
		fsp.commit_lhs();
		fsp.solve("umfpack");

		Eigen::VectorX<long double> z1 = fsp.zs_sol;
		error = (z1 - z0).norm();
		z0 = z1;
		if (error < picard_fs_atol) {
			print_msg((boost::format{
				"Midpoint :: Picard iterations converged to specified tolerance %g after %d iterations."
				} %picard_fs_atol %p).str(),
				INFO
			);
			break;
		}
	}

	if (p == picard_fs_max_iter) {
		print_msg((boost::format{
			"Midpoint :: Picard iterations failed to converge to specified tolerance %g (error = %g)."
			} %picard_fs_atol %error).str(),
			WARN
		);
	}
	Eigen::VectorX<long double> dh_vec = z0 - svec_mid_fi;
	return dh_vec;
}

Eigen::VectorX<long double> Biceps::step_crank_nicolson(
	long double dt,
	long double picard_ps_rtol,
	long double picard_fs_atol,
	int picard_ps_max_iter,
	int picard_fs_max_iter,
	int gp_ps,
	int gp_fs
) {
	std::vector<int> svinds = psp.u_mesh.extract_vertex_dof_inds(SURFACE_ID);
	Eigen::VectorX<long double> hvec_cni_fi = psp.u_mesh.pmat(svinds, 1);
	Eigen::VectorX<long double> z0 = hvec_cni_fi;
	assemble_and_solve_pstokes(
		FSSA_NONE, 0.0, picard_ps_rtol, picard_ps_max_iter, gp_ps, true
	);
	Eigen::VectorX<long double> ux_0_vec = psp.velocity_x().extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);
	Eigen::VectorX<long double> uz_0_vec = psp.velocity_z().extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);

	FEMFunction1D
		h0_fem_func(fsp.h_mesh),
		ux_0_fem_func(fsp.u_mesh),
		ux_1_fem_func(fsp.u_mesh),
		uz_av_fem_func(fsp.u_mesh), 
		ac_av_fem_func(fsp.h_mesh);
	h0_fem_func.vals = hvec_cni_fi;
	ux_0_fem_func.vals = ux_0_vec;

	int p;
	long double error = 0;
	for (p = 0; p < picard_fs_max_iter; p++) {
		extrude_mesh_z(z0);
		assemble_and_solve_pstokes(
			FSSA_NONE, 0.0, picard_ps_rtol, picard_ps_max_iter, gp_ps, true
		);

		// Solve the free-surface equation for new height
		Eigen::VectorX<long double> ux_1_vec = psp.velocity_x().extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);
		Eigen::VectorX<long double> uz_1_vec = psp.velocity_z().extract_dof_subvec(SURFACE_ID)(Eigen::all, 2);
		ux_1_fem_func.vals = ux_1_vec;
		uz_av_fem_func.vals = 0.5*(uz_0_vec + uz_1_vec);

		fsp.assemble_lhs_crank_nicolson(ux_1_fem_func, dt, 3);
		fsp.assemble_rhs_crank_nicolson(
			h0_fem_func, ux_0_fem_func, uz_av_fem_func, ac_av_fem_func, dt, 3
		);
		fsp.commit_lhs();
		fsp.solve("umfpack");

		Eigen::VectorX<long double> z1 = fsp.zs_sol;
		error = (z1 - z0).norm();
		z0 = z1;
		if (error < picard_fs_atol) {
			print_msg((boost::format{
				"Crank-Nicolson :: Picard iterations converged to specified tolerance %g after %d iterations."
				} %picard_fs_atol %p).str(),
				INFO
			);
			break;
		}
	}
	if (p == picard_fs_max_iter) {
		print_msg((boost::format{
			"Crank-Nicolson :: Picard iterations failed to converge to specified tolerance %g (error = %g)."
		} %picard_fs_atol %error).str(), ERROR);
	}
	Eigen::VectorX<long double> dh_vec = z0 - hvec_cni_fi;
	return dh_vec;
}

void Biceps::extrude_mesh_x(Eigen::VectorX<long double> xvec_p1) {
	psp.u_mesh.extrude_x(xvec_p1);
	psp.p_mesh.extrude_x(xvec_p1);
}

void Biceps::extrude_mesh_x(long double x0, long double x1) {
	psp.u_mesh.extrude_x(x0, x1);
	psp.p_mesh.extrude_x(x0, x1);
}

void Biceps::extrude_mesh_z(Eigen::VectorX<long double> svec_p1) {
	psp.u_mesh.extrude_z(svec_p1);
	psp.p_mesh.extrude_z(svec_p1);
}

void Biceps::extrude_mesh_z(Eigen::VectorX<long double> svec_p1, Eigen::VectorX<long double> bvec_p1) {
	psp.u_mesh.extrude_z(svec_p1, bvec_p1);
	psp.p_mesh.extrude_z(svec_p1, bvec_p1);
}

void Biceps::extrude_mesh_z(
	std::function<long double(long double)> s_expr, std::function<long double(long double)> b_expr
) {
	psp.u_mesh.extrude_z(s_expr, b_expr);
	psp.p_mesh.extrude_z(s_expr, b_expr);
}

void Biceps::print_ps_summary(int iters, long double error_rel)
{
	print_msg("-------------------------------------", INFO, false);
	print_msg("********** pStokes Summary **********", INFO, false);
	print_msg("-------------------------------------", INFO, false);

	print_msg((boost::format{"Picard solver finished in %d iterations"} % iters).str(), INFO, false);
	print_msg((boost::format{"r (rel): %.4g"} % (error_rel)).str(), INFO, false);
	print_msg((boost::format{"min ux, max ux: %.4g, %.4g m/yr"} % (1e3*psp.w_vec(psp.ux_v2d).minCoeff()) % (1e3*psp.w_vec(psp.ux_v2d).maxCoeff())).str(), INFO, false);
	print_msg((boost::format{"min uz, max uz: %.4g, %.4g m/yr"} % (1e3*psp.w_vec(psp.uz_v2d).minCoeff()) % (1e3*psp.w_vec(psp.uz_v2d).maxCoeff())).str(), INFO, false);
	print_msg((boost::format{"min p, max p: %.4g, %.4g MPa"} % (psp.w_vec(psp.p_v2d).minCoeff()) % (psp.w_vec(psp.p_v2d).maxCoeff())).str(), INFO, false);
	print_msg("-------------------------------------", INFO, false);
}

void Biceps::init_xdmf(std::string xdmf_filename) {
	io_handler.xdmf_init_time_series(xdmf_filename);
}

void Biceps::write_xdmf() {
	FEMFunction2D ux_fem_func_p1(psp.p_mesh), uz_fem_func_p1(psp.p_mesh);
	FEMFunction2D ux_fem_func_p2 = psp.velocity_x();
	FEMFunction2D uz_fem_func_p2 = psp.velocity_z();
	FEMFunction2D p_fem_func_p1 = psp.pressure();
	ux_fem_func_p1.vals = 
		ux_fem_func_p2.extract_vertex_subvec(DOMAIN_ID)(Eigen::all, 2);
	uz_fem_func_p1.vals = 
		uz_fem_func_p2.extract_vertex_subvec(DOMAIN_ID)(Eigen::all, 2);
	io_handler.xdmf_prepare_data_insertion(psp.p_mesh, time);
	io_handler.xdmf_insert_data(ux_fem_func_p1, "Velocity X");
	io_handler.xdmf_insert_data(uz_fem_func_p1, "Velocity Z");
	io_handler.xdmf_insert_data(p_fem_func_p1, "Pressure");
	io_handler.xdmf_insert_marker_data(psp.p_mesh);
	io_handler.xdmf_finalize_data_insertion();
}
