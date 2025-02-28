#include <pstokes_fem.hpp>

pStokesProblem::pStokesProblem(
	int nx, int nz, int deg_u, int deg_p, int cell_type
) : u_mesh(nx, nz, deg_u, cell_type), p_mesh(nx, nz, deg_p, cell_type) 
{
	nv_dofs = u_mesh.nof_dofs();
	np_dofs = p_mesh.nof_dofs();
	n_dofs = 2*nv_dofs + np_dofs; 

	rhs_vec = Eigen::VectorX<long double>::Zero(n_dofs);

	ux_v2d = Eigen::VectorXi::Zero(nv_dofs);
	uz_v2d = Eigen::VectorXi::Zero(nv_dofs);	
	u_v2d = Eigen::VectorXi::Zero(2*nv_dofs);
	p_v2d = Eigen::VectorXi::Zero(np_dofs);

	ux_d2v = Eigen::VectorXi::Constant(n_dofs, -1);
	uz_d2v = Eigen::VectorXi::Constant(n_dofs, -1);

	int dof = 0;
	for (int vi = 0; vi < nv_dofs; vi++) {
		ux_v2d(vi) = dof;
		ux_d2v(dof) = vi;
		dof += 2;
	}
	dof = 1;
	for (int vi = 0; vi < nv_dofs; vi++) {
		uz_v2d(vi) = dof;
		uz_d2v(dof) = vi;
		dof += 2;
	}
	dof = 2*nv_dofs;
	for (int vi = 0; vi < np_dofs; vi++) {
		p_v2d(vi) = dof;
		dof += 1;
	}
	u_v2d(Eigen::seqN(0, nv_dofs, 2)) = ux_v2d;
	u_v2d(Eigen::seq(1, 2*nv_dofs, 2)) = uz_v2d;

	w_vec = Eigen::VectorX<long double>::Zero(n_dofs);

	lhs_mat = Eigen::SparseMatrix<long double>(n_dofs, n_dofs);
	// FIXME: What is an appropriate number of nnz?
	lhs_coeffs.reserve(12*u_mesh.lhs_nnz());
	reset_lhs();
}

// FIXME: Is there a more efficient way of inserting elements directly into the matrix?
void pStokesProblem::assemble_stress_block(
	long double A, long double n_i, long double eps_reg_2, int gp
) {
	Eigen::MatrixX<long double> node_coords_v, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<long double> qweights, detJ_rs;
	Eigen::VectorXi element_v;
	
	FEM2D::gauss_legendre_quadrature(
		gp, u_mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		u_mesh.degree(),
		u_mesh.cell_type(),
		qpoints_rs,
		phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<long double>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<long double>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<long double>::Zero(2*qpoints_rs.rows(), u_mesh.dofs_per_cell());

	Eigen::MatrixX<long double> A_xx = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_cell(), u_mesh.dofs_per_cell()
	);
	Eigen::MatrixX<long double> A_xz = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_cell(), u_mesh.dofs_per_cell()
	);
	Eigen::MatrixX<long double> A_zz = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_cell(), u_mesh.dofs_per_cell()
	);

	Eigen::MatrixX<long double> A_block = Eigen::MatrixX<long double>::Zero(
		2*u_mesh.dofs_per_cell(), 2*u_mesh.dofs_per_cell()
	);

	for (int k = 0; k < u_mesh.cmat.rows(); k++) {
		element_v = u_mesh.cmat.row(k);
		node_coords_v = u_mesh.pmat(element_v, Eigen::all);
		FEM2D::map_to_reference_cell(
			u_mesh.degree(), u_mesh.cell_type(), node_coords_v, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		Eigen::VectorX<long double> ux_vec = w_vec(ux_v2d(element_v));
		Eigen::VectorX<long double> uz_vec = w_vec(uz_v2d(element_v));
		Eigen::MatrixX<long double> Sx = dphi_xz(Eigen::seq(0, Eigen::last, 2), Eigen::all);
		Eigen::MatrixX<long double> Sz = dphi_xz(Eigen::seq(1, Eigen::last, 2), Eigen::all);
		Eigen::VectorX<long double> ddx_ux = Sx*ux_vec;
		Eigen::VectorX<long double> ddz_ux = Sz*ux_vec;
		Eigen::VectorX<long double> ddx_uz = Sx*uz_vec;
		Eigen::VectorX<long double> ddz_uz = Sz*uz_vec;

		for (int q = 0; q < qpoints_rs.rows(); q++) {
			int w = 2*q;
			long double eff_strain_rate_2 = ddx_ux(q)*ddx_ux(q) + 0.25*(
				ddz_ux(q)*ddz_ux(q)
				+ 2.0*ddz_ux(q)*ddx_uz(q)
				+ ddx_uz(q)*ddx_uz(q)
			);
			// Glen's flow law rheology
			long double eta = glen_flow_viscosity(
				A, n_i, eps_reg_2, eff_strain_rate_2
			);
			long double detJxWxEta = fabs(detJ_rs(q))*qweights(q)*eta;
			// TODO: tabulate dphi*dphi and calculate in advance 
			for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < u_mesh.dofs_per_cell(); j++) {
					A_xx(i, j) += dphi_xz(w, i)*dphi_xz(w, j)*detJxWxEta;
					A_xz(i, j) += dphi_xz(w+1, i)*dphi_xz(w, j)*detJxWxEta;
					A_zz(i, j) += dphi_xz(w+1, i)*dphi_xz(w+1, j)*detJxWxEta;
				}
			}
		}

		for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < u_mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(TripletXd(
					ux_v2d(element_v(i)),
					ux_v2d(element_v(j)),
					2.0*A_xx(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					ux_v2d(element_v(i)),
					ux_v2d(element_v(j)),
					A_zz(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					ux_v2d(element_v(i)),
					uz_v2d[element_v(j)],
					A_xz(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					uz_v2d(element_v(i)),
					ux_v2d(element_v(j)),
					A_xz(j, i)
				));
				lhs_coeffs.push_back(TripletXd(
					uz_v2d(element_v(i)),
					uz_v2d(element_v(j)),
					A_xx(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					uz_v2d(element_v(i)),
					uz_v2d(element_v(j)),
					2.0*A_zz(i, j)
				));
			}
		}
		
		A_xx.setZero();
		A_xz.setZero();
		A_zz.setZero();
	}
	nof_pushed_elements = lhs_coeffs.size();
}

void pStokesProblem::assemble_incomp_block(int gp)
{
	Eigen::MatrixX<long double> node_coords_v, node_coords_p, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz, psi_rs, dpsi_rs, dpsi_xz;
	Eigen::VectorX<long double> qweights, detJ_rs;
	Eigen::VectorXi element_v, element_p;
	
	FEM2D::gauss_legendre_quadrature(
		gp, u_mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		p_mesh.degree(), p_mesh.cell_type(), qpoints_rs, psi_rs, dpsi_rs
	);

	FEM2D::lagrange_basis(
		u_mesh.degree(), u_mesh.cell_type(), qpoints_rs, phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<long double>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<long double>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<long double>::Zero(2*qpoints_rs.rows(), u_mesh.dofs_per_cell());
	dpsi_xz = Eigen::MatrixX<long double>::Zero(2*qpoints_rs.rows(), p_mesh.dofs_per_cell());

	Eigen::MatrixX<long double> B_xx = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_cell(), p_mesh.dofs_per_cell()
	);
	Eigen::MatrixX<long double> B_zz = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_cell(), p_mesh.dofs_per_cell()
	);

	for (int k = 0; k < u_mesh.cmat.rows(); k++) {
		element_v = u_mesh.cmat.row(k);
		element_p = p_mesh.cmat.row(k);
		node_coords_v = u_mesh.pmat(element_v, Eigen::all);
		node_coords_p = p_mesh.pmat(element_p, Eigen::all);

		// TODO: only psi_rs is need; speed up by discard calculating everything else
		FEM2D::map_to_reference_cell(
			p_mesh.degree(), p_mesh.cell_type(), node_coords_p, qpoints_rs,
			psi_rs, dpsi_rs, detJ_rs, qpoints_xz, dpsi_xz
		);

		FEM2D::map_to_reference_cell(
			u_mesh.degree(), u_mesh.cell_type(), node_coords_v, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		for (int q = 0; q < qpoints_rs.rows(); q++) {
			int w = 2*q;
			long double detJxW = fabs(detJ_rs(q))*qweights(q);
			for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < p_mesh.dofs_per_cell(); j++) {
					B_xx(i, j) -= dphi_xz(w, i)*psi_rs(q, j)*detJxW;
					B_zz(i, j) -= dphi_xz(w+1, i)*psi_rs(q, j)*detJxW;
				}
			}
		}

		for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < p_mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(TripletXd(
					ux_v2d(element_v(i)),
					p_v2d(element_p(j)),
					B_xx(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					uz_v2d(element_v(i)),
					p_v2d(element_p(j)),
					B_zz(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					p_v2d(element_p(j)),
					ux_v2d(element_v(i)),
					B_xx(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					p_v2d(element_p(j)),
					uz_v2d(element_v(i)),
					B_zz(i, j)
				));
			}
		}

		B_xx.setZero();
		B_zz.setZero();
	}
	nof_pushed_elements = lhs_coeffs.size();
}

// TODO: add fx component to the weak form
void pStokesProblem::assemble_fssa_block(
	std::function<long double(long double, long double)> force_x,
	std::function<long double(long double, long double)> force_z,
	long double dt, long double theta, int gp
) {
	Eigen::MatrixX<long double> node_coords, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<long double> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi edge_vi;

	FEM1D::gauss_legendre_quadrature(
		gp, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		u_mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<long double>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<long double>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<long double>::Zero(
		2*qpoints_r.rows(), u_mesh.dofs_per_edge()
	);

	Eigen::MatrixX<long double> A_xz = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);
	Eigen::MatrixX<long double> A_zz = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);

	std::vector<int> surf_edge_inds = u_mesh.extract_edge_inds(SURFACE_ID);
	for (int si: surf_edge_inds) {
		edge_vi = u_mesh.emat(si, Eigen::all);
		node_coords = u_mesh.pmat(edge_vi, Eigen::all);
		FEM1D::map_to_reference_cell(
			u_mesh.degree(), node_coords, qpoints_r,
			phi_r, dphi_r, detJ_r, qpoints_x, dphi_x
		);

		long double nx = u_mesh.edge_normals(si, 0);
		long double nz = u_mesh.edge_normals(si, 1);
		for (int q = 0; q < qpoints_r.rows(); q++) {
			long double fx = force_x(qpoints_x(q, 0), qpoints_x(q, 1));
			long double fz = force_z(qpoints_x(q, 0), qpoints_x(q, 1));
			long double dFxW = fabs(detJ_r(q))*qweights(q);
			// TODO: optimize by precalculating thetaxdtxfzx...
			for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
				for (int j = 0; j < u_mesh.dofs_per_edge(); j++) {
					A_xz(i, j) -= theta*dt*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nx;
					A_zz(i, j) -= theta*dt*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nz;
				}
			}
		}

		for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
			for (int j = 0; j < u_mesh.dofs_per_edge(); j++) {
				lhs_coeffs.push_back(TripletXd(
					uz_v2d(edge_vi(i)),
					ux_v2d(edge_vi(j)),
					A_xz(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					uz_v2d(edge_vi(i)),
					uz_v2d(edge_vi(j)),
					A_zz(i, j)
				));
			}
		}

		A_xz.setZero();
		A_zz.setZero();
	}
	nof_pushed_elements = lhs_coeffs.size();
}

void pStokesProblem::assemble_fssa_normal_block(
	std::function<long double(long double, long double)> force_z,
	long double dt, long double theta, int gp
) {
	Eigen::MatrixX<long double> node_coords, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<long double> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi edge_vi;

	FEM1D::gauss_legendre_quadrature(
		gp, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		u_mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<long double>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<long double>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<long double>::Zero(
		2*qpoints_r.rows(), u_mesh.dofs_per_edge()
	);

	Eigen::MatrixX<long double> A_xx = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);
	Eigen::MatrixX<long double> A_xz = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);
	Eigen::MatrixX<long double> A_zz = Eigen::MatrixX<long double>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);

	std::vector<int> surf_edge_inds = u_mesh.extract_edge_inds(SURFACE_ID);
	for (int si: surf_edge_inds) {
		edge_vi = u_mesh.emat(si, Eigen::all);
		node_coords = u_mesh.pmat(edge_vi, Eigen::all);
		FEM1D::map_to_reference_cell(
			u_mesh.degree(), node_coords, qpoints_r,
			phi_r, dphi_r, detJ_r, qpoints_x, dphi_x
		);
		long double nx = u_mesh.edge_normals(si, 0);
		long double nz = u_mesh.edge_normals(si, 1);
		nx = nx/sqrt(nz);
		nz = sqrt(nz);
		for (int q = 0; q < qpoints_r.rows(); q++) {
			long double fz = force_z(qpoints_x(q, 0), qpoints_x(q, 1));
			long double dFxW = fabs(detJ_r(q))*qweights(q);
			// TODO: optimize by precalculating thetaxdtxfzx...
			for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
				for (int j = 0; j < u_mesh.dofs_per_edge(); j++) {
					A_xx(i, j) -= theta*dt*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nx*nx;
					A_xz(i, j) -= theta*dt*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nx*nz;
					A_zz(i, j) -= theta*dt*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nz*nz;
				}
			}
		}

		for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
			for (int j = 0; j < u_mesh.dofs_per_edge(); j++) {
				lhs_coeffs.push_back(TripletXd(
					ux_v2d(edge_vi(i)),
					uz_v2d(edge_vi(j)),
					A_xx(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					ux_v2d(edge_vi(i)),
					uz_v2d(edge_vi(j)),
					A_xz(i, j)
				));
				lhs_coeffs.push_back(TripletXd(
					uz_v2d(edge_vi(i)),
					ux_v2d(edge_vi(j)),
					A_xz(j, i)
				));
				lhs_coeffs.push_back(TripletXd(
					uz_v2d(edge_vi(i)),
					uz_v2d(edge_vi(j)),
					A_zz(i, j)
				));
			}
		}

		A_xx.setZero();
		A_xz.setZero();
		A_zz.setZero();
	}
	nof_pushed_elements = lhs_coeffs.size();
}

// TODO: add fx component to the weak form
void pStokesProblem::assemble_fssa_rhs(
	std::function<long double(long double, long double)> force_x,
	std::function<long double(long double, long double)> force_z,
	std::function<long double(long double, long double)> accumulation,
	long double dt, long double theta, int gp
) {
	Eigen::MatrixX<long double> node_coords, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<long double> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi edge_vi;

	FEM1D::gauss_legendre_quadrature(
		gp, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		u_mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<long double>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<long double>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<long double>::Zero(
		2*qpoints_r.rows(), u_mesh.dofs_per_edge()
	);

	std::vector<int> surf_edge_inds = u_mesh.extract_edge_inds(SURFACE_ID);
	for (int si: surf_edge_inds) {
		edge_vi = u_mesh.emat(si, Eigen::all);
		node_coords = u_mesh.pmat(edge_vi, Eigen::all);
		FEM1D::map_to_reference_cell(
			u_mesh.degree(), node_coords, qpoints_r,
			phi_r, dphi_r, detJ_r, qpoints_x, dphi_x
		);
		long double nx = u_mesh.edge_normals(si, 0);
		long double nz = u_mesh.edge_normals(si, 1);
		for (int q = 0; q < qpoints_r.rows(); q++) {
			long double fx = force_x(qpoints_x(q, 0), qpoints_x(q, 1));
			long double fz = force_z(qpoints_x(q, 0), qpoints_x(q, 1));
			long double ac = accumulation(qpoints_x(q, 0), qpoints_x(q, 1));
			long double dFxW = fabs(detJ_r(q))*qweights(q);
			// TODO: optimize by precalculating thetaxdtxfzx...
			for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
				rhs_vec(ux_v2d(edge_vi(i))) += theta*dt*phi_r(q, i)*ac*fx*dFxW*nz;
				rhs_vec(uz_v2d(edge_vi(i))) += theta*dt*phi_r(q, i)*ac*fz*dFxW*nz;
			}
		}
	}
}

void pStokesProblem::assemble_rhs_vec(
	std::function<long double(long double, long double)> fx,
	std::function<long double(long double, long double)> fz,
	int gp
) {
	Eigen::MatrixX<long double> node_coords_v, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<long double> qweights, detJ_rs;
	Eigen::VectorXi element_v;

	FEM2D::gauss_legendre_quadrature(
		gp, u_mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		u_mesh.degree(), u_mesh.cell_type(), qpoints_rs, phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<long double>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<long double>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<long double>::Zero(2*qpoints_rs.rows(), u_mesh.dofs_per_cell());

	for (int k = 0; k < u_mesh.cmat.rows(); k++) {
		element_v = u_mesh.cmat.row(k);
		node_coords_v = u_mesh.pmat(element_v, Eigen::all);
		// FIXME: only phi_rs is need; speed up by discard calculating everything else
		FEM2D::map_to_reference_cell(
			u_mesh.degree(), u_mesh.cell_type(), node_coords_v, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		for (int q = 0; q < qpoints_rs.rows(); q++) {
			long double detJxW = fabs(detJ_rs(q))*qweights(q);
			for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
				rhs_vec(ux_v2d(element_v(i))) += (
					fx(qpoints_xz(q, 0), qpoints_xz(q, 1))*phi_rs(q, i)
				)*detJxW;
				rhs_vec(uz_v2d(element_v(i))) += (
					fz(qpoints_xz(q, 0), qpoints_xz(q, 1))*phi_rs(q, i)
				)*detJxW;
			}
		}
	}
}

// TODO: Make u a function.
// TODO: Apply bc implicitly by modifying the rhs.
// This reduces the linear system size and preserves the symmetric form,
// which means we can utilize the Cholesky factorization (I think...)
// and other linear solvers which requires a symmetric linear system.
void pStokesProblem::apply_dirichlet_bc(
	const int boundary_part, const int velocity_component, long double u
) {
	// Loop over all nnz elements and set off-diagonals corresponding to boundary nodes to zero
	int *outer_end = lhs_mat.outerIndexPtr() + lhs_mat.rows();
	int col = 0;
	for (
		int *outer_p = lhs_mat.outerIndexPtr();
		outer_p != outer_end;
		outer_p++
	) {
		int nnz = *(outer_p + 1) - *outer_p;
		for (int i = 0; i < nnz; i++) {
			int row = *(lhs_mat.innerIndexPtr() + *outer_p + i);
			if (velocity_component == HORIZONTAL && ux_d2v(row) != -1) {
				if (u_mesh.dimat(u_mesh.d2v(ux_d2v(row))) & boundary_part) {
					long double *val_ptr = lhs_mat.valuePtr() + *outer_p + i;
					*val_ptr = 0.0;
				}
			}
			else if (velocity_component == VERTICAL && uz_d2v(row) != -1) {
				if (u_mesh.dimat(u_mesh.d2v(uz_d2v(row))) & boundary_part) {
					long double *val_ptr = lhs_mat.valuePtr() + *outer_p + i;
					*val_ptr = 0.0;
				}
			}
		}
		col++;
	}

	// Loop over all nnz elements and set diagonals corresponding to boundary nodes to one
	col = 0;
	for (
		int *outer_p = lhs_mat.outerIndexPtr();
		outer_p != outer_end;
		outer_p++
	) {
		if (velocity_component == HORIZONTAL && ux_d2v(col) != -1) {
			if (u_mesh.dimat(u_mesh.d2v(ux_d2v(col))) & boundary_part) {
				lhs_mat.coeffRef(col, col) = 1.0;
				rhs_vec(col) = u;  // Dirichlet bc
			}
		}
		else if (velocity_component == VERTICAL && uz_d2v(col) != -1) {	
			if (u_mesh.dimat(u_mesh.d2v(uz_d2v[col])) & boundary_part) {
				lhs_mat.coeffRef(col, col) = 1.0;
				rhs_vec(col) = u;  // Dirichlet bc
			}
		}
		col++;
	}
}

// TODO: Apply bc implicitly by modifying the rhs.
// TODO: use boolean operators to select subdomains.
// This reduces the linear system size and preserves the symmetric form,
// which means we can utilize the Cholesky factorization (I think...)
// and other linear solvers which requires a symmetric linear system.
void pStokesProblem::apply_dirichlet_bc(
	const int boundary_part,
	const int velocity_component,
	std::function<long double(long double, long double)> u_func
) {
	// Loop over all nnz elements and set off-diagonals corresponding to boundary nodes to zero
	int *outer_end = lhs_mat.outerIndexPtr() + lhs_mat.rows();
	int col = 0;
	for (
		int *outer_p = lhs_mat.outerIndexPtr();
		outer_p != outer_end;
		outer_p++
	) {
		int nnz = *(outer_p + 1) - *outer_p;
		for (int i = 0; i < nnz; i++) {
			int row = *(lhs_mat.innerIndexPtr() + *outer_p + i);
			if (velocity_component == HORIZONTAL && ux_d2v(row) != -1) {
				if (u_mesh.dimat(u_mesh.d2v(ux_d2v(row))) & boundary_part) {
					long double *val_ptr = lhs_mat.valuePtr() + *outer_p + i;
					*val_ptr = 0.0;
				}
			}
			else if (velocity_component == VERTICAL && uz_d2v(row) != -1) {
				if (u_mesh.dimat(u_mesh.d2v(uz_d2v(row))) & boundary_part) {
					long double *val_ptr = lhs_mat.valuePtr() + *outer_p + i;
					*val_ptr = 0.0;
				}
			}
		}
		col++;
	}

	// Loop over all nnz elements and set diagonals corresponding to boundary nodes to one
	col = 0;
	for (
		int *outer_p = lhs_mat.outerIndexPtr();
		outer_p != outer_end;
		outer_p++
	) {
		if (velocity_component == HORIZONTAL && ux_d2v(col) != -1) {
			if (u_mesh.dimat(u_mesh.d2v(ux_d2v(col))) & boundary_part) {
				long double x = u_mesh.pmat(u_mesh.d2v(ux_d2v(col)), 0);
				long double z = u_mesh.pmat(u_mesh.d2v(ux_d2v(col)), 1);
				lhs_mat.coeffRef(col, col) = 1.0;  // FIXME: more efficient to use insert?
				rhs_vec(col) = u_func(x, z);  // Dirichlet bc
			}
		}
		else if (velocity_component == VERTICAL && uz_d2v(col) != -1) {	
			if (u_mesh.dimat(u_mesh.d2v(uz_d2v[col])) & boundary_part) {
				long double x = u_mesh.pmat(u_mesh.d2v(uz_d2v[col]), 0);
				long double z = u_mesh.pmat(u_mesh.d2v(uz_d2v[col]), 1);
				lhs_mat.coeffRef(col, col) = 1.0; // FIXME: more efficient to use insert?
				rhs_vec(col) = u_func(x, z);  // Dirichlet bc
			}
		}
		col++;
	}
}
// TODO: Apply bc implicitly by modifying the rhs.
// TODO: use boolean operators to select subdomains.
// This reduces the linear system size and preserves the symmetric form,
// which means we can utilize the Cholesky factorization (I think...)
// and other linear solvers which requires a symmetric linear system.
void pStokesProblem::apply_homogenous_dirichlet_bc(const Eigen::VectorXi &constrained_dofs)
{
	// Loop over all nnz elements and set off-diagonals corresponding to boundary nodes to zero
	int *outer_end = lhs_mat.outerIndexPtr() + lhs_mat.rows();
	int col = 0;
	for (
		int *outer_p = lhs_mat.outerIndexPtr();
		outer_p != outer_end;
		outer_p++
	) {
		int nnz = *(outer_p + 1) - *outer_p;
		if (std::binary_search(constrained_dofs.begin(), constrained_dofs.end(), col)) {
			for (int i = 0; i < nnz; i++) {
				long double *val_ptr = lhs_mat.valuePtr() + *outer_p + i;
				*val_ptr = 0.0;
			}
		} else {
			for (int i = 0; i < nnz; i++) {
				int row = *(lhs_mat.innerIndexPtr() + *outer_p + i);
				if (std::binary_search(constrained_dofs.begin(), constrained_dofs.end(), row)) {
					long double *val_ptr = lhs_mat.valuePtr() + *outer_p + i;
					*val_ptr = 0.0;
				}
			}
		}
		col++;
	}

	// Loop over all nnz elements and set diagonals corresponding to boundary nodes to one
	col = 0;
	for (
		int *outer_p = lhs_mat.outerIndexPtr();
		outer_p != outer_end;
		outer_p++
	) {
		if (std::binary_search(constrained_dofs.begin(), constrained_dofs.end(), col)) {
			lhs_mat.coeffRef(col, col) = 1.0;  // FIXME: more efficient to use insert?
			rhs_vec(col) = 0.0;  // Dirichlet bc
		}
		col++;
	}
}

void pStokesProblem::prune_lhs(long double ref, long double prec) {
	// NOTE: Pruning zeros is a costly operation and should only be done once.
	lhs_mat.prune(ref, prec);
}

void pStokesProblem::apply_noslip_bc(const int boundary_part)
{
	apply_dirichlet_bc(boundary_part, HORIZONTAL, 0.0);
	apply_dirichlet_bc(boundary_part, VERTICAL, 0.0);
}

void pStokesProblem::apply_impenetrability_bc(const int boundary_part)
{
	if (boundary_part & WEST_ID || boundary_part == EAST_ID)
		apply_dirichlet_bc(boundary_part, HORIZONTAL, 0.0);
	else if (boundary_part == BED_ID || boundary_part == SURFACE_ID)
		apply_dirichlet_bc(boundary_part, VERTICAL, 0.0);
}

void pStokesProblem::commit_lhs_mat()
{
	lhs_mat.setFromTriplets(lhs_coeffs.begin(), lhs_coeffs.end());
}

// TODO: Make solve step backend agnostic
void pStokesProblem::solve_linear(std::string sp_solver_name)
{
	if (sp_solver_name.compare("sparselu") == 0) {
		Eigen::SparseLU<Eigen::SparseMatrix<long double>> sp_solver;
		sp_solver.analyzePattern(lhs_mat);
		sp_solver.factorize(lhs_mat);
		w_vec = sp_solver.solve(rhs_vec);
	} else if (sp_solver_name.compare("umfpack") == 0) {
		Eigen::UmfPackLU<Eigen::SparseMatrix<long double>> sp_solver;
		sp_solver.analyzePattern(lhs_mat);
		sp_solver.factorize(lhs_mat);
		w_vec = sp_solver.solve(rhs_vec);
	}
}

// TODO: Make solve step backend agnostic
void pStokesProblem::solve_linear_iterative(std::string sp_solver_name)
{
	if (sp_solver_name.compare("BiCGSTAB") == 0) {
		Eigen::BiCGSTAB<Eigen::SparseMatrix<long double>> solver;
		solver.compute(lhs_mat);
		w_vec = solver.solve(rhs_vec);
	}
}

bool pStokesProblem::is_unstable(long double uz_tol_myr)
{
	long double uz_max = abs(w_vec(uz_v2d).maxCoeff());
	long double uz_min = abs(w_vec(uz_v2d).minCoeff());
	long double uz_max_norm = uz_max > uz_min ? uz_max : uz_min;
	return uz_max_norm > 1e-3*uz_tol_myr;
}

Eigen::SparseMatrix<long double> pStokesProblem::extract_A_block() 
{
	return lhs_mat.extract_block(u_v2d, u_v2d);
}

void pStokesProblem::reset_lhs()
{
	lhs_coeffs.clear();
	nof_pushed_elements = 0;
}

void pStokesProblem::reset_rhs()
{
	rhs_vec.setZero();
}

void pStokesProblem::reset_system()
{
	reset_lhs();
	reset_rhs();
}

FEMFunction2D pStokesProblem::velocity_x()
{
	FEMFunction2D ux_fem_func(u_mesh);
	ux_fem_func.vals = w_vec(ux_v2d);
	return ux_fem_func;
}

FEMFunction2D pStokesProblem::velocity_z()
{
	FEMFunction2D uz_fem_func(u_mesh);
	uz_fem_func.vals = w_vec(uz_v2d);
	return uz_fem_func;
}

FEMFunction2D pStokesProblem::pressure()
{
	FEMFunction2D p_fem_func(p_mesh);
	p_fem_func.vals = w_vec(p_v2d);
	return p_fem_func;
}
