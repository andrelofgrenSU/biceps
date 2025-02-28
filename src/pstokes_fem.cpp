#include <boost/format.hpp>
#include <pstokes_fem.hpp>
#include <fem_1d.hpp>
#include <fem_2d.hpp>

static inline FloatType glen_flow_viscosity(
	FloatType A, FloatType n_i, FloatType eps_reg_2, FloatType eff_strain_rate_2
) {
	return 0.5*POW_FUNC(A, -1.0/n_i)*POW_FUNC(
		eff_strain_rate_2 + eps_reg_2, (1.0-n_i)/(2.0*n_i)
	);
}

pStokesProblem::pStokesProblem(
	FloatType rate_factor,
	FloatType glen_exp,
	FloatType eps_reg_2,
	std::function<FloatType(FloatType, FloatType)> force_x,
	std::function<FloatType(FloatType, FloatType)> force_z,
	StructuredMesh &u_mesh,
	StructuredMesh &p_mesh
) : 
	A(rate_factor), n_i(glen_exp), eps_reg_2(eps_reg_2),
	force_x(force_x), force_z(force_z),
	u_mesh(u_mesh), p_mesh(p_mesh), logger(WARN, std::cout)
{
	nv_dofs = u_mesh.nof_dofs();
	np_dofs = p_mesh.nof_dofs();
	n_dofs = 2*nv_dofs + np_dofs; 

	rhs_vec = Eigen::VectorX<FloatType>::Zero(n_dofs);

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

	w_vec = Eigen::VectorX<FloatType>::Zero(n_dofs);
	lhs_mat = Eigen::SparseMatrix<FloatType>(n_dofs, n_dofs);

	// Estimate max nnz
	int u_nnz_per_dof = 2*(2*u_mesh.degree() + 1)*(2*u_mesh.degree() + 1);
	int p_nnz_per_dof = (2*p_mesh.degree() + 1)*(2*p_mesh.degree() + 1);
	int nnz_per_dof = u_nnz_per_dof + p_nnz_per_dof;
	lhs_coeffs.reserve(nnz_per_dof*n_dofs);
}

// FIXME: Is there a more efficient way of inserting elements directly into the matrix?
void pStokesProblem::assemble_stress_block()
{
	Eigen::MatrixX<FloatType> node_coords_u, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> qweights, detJ_rs;
	Eigen::VectorXi element_u;
	
	FEM2D::gauss_legendre_quadrature(
		gp_stress, u_mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		u_mesh.degree(),
		u_mesh.cell_type(),
		qpoints_rs,
		phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), u_mesh.dofs_per_cell());

	Eigen::MatrixX<FloatType> A_xx = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_cell(), u_mesh.dofs_per_cell()
	);
	Eigen::MatrixX<FloatType> A_xz = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_cell(), u_mesh.dofs_per_cell()
	);
	Eigen::MatrixX<FloatType> A_zz = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_cell(), u_mesh.dofs_per_cell()
	);

	Eigen::MatrixX<FloatType> A_block = Eigen::MatrixX<FloatType>::Zero(
		2*u_mesh.dofs_per_cell(), 2*u_mesh.dofs_per_cell()
	);

	for (int k = 0; k < u_mesh.cmat.rows(); k++) {
		element_u = u_mesh.cmat.row(k);
		node_coords_u = u_mesh.pmat(element_u, Eigen::all);
		FEM2D::map_to_reference_cell(
			u_mesh.degree(), u_mesh.cell_type(), node_coords_u, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		Eigen::VectorX<FloatType> ux_vec = w_vec(ux_v2d(element_u));
		Eigen::VectorX<FloatType> uz_vec = w_vec(uz_v2d(element_u));
		Eigen::MatrixX<FloatType> Sx = dphi_xz(Eigen::seq(0, Eigen::last, 2), Eigen::all);
		Eigen::MatrixX<FloatType> Sz = dphi_xz(Eigen::seq(1, Eigen::last, 2), Eigen::all);
		Eigen::VectorX<FloatType> ddx_ux = Sx*ux_vec;
		Eigen::VectorX<FloatType> ddz_ux = Sz*ux_vec;
		Eigen::VectorX<FloatType> ddx_uz = Sx*uz_vec;
		Eigen::VectorX<FloatType> ddz_uz = Sz*uz_vec;

		for (int q = 0; q < qpoints_rs.rows(); q++) {
			int w = 2*q;
			FloatType eff_strain_rate_2 = 0.5*(
				POW_FUNC(ddx_ux(q), 2)
				+ 0.5*POW_FUNC(ddz_ux(q) + ddx_uz(q), 2)
				+ POW_FUNC(ddz_uz(q), 2)
			);
			// Glen's flow law rheology
			FloatType eta = glen_flow_viscosity(
				A, n_i, eps_reg_2, eff_strain_rate_2
			);
			FloatType detJxWxEta = ABS_FUNC(detJ_rs(q))*qweights(q)*eta;
			// TODO: tabulate dphi*dphi and calculate in advance
			for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < u_mesh.dofs_per_cell(); j++) {
					A_xx(i, j) += (2.0*dphi_xz(w, i)*dphi_xz(w, j) + dphi_xz(w+1, i)*dphi_xz(w+1, j))*detJxWxEta;
					A_xz(i, j) += dphi_xz(w+1, i)*dphi_xz(w, j)*detJxWxEta;
					A_zz(i, j) += (dphi_xz(w, i)*dphi_xz(w, j) + 2.0*dphi_xz(w+1, i)*dphi_xz(w+1, j))*detJxWxEta;
				}
			}
		}

		for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < u_mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					ux_v2d(element_u(i)),
					ux_v2d(element_u(j)),
					A_xx(i, j)
				));
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					ux_v2d(element_u(i)),
					uz_v2d(element_u(j)),
					A_xz(i, j)
				));
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					uz_v2d(element_u(i)),
					ux_v2d(element_u(j)),
					A_xz(j, i)
				));
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					uz_v2d(element_u(i)),
					uz_v2d(element_u(j)),
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

void pStokesProblem::assemble_incomp_block()
{
	Eigen::MatrixX<FloatType> node_coords_u, node_coords_p, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz, psi_rs, dpsi_rs, dpsi_xz;
	Eigen::VectorX<FloatType> qweights, detJ_rs;
	Eigen::VectorXi element_u, element_p;
	
	FEM2D::gauss_legendre_quadrature(
		gp_incomp, u_mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		p_mesh.degree(), p_mesh.cell_type(), qpoints_rs, psi_rs, dpsi_rs
	);

	FEM2D::lagrange_basis(
		u_mesh.degree(), u_mesh.cell_type(), qpoints_rs, phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), u_mesh.dofs_per_cell());
	dpsi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), p_mesh.dofs_per_cell());

	Eigen::MatrixX<FloatType> B_px = Eigen::MatrixX<FloatType>::Zero(
		p_mesh.dofs_per_cell(), u_mesh.dofs_per_cell()
	);
	Eigen::MatrixX<FloatType> B_pz = Eigen::MatrixX<FloatType>::Zero(
		p_mesh.dofs_per_cell(), u_mesh.dofs_per_cell()
	);
	Eigen::MatrixX<FloatType> B_xp = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_cell(), p_mesh.dofs_per_cell()
	);
	Eigen::MatrixX<FloatType> B_zp = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_cell(), p_mesh.dofs_per_cell()
	);

	for (int k = 0; k < u_mesh.cmat.rows(); k++) {
		element_u = u_mesh.cmat.row(k);
		element_p = p_mesh.cmat.row(k);
		node_coords_u = u_mesh.pmat(element_u, Eigen::all);
		node_coords_p = p_mesh.pmat(element_p, Eigen::all);

		// TODO: only psi_rs is need; speed up by discard calculating everything else
		FEM2D::map_to_reference_cell(
			p_mesh.degree(), p_mesh.cell_type(), node_coords_p, qpoints_rs,
			psi_rs, dpsi_rs, detJ_rs, qpoints_xz, dpsi_xz
		);

		FEM2D::map_to_reference_cell(
			u_mesh.degree(), u_mesh.cell_type(), node_coords_u, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		for (int q = 0; q < qpoints_rs.rows(); q++) {
			int w = 2*q;
			FloatType detJxW = ABS_FUNC(detJ_rs(q))*qweights(q);
			for (int i = 0; i < p_mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < u_mesh.dofs_per_cell(); j++) {
					B_px(i, j) += -psi_rs(q, i)*dphi_xz(w, j)*detJxW;
					B_pz(i, j) += -psi_rs(q, i)*dphi_xz(w+1, j)*detJxW;
				}
			}
			for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < p_mesh.dofs_per_cell(); j++) {
					B_xp(i, j) += -dphi_xz(w, i)*psi_rs(q, j)*detJxW;
					B_zp(i, j) += -dphi_xz(w+1, i)*psi_rs(q, j)*detJxW;
				}
			}
		}

		for (int i = 0; i < p_mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < u_mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					p_v2d(element_p(i)),
					ux_v2d(element_u(j)),
					B_px(i, j)
				));
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					p_v2d(element_p(i)),
					uz_v2d(element_u(j)),
					B_pz(i, j)
				));
			}
		}
		for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < p_mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					ux_v2d(element_u(i)),
					p_v2d(element_p(j)),
					B_xp(i, j)
				));
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					uz_v2d(element_u(i)),
					p_v2d(element_p(j)),
					B_zp(i, j)
				));
			}
		}

		B_px.setZero();
		B_pz.setZero();
		B_xp.setZero();
		B_zp.setZero();
	}
	nof_pushed_elements = lhs_coeffs.size();
}

// TODO: add fx component to the weak form
void pStokesProblem::assemble_fssa_vertical_block()
{
	Eigen::MatrixX<FloatType> node_coords, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi edge_vi;

	FEM1D::gauss_legendre_quadrature(
		gp_fssa, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		u_mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<FloatType>::Zero(
		2*qpoints_r.rows(), u_mesh.dofs_per_edge()
	);

	Eigen::MatrixX<FloatType> A_xz = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);
	Eigen::MatrixX<FloatType> A_zz = Eigen::MatrixX<FloatType>::Zero(
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

		FloatType nx = u_mesh.edge_normals(si, 0);
		FloatType nz = u_mesh.edge_normals(si, 1);
		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType fz = force_z(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType dFxW = ABS_FUNC(detJ_r(q))*qweights(q);
			for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
				for (int j = 0; j < u_mesh.dofs_per_edge(); j++) {
					A_xz(i, j) -= fssa_param*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nx;
					A_zz(i, j) -= fssa_param*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nz;
				}
			}
		}

		for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
			for (int j = 0; j < u_mesh.dofs_per_edge(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					uz_v2d(edge_vi(i)),
					ux_v2d(edge_vi(j)),
					A_xz(i, j)
				));
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
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

// TODO: add fx component to the weak form
void pStokesProblem::assemble_fssa_normal_block()
{
	Eigen::MatrixX<FloatType> node_coords_u, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi edge_vi;

	FEM1D::gauss_legendre_quadrature(
		gp_fssa, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		u_mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<FloatType>::Zero(
		2*qpoints_r.rows(), u_mesh.dofs_per_edge()
	);

	Eigen::MatrixX<FloatType> A_xx = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);
	Eigen::MatrixX<FloatType> A_xz = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);
	Eigen::MatrixX<FloatType> A_zx = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);
	Eigen::MatrixX<FloatType> A_zz = Eigen::MatrixX<FloatType>::Zero(
		u_mesh.dofs_per_edge(), u_mesh.dofs_per_edge()
	);

	std::vector<int> surf_edge_inds = u_mesh.extract_edge_inds(SURFACE_ID);
	for (int si: surf_edge_inds) {
		edge_vi = u_mesh.emat(si, Eigen::all);
		node_coords_u = u_mesh.pmat(edge_vi, Eigen::all);
		FEM1D::map_to_reference_cell(
			u_mesh.degree(), node_coords_u, qpoints_r,
			phi_r, dphi_r, detJ_r, qpoints_x, dphi_x
		);
		FloatType nz = sqrt(u_mesh.edge_normals(si, 1));
		FloatType nx = u_mesh.edge_normals(si, 0)/nz;


		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType fz = force_z(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType dFxW = ABS_FUNC(detJ_r(q))*qweights(q);
			for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
				for (int j = 0; j < u_mesh.dofs_per_edge(); j++) {
					A_xx(i, j) -= fssa_param*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nx*nx;
					A_xz(i, j) -= fssa_param*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nx*nz;
					A_zx(i, j) -= fssa_param*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nz*nx;
					A_zz(i, j) -= fssa_param*phi_r(q, i)*phi_r(q, j)*fz*dFxW*nz*nz;
				}
			}
		}

		for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
			for (int j = 0; j < u_mesh.dofs_per_edge(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					ux_v2d(edge_vi(i)),
					ux_v2d(edge_vi(j)),
					A_xx(i, j)
				));
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					ux_v2d(edge_vi(i)),
					uz_v2d(edge_vi(j)),
					A_xz(i, j)
				));
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					uz_v2d(edge_vi(i)),
					ux_v2d(edge_vi(j)),
					A_zx(i, j)
				));
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					uz_v2d(edge_vi(i)),
					uz_v2d(edge_vi(j)),
					A_zz(i, j)
				));
			}
		}

		A_xx.setZero();
		A_xz.setZero();
		A_zx.setZero();
		A_zz.setZero();
	}
	nof_pushed_elements = lhs_coeffs.size();
}

void pStokesProblem::assemble_fssa_vertical_rhs() {
	Eigen::MatrixX<FloatType> node_coords_u, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi edge_vi;

	FEM1D::gauss_legendre_quadrature(
		gp_rhs, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		u_mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<FloatType>::Zero(
		2*qpoints_r.rows(), u_mesh.dofs_per_edge()
	);

	std::vector<int> surf_edge_inds = u_mesh.extract_edge_inds(SURFACE_ID);
	for (int si: surf_edge_inds) {
		edge_vi = u_mesh.emat(si, Eigen::all);
		node_coords_u = u_mesh.pmat(edge_vi, Eigen::all);
		FEM1D::map_to_reference_cell(
			u_mesh.degree(), node_coords_u, qpoints_r,
			phi_r, dphi_r, detJ_r, qpoints_x, dphi_x
		);
		FloatType nz = u_mesh.edge_normals(si, 1);
		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType fx = force_x(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType fz = force_z(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType ac = fssa_accum(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType dFxW = ABS_FUNC(detJ_r(q))*qweights(q);
			for (int i = 0; i < u_mesh.dofs_per_edge(); i++) {
				rhs_vec(ux_v2d(edge_vi(i))) += fssa_param*phi_r(q, i)*ac*fx*dFxW*nz;
				rhs_vec(uz_v2d(edge_vi(i))) += fssa_param*phi_r(q, i)*ac*fz*dFxW*nz;
			}
		}
	}
}

void pStokesProblem::assemble_rhs_vec() {
	Eigen::MatrixX<FloatType> node_coords_u, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> qweights, detJ_rs;
	Eigen::VectorXi element_u;

	FEM2D::gauss_legendre_quadrature(
		gp_rhs, u_mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		u_mesh.degree(), u_mesh.cell_type(), qpoints_rs, phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), u_mesh.dofs_per_cell());

	for (int k = 0; k < u_mesh.cmat.rows(); k++) {
		element_u = u_mesh.cmat.row(k);
		node_coords_u = u_mesh.pmat(element_u, Eigen::all);
		FEM2D::map_to_reference_cell(
			u_mesh.degree(), u_mesh.cell_type(), node_coords_u, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		for (int q = 0; q < qpoints_rs.rows(); q++) {
			FloatType detJxW = ABS_FUNC(detJ_rs(q))*qweights(q);
			for (int i = 0; i < u_mesh.dofs_per_cell(); i++) {
				rhs_vec(ux_v2d(element_u(i))) += (
					force_x(qpoints_xz(q, 0), qpoints_xz(q, 1))*phi_rs(q, i)
				)*detJxW;
				rhs_vec(uz_v2d(element_u(i))) += (
					force_z(qpoints_xz(q, 0), qpoints_xz(q, 1))*phi_rs(q, i)
				)*detJxW;
			}
		}
	}
}

void pStokesProblem::commit_lhs_mat()
{
	lhs_mat.setFromTriplets(lhs_coeffs.begin(), lhs_coeffs.end());
}

void pStokesProblem::apply_standard_dirichlet_bc()
{
	// Constrained dofs for standard bc case (no slip on bed and impenetrability on sides)
	int boundary_ux = NORTH_WEST_ID | WEST_ID | BED_ID | EAST_ID | NORTH_EAST_ID;
	int boundary_uz = BED_ID;

	// Loop over all nnz elements and set off-diagonals corresponding to boundary nodes to zero
	for (int col = 0; col < lhs_mat.outerSize(); ++col)
		if (
			ux_d2v(col) != -1 && u_mesh.dimat(u_mesh.d2v(ux_d2v(col))) & boundary_ux ||
			uz_d2v(col) != -1 && u_mesh.dimat(u_mesh.d2v(uz_d2v(col))) & boundary_uz
		) {
			// Eliminate column
			for (Eigen::SparseMatrix<FloatType>::InnerIterator it(lhs_mat, col); it; ++it) {
				it.valueRef() = 0.0;
			}
			lhs_mat.coeffRef(col, col) = 1.0;
			rhs_vec(col) = 0.0; // Dirichlet bc
		} else {
			for (Eigen::SparseMatrix<FloatType>::InnerIterator it(lhs_mat, col); it; ++it) {
				int row = it.row();
				if (
					ux_d2v(row) != -1 && u_mesh.dimat(u_mesh.d2v(ux_d2v(row))) & boundary_ux ||
					uz_d2v(row) != -1 && u_mesh.dimat(u_mesh.d2v(uz_d2v(row))) & boundary_uz
				) {
					// Eliminate row
					it.valueRef() = 0.0;
				}
			}
		}
	prune_lhs();
}

void pStokesProblem::apply_dirichlet_bc(
	const int boundary_part,
	const int velocity_component,
	std::function<FloatType(FloatType, FloatType)> u_func
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
					FloatType *val_ptr = lhs_mat.valuePtr() + *outer_p + i;
					*val_ptr = 0.0;
				}
			}
			else if (velocity_component == VERTICAL && uz_d2v(row) != -1) {
				if (u_mesh.dimat(u_mesh.d2v(uz_d2v(row))) & boundary_part) {
					FloatType *val_ptr = lhs_mat.valuePtr() + *outer_p + i;
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
				FloatType x = u_mesh.pmat(u_mesh.d2v(ux_d2v(col)), 0);
				FloatType z = u_mesh.pmat(u_mesh.d2v(ux_d2v(col)), 1);
				lhs_mat.coeffRef(col, col) = 1.0;  // FIXME: more efficient to use insert?
				rhs_vec(col) = u_func(x, z);  // Dirichlet bc
			}
		}
		else if (velocity_component == VERTICAL && uz_d2v(col) != -1) {	
			if (u_mesh.dimat(u_mesh.d2v(uz_d2v[col])) & boundary_part) {
				FloatType x = u_mesh.pmat(u_mesh.d2v(uz_d2v(col)), 0);
				FloatType z = u_mesh.pmat(u_mesh.d2v(uz_d2v(col)), 1);
				lhs_mat.coeffRef(col, col) = 1.0; // FIXME: more efficient to use insert?
				rhs_vec(col) = u_func(x, z);  // Dirichlet bc
			}
		}
		col++;
	}
}

void pStokesProblem::prune_lhs() {
	// NOTE: Pruning zeros is a costly operation and should only be done once.
	lhs_mat.prune(prune_threshold, 1);
}

void pStokesProblem::solve_linear_system()
{
	if (linear_solver == SPARSE_LU) {
		Eigen::SparseLU<Eigen::SparseMatrix<FloatType>> solver;
		solver.analyzePattern(lhs_mat);
		solver.factorize(lhs_mat);
		w_vec = solver.solve(rhs_vec);
	} else if (linear_solver == BICGSTAB) {
		Eigen::BiCGSTAB<Eigen::SparseMatrix<FloatType>> solver;
		solver.compute(lhs_mat);
		w_vec = solver.solve(rhs_vec);
	} else {
		// TODO: throw exception
	}
}

void pStokesProblem::solve_picard(bool reset_linear_system)
{
	Eigen::SparseMatrix<FloatType> M = FEM2D::assemble_mass_matrix(u_mesh, u_mesh.degree() + 1);
	assemble_incomp_block();
	if (fssa_version == FSSA_VERTICAL)
		assemble_fssa_vertical_block();
	else if (fssa_version == FSSA_NORMAL)
		assemble_fssa_normal_block();
	assemble_rhs_vec();
	int nof_pushed_elements_before_stress_assembly = nof_pushed_elements;

	FloatType error_abs;
	FloatType area_sqrt = sqrt(FEM2D::calculate_area(u_mesh, u_mesh.degree() + 1));

	FEMFunction2D ux_old = velocity_x();
	FEMFunction2D uz_old = velocity_z();

	// Picard loop
	int ip = 0;
	do {
		lhs_coeffs.resize(nof_pushed_elements_before_stress_assembly);
		nof_pushed_elements = nof_pushed_elements_before_stress_assembly;
		assemble_stress_block();
		commit_lhs_mat();

		apply_standard_dirichlet_bc();
		solve_linear_system();

		// Estimate the error
		FEMFunction2D ux = velocity_x();
		FEMFunction2D uz = velocity_z();
		FEMFunction2D ex_func = ux - ux_old;
		FEMFunction2D ez_func = uz - uz_old;
		FloatType ex = FEM2D::L2_norm(ex_func, M);
		FloatType ez = FEM2D::L2_norm(ez_func, M);
		error_abs = sqrt(ex*ex + ez*ez)/area_sqrt;
		logger.log_msg((boost::format{"<pStokes> Picard iteration %d: e(abs) = %g"} %(ip+1) %error_abs).str(), TRACE);

		ux_old.vals = ux.vals;
		uz_old.vals = uz.vals;
		ip++;
	} while (error_abs > picard_atol && ip < picard_max_iter);

	logger.log_msg("-------------------------------------", INFO);
	logger.log_msg("********** pStokes Summary **********", INFO);
	logger.log_msg("-------------------------------------", INFO);

	logger.log_msg((boost::format{"Picard solver finished in %d iterations"} %ip).str(), INFO);
	logger.log_msg((boost::format{"e (abs): %.4g"} % (error_abs)).str(), INFO);
	logger.log_msg((boost::format{"min ux, max ux: %.4g, %.4g m/yr"} % (1e3*w_vec(ux_v2d).minCoeff()) % (1e3*w_vec(ux_v2d).maxCoeff())).str(), INFO);
	logger.log_msg((boost::format{"min uz, max uz: %.4g, %.4g m/yr"} % (1e3*w_vec(uz_v2d).minCoeff()) % (1e3*w_vec(uz_v2d).maxCoeff())).str(), INFO);
	logger.log_msg((boost::format{"min p, max p: %.4g, %.4g MPa"} % (w_vec(p_v2d).minCoeff()) % (w_vec(p_v2d).maxCoeff())).str(), INFO);
	logger.log_msg("-------------------------------------", INFO);

	if (reset_linear_system)
		reset_system();
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
