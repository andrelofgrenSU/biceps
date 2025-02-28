#include <poisson_fem.hpp>
#include <fem_1d.hpp>
#include <fem_2d.hpp>

PoissonProblem::PoissonProblem(int nx, int nz, int deg, int cell_type) : 
	mesh(nx, nz, deg, cell_type)
{
	n_dofs = mesh.nof_dofs();
	lhs_mat = Eigen::SparseMatrix<FloatType>(n_dofs, n_dofs);
	rhs_vec = Eigen::VectorX<FloatType>::Zero(n_dofs);
	sol_vec = Eigen::VectorX<FloatType>::Zero(n_dofs);
	lhs_coeffs.reserve(6*mesh.lhs_nnz());
}

void PoissonProblem::reset()
{
	lhs_coeffs.clear();
	lhs_coeffs.reserve(6*mesh.lhs_nnz());
	rhs_vec.setZero();
	sol_vec.setZero();
}

void PoissonProblem::assemble_stiffness_block(
	std::function<FloatType(FloatType, FloatType)> alpha,
	std::function<FloatType(FloatType, FloatType)> beta,
	int gp
) {
	Eigen::MatrixX<FloatType> node_coords_v, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> qweights, detJ_rs;
	Eigen::VectorXi element;
	
	FEM2D::gauss_legendre_quadrature(
		gp, mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		mesh.degree(),
		mesh.cell_type(),
		qpoints_rs,
		phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), mesh.dofs_per_cell());

	Eigen::MatrixX<FloatType> A_block = Eigen::MatrixX<FloatType>::Zero(
		mesh.dofs_per_cell(), mesh.dofs_per_cell()
	);

	for (int k = 0; k < mesh.cmat.rows(); k++) {
		element = mesh.cmat.row(k);
		node_coords_v = mesh.pmat(element, Eigen::all);
		FEM2D::map_to_reference_cell(
			mesh.degree(), mesh.cell_type(), node_coords_v, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		for (int q = 0; q < qpoints_rs.rows(); q++) {
			int w = 2*q;
			FloatType alpha_xz = alpha(qpoints_xz(q, 0), qpoints_xz(q, 1));
			FloatType beta_xz = beta(qpoints_xz(q, 0), qpoints_xz(q, 1));
			FloatType detJxW = ABS_FUNC(detJ_rs(q))*qweights(q);
			for (int i = 0; i < mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < mesh.dofs_per_cell(); j++) {
					A_block(i, j) += (
						alpha_xz*dphi_xz(w, i)*dphi_xz(w, j) + beta_xz*dphi_xz(w+1, i)*dphi_xz(w+1, j)
					)*detJxW;
				}
			}
		}

		for (int i = 0; i < mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					element(i),
					element(j),
					A_block(i, j)
				));
			}
		}	
		A_block.setZero();
	}
}

void PoissonProblem::assemble_mass_block(
	std::function<FloatType(FloatType, FloatType)> gamma, int gp
) {
	Eigen::MatrixX<FloatType> node_coords_v, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> qweights, detJ_rs;
	Eigen::VectorXi element;
	
	FEM2D::gauss_legendre_quadrature(
		gp, mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		mesh.degree(),
		mesh.cell_type(),
		qpoints_rs,
		phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), mesh.dofs_per_cell());

	Eigen::MatrixX<FloatType> M_block = Eigen::MatrixX<FloatType>::Zero(
		mesh.dofs_per_cell(), mesh.dofs_per_cell()
	);

	for (int k = 0; k < mesh.cmat.rows(); k++) {
		element = mesh.cmat.row(k);
		node_coords_v = mesh.pmat(element, Eigen::all);
		FEM2D::map_to_reference_cell(
			mesh.degree(), mesh.cell_type(), node_coords_v, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		for (int q = 0; q < qpoints_rs.rows(); q++) {
			int w = 2*q;

			FloatType detJxWxGamma = gamma(
				qpoints_xz(q, 0), qpoints_xz(q, 1)
			)*ABS_FUNC(detJ_rs(q))*qweights(q);
			for (int i = 0; i < mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < mesh.dofs_per_cell(); j++) {
					M_block(i, j) += phi_rs(q, i)*phi_rs(q, j)*detJxWxGamma;
				}
			}
		}

		for (int i = 0; i < mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					element(i),
					element(j),
					M_block(i, j)
				));
			}
		}	
		M_block.setZero();
	}
}

void PoissonProblem::assemble_FSSA_block(
	FEMFunction2D &f,
	FloatType theta,
	FloatType dt,
	int gp
) {
	Eigen::MatrixX<FloatType> node_coords, qpoints_xz, qpoints_rs, phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> qweights, detJ_rs, f_cell;
	Eigen::VectorXi cell_vi;
	
	FEM2D::gauss_legendre_quadrature_edge_quadrilateral(
		gp, NORTH_ID, qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		mesh.degree(), mesh.cell_type(), qpoints_rs, phi_rs, dphi_rs
	);

	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(
		2*qpoints_rs.rows(), mesh.dofs_per_cell()
	);

	Eigen::MatrixX<FloatType> A_block = Eigen::MatrixX<FloatType>::Zero(
		mesh.dofs_per_cell(), mesh.dofs_per_cell()
	);

	std::vector<int> surf_cell_inds = mesh.extract_cell_inds(SURFACE_ID);
	for (int si: surf_cell_inds) {
		cell_vi = mesh.cmat(si, Eigen::all);
		node_coords = mesh.pmat(cell_vi, Eigen::all);
		FEM2D::map_to_reference_cell(
			mesh.degree(), mesh.cell_type(), node_coords, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		Eigen::VectorX<FloatType> x0 = mesh.pmat(cell_vi(3*mesh.degree()), Eigen::all);
		Eigen::VectorX<FloatType> x1 = mesh.pmat(cell_vi(2*mesh.degree()), Eigen::all);
		FloatType dFdr = (x1-x0).norm();
		f_cell = f.eval_cell(si);
		for (int q = 0; q < qpoints_rs.rows(); q++) {
			int w = 2*q;
			FloatType dFdrxWxfxThetaxDT =
				f_cell.dot(phi_rs(q, Eigen::all))*qweights(q)*dFdr*dt*theta;
			for (int i = 0; i < mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < mesh.dofs_per_cell(); j++) {
					A_block(i, j) += 
						phi_rs(q, i)*dphi_xz(w, j)*dFdrxWxfxThetaxDT;
				}
			}
		}

		for (int i = 0; i < mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					cell_vi(i),
					cell_vi(j),
					A_block(i, j)
				));
			}
		}	
		A_block.setZero();
	}
}

void PoissonProblem::assemble_robin_block(
	std::function<FloatType(FloatType, FloatType)> g_robin,
	std::function<FloatType(FloatType, FloatType)> a_robin,
	std::function<FloatType(FloatType, FloatType)> b_robin,
	int boundary_id,
	int gp
) {
	Eigen::MatrixX<FloatType> node_coords, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi edge_vi;
	
	FEM1D::gauss_legendre_quadrature(
		gp, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<FloatType>::Zero(
		2*qpoints_r.rows(), mesh.dofs_per_cell()
	);

	Eigen::MatrixX<FloatType> M_block = Eigen::MatrixX<FloatType>::Zero(
		mesh.dofs_per_edge(), mesh.dofs_per_edge()
	);

	std::vector<int> surf_edge_inds = mesh.extract_edge_inds(boundary_id);
	for (int si: surf_edge_inds) {
		edge_vi = mesh.emat(si, Eigen::all);
		node_coords = mesh.pmat(edge_vi, Eigen::all);
		FEM1D::map_to_reference_cell(
			mesh.degree(), node_coords, qpoints_r,
			phi_r, dphi_r, detJ_r, qpoints_x, dphi_x
		);

		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType ar = a_robin(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType br = b_robin(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType gr = g_robin(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType cr = ar/br;
			FloatType dFxWxCr = cr*ABS_FUNC(detJ_r(q))*qweights(q);
			// TODO: optimize by precalculating thetaxdtxfzx...
			for (int i = 0; i < mesh.dofs_per_edge(); i++) {
				for (int j = 0; j < mesh.dofs_per_edge(); j++) {
					M_block(i, j) += phi_r(q, i)*phi_r(q, j)*dFxWxCr;
				}
			}
		}

		for (int i = 0; i < mesh.dofs_per_edge(); i++) {
			for (int j = 0; j < mesh.dofs_per_edge(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					edge_vi(i),
					edge_vi(j),
					M_block(i, j)
				));
			}
		}	
		M_block.setZero();
	}
}

void PoissonProblem::commit_lhs_mat()
{
	lhs_mat.setFromTriplets(lhs_coeffs.begin(), lhs_coeffs.end());
}

void PoissonProblem::assemble_neumann_rhs(
	std::function<FloatType(FloatType, FloatType)> g_neumann,
	int boundary_id,
	int gp
) {
	Eigen::MatrixX<FloatType> node_coords, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi edge_vi;
	
	FEM1D::gauss_legendre_quadrature(
		gp, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<FloatType>::Zero(
		2*qpoints_r.rows(), mesh.dofs_per_cell()
	);

	std::vector<int> surf_edge_inds = mesh.extract_edge_inds(boundary_id);
	for (int si: surf_edge_inds) {
		edge_vi = mesh.emat(si, Eigen::all);
		node_coords = mesh.pmat(edge_vi, Eigen::all);
		FEM1D::map_to_reference_cell(
			mesh.degree(), node_coords, qpoints_r,
			phi_r, dphi_r, detJ_r, qpoints_x, dphi_x
		);

		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType gn = g_neumann(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType detJxWxGn = gn*ABS_FUNC(detJ_r(q))*qweights(q);
			for (int i = 0; i < mesh.dofs_per_cell(); i++) {
				 rhs_vec(edge_vi(i)) += phi_r(q, i)*detJxWxGn;
			}
		}
	}	
}

void PoissonProblem::assemble_force_rhs(
	std::function<FloatType(FloatType, FloatType)> f, int gp
) {
	Eigen::MatrixX<FloatType> node_coords_v, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> qweights, detJ_rs;
	Eigen::VectorXi element;
	
	FEM2D::gauss_legendre_quadrature(
		gp, mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		mesh.degree(), mesh.cell_type(), qpoints_rs, phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), mesh.dofs_per_cell());

	for (int k = 0; k < mesh.cmat.rows(); k++) {
		element = mesh.cmat.row(k);
		node_coords_v = mesh.pmat(element, Eigen::all);
		// FIXME: only phi_rs is need; speed up by discard calculating everything else
		FEM2D::map_to_reference_cell(
			mesh.degree(), mesh.cell_type(), node_coords_v, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);

		for (int q = 0; q < qpoints_rs.rows(); q++) {
			FloatType detJxWxf = f(
				qpoints_xz(q, 0), qpoints_xz(q, 1)
			)*ABS_FUNC(detJ_rs(q))*qweights(q);
			for (int i = 0; i < mesh.dofs_per_cell(); i++) {
				 rhs_vec(element(i)) += phi_rs(q, i)*detJxWxf;
			}
		}
	}	
}

void PoissonProblem::assemble_force_rhs(
	FEMFunction2D &f, int gp
) {
	Eigen::MatrixX<FloatType> node_coords_v, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> qweights, detJ_rs, f_cell;
	Eigen::VectorXi element;
	
	FEM2D::gauss_legendre_quadrature(
		gp, mesh.cell_type(), qpoints_rs, qweights
	);

	FEM2D::lagrange_basis(
		mesh.degree(), mesh.cell_type(), qpoints_rs, phi_rs, dphi_rs
	);

	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), mesh.dofs_per_cell());

	for (int k = 0; k < mesh.cmat.rows(); k++) {
		element = mesh.cmat.row(k);
		node_coords_v = mesh.pmat(element, Eigen::all);
		// FIXME: only phi_rs is need; speed up by discard calculating everything else
		FEM2D::map_to_reference_cell(
			mesh.degree(), mesh.cell_type(), node_coords_v, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);
		
		f_cell = f.eval_cell(k);
		for (int q = 0; q < qpoints_rs.rows(); q++) {
			FloatType detJxWxf = 
				f_cell.dot(phi_rs(q, Eigen::all))*ABS_FUNC(detJ_rs(q))*qweights(q);
			for (int i = 0; i < mesh.dofs_per_cell(); i++) {
				 rhs_vec(element(i)) += phi_rs(q, i)*detJxWxf;
			}
		}
	}	
}

void PoissonProblem::assemble_robin_rhs(
	std::function<FloatType(FloatType, FloatType)> b_robin,
	std::function<FloatType(FloatType, FloatType)> g_robin,
	int boundary_id,
	int gp
) {
	Eigen::MatrixX<FloatType> node_coords, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi edge_vi;
	
	FEM1D::gauss_legendre_quadrature(
		gp, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<FloatType>::Zero(
		2*qpoints_r.rows(), mesh.dofs_per_cell()
	);

	std::vector<int> surf_edge_inds = mesh.extract_edge_inds(boundary_id);
	for (int si: surf_edge_inds) {
		edge_vi = mesh.emat(si, Eigen::all);
		node_coords = mesh.pmat(edge_vi, Eigen::all);
		FEM1D::map_to_reference_cell(
			mesh.degree(), node_coords, qpoints_r,
			phi_r, dphi_r, detJ_r, qpoints_x, dphi_x
		);

		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType br = b_robin(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType gr = g_robin(qpoints_x(q, 0), qpoints_x(q, 1));
			FloatType cr = gr/br;
			FloatType detJxWxCr = cr*ABS_FUNC(detJ_r(q))*qweights(q);
			for (int i = 0; i < mesh.dofs_per_cell(); i++) {
				 rhs_vec(edge_vi(i)) += phi_r(q, i)*detJxWxCr;
			}
		}
	}	
}

void PoissonProblem::apply_zero_dirichlet_bc(int boundary_part) {
	std::vector<int> interior_dofs = mesh.extract_dof_inds(INTERIOR_ID);
	std::vector<int> boundary_dofs = mesh.extract_dof_inds(BOUNDARY_ID);
	Eigen::Map<Eigen::VectorXi> free_dofs(interior_dofs.data(), interior_dofs.size());
	lhs_mat_free = lhs_mat.extract_block(free_dofs, free_dofs);
	rhs_vec_free = rhs_vec(free_dofs);
}

void PoissonProblem::apply_robin_neumann_bc()
{
	free_dofs = mesh.extract_dof_inds(DOMAIN_ID);
	lhs_mat_free = lhs_mat.extract_block(free_dofs, free_dofs);
	rhs_vec_free = rhs_vec(free_dofs);
}

void PoissonProblem::apply_dirichlet_bc(
	int boundary_part, std::function<FloatType(FloatType, FloatType)> bc_func
) {

	free_dofs = mesh.extract_dof_inds(DOMAIN_ID & ~boundary_part);
	fixed_dofs = mesh.extract_dof_inds(boundary_part);

	lhs_mat_free = lhs_mat.extract_block(free_dofs, free_dofs);
	lhs_mat_fixed = lhs_mat.extract_block(free_dofs, fixed_dofs);

	bc_vec = Eigen::VectorX<FloatType>::Zero(fixed_dofs.size());
	int i = 0;
	for (int dof: fixed_dofs) {
		FloatType x = mesh.pmat(dof, 0);
		FloatType z = mesh.pmat(dof, 1);
		bc_vec(i) = bc_func(x, z);
		i++;
	}
	rhs_vec_free = rhs_vec(free_dofs) - lhs_mat_fixed*bc_vec;
}

// TODO: Make solve step backend agnostic
void PoissonProblem::solve_linear_system(std::string sp_solver_name)
{
	if (sp_solver_name.compare("sparselu") == 0) {
		Eigen::SparseLU<Eigen::SparseMatrix<FloatType>> sp_solver;
		sp_solver.analyzePattern(lhs_mat_free);
		sp_solver.factorize(lhs_mat_free);
		sol_vec_free = sp_solver.solve(rhs_vec_free);
	} 
	// else if (sp_solver_name.compare("umfpack") == 0) {
	// 	Eigen::UmfPackLU<Eigen::SparseMatrix<FloatType>> sp_solver;
	// 	sp_solver.analyzePattern(lhs_mat_free);
	// 	sp_solver.factorize(lhs_mat_free);
	// 	sol_vec_free = sp_solver.solve(rhs_vec_free);
	// }

	sol_vec(free_dofs) = sol_vec_free;
	sol_vec(fixed_dofs) = bc_vec;
}

FEMFunction2D PoissonProblem::solution()
{
	FEMFunction2D sol_function(mesh);
	sol_function.assign(sol_vec);
	return sol_function;
}
