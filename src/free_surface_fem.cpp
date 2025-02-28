#include <free_surface_fem.hpp>
#include <fem_1d.hpp>

FreeSurfaceProblem::FreeSurfaceProblem(
	IntervalMesh &h_mesh, IntervalMesh &u_mesh
) : h_mesh(h_mesh), u_mesh(u_mesh), ac_fem_func(h_mesh)
{
	lhs_mat = Eigen::SparseMatrix<FloatType>(
		h_mesh.nof_dofs(), h_mesh.nof_dofs()
	);
	rhs_vec = Eigen::VectorX<FloatType>::Zero(h_mesh.nof_dofs());
	zs_vec = h_mesh.pmat(Eigen::all, 1);
}

void FreeSurfaceProblem::set_accumulation(const Eigen::VectorX<FloatType> &ac_vec)
{
	ac_fem_func.vals = ac_vec;
}

void FreeSurfaceProblem::assemble_lhs_explicit() {
	lhs_coeffs.clear();
	lhs_coeffs.reserve((2*h_mesh.degree()+1)*h_mesh.nof_dofs());

	Eigen::MatrixX<FloatType> node_coords, qpoints_x, phi_r, dphi_r, dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi element;
	
	FEM1D::gauss_legendre_quadrature(
		gp_lhs, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		h_mesh.degree(), qpoints_r, phi_r, dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	dphi_x = Eigen::MatrixX<FloatType>::Zero(
		qpoints_r.rows(), h_mesh.dofs_per_cell()
	);

	Eigen::MatrixX<FloatType> A = Eigen::MatrixX<FloatType>::Zero(
		h_mesh.dofs_per_cell(), h_mesh.dofs_per_cell()
	);

	for (int ci = 0; ci < h_mesh.nof_cells(); ci++) {
		element = h_mesh.cmat(ci, Eigen::all);
		node_coords = h_mesh.pmat(element, Eigen::all);
		FEM1D::map_to_reference_cell(
			h_mesh.degree(), node_coords, qpoints_r,
			phi_r, dphi_r, detJ_r, qpoints_x, dphi_x
		);

		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType detJxW = qweights(q)*detJ_r(q);
			for (int i = 0; i < h_mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < h_mesh.dofs_per_cell(); j++) {
					A(i, j) += phi_r(q, i)*phi_r(q, j)*detJxW;
				}
			}
		}

		for (int i = 0; i < h_mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < h_mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					element(i),
					element(j),
					A(i, j)
				));
			}
		}	
		A.setZero();
	}
}

void FreeSurfaceProblem::assemble_rhs_explicit(
	FEMFunction1D &h_fem_func,
	FEMFunction1D &ux_fem_func,
	FEMFunction1D &uz_fem_func,
	FloatType dt
) {
	rhs_vec.setZero();

	Eigen::MatrixX<FloatType> node_coords_h, node_coords_u, qpoints_x, h_phi_r, h_dphi_r,
		h_dphi_x, u_phi_r, u_dphi_r, u_dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi element_h, element_u;

	FEM1D::gauss_legendre_quadrature(
		gp_rhs, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		h_mesh.degree(), qpoints_r, h_phi_r, h_dphi_r
	);

	FEM1D::lagrange_basis(
		u_mesh.degree(), qpoints_r, u_phi_r, u_dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	h_dphi_x = Eigen::MatrixX<FloatType>::Zero(
		qpoints_r.rows(), h_mesh.dofs_per_cell()
	);
	u_dphi_x = Eigen::MatrixX<FloatType>::Zero(
		qpoints_r.rows(), u_mesh.dofs_per_cell()
	);

	for (int ci = 0; ci < h_mesh.nof_cells(); ci++) {
		element_h = h_mesh.cmat(ci, Eigen::all);
		element_u = u_mesh.cmat(ci, Eigen::all);
		node_coords_h = h_mesh.pmat(element_h, Eigen::all);
		node_coords_u = u_mesh.pmat(element_u, Eigen::all);
		FEM1D::map_to_reference_cell(
			u_mesh.degree(), node_coords_u, qpoints_r,
			u_phi_r, u_dphi_r, detJ_r, qpoints_x, u_dphi_x
		);
		FEM1D::map_to_reference_cell(
			h_mesh.degree(), node_coords_h, qpoints_r,
			h_phi_r, h_dphi_r, detJ_r, qpoints_x, h_dphi_x
		);

		Eigen::VectorX<FloatType> h_vec = h_phi_r*h_fem_func.eval_cell(ci);
		Eigen::VectorX<FloatType> dh_vec = h_dphi_x*h_fem_func.eval_cell(ci);
		Eigen::VectorX<FloatType> ux_vec = u_phi_r*ux_fem_func.eval_cell(ci);
		Eigen::VectorX<FloatType> uz_vec = u_phi_r*uz_fem_func.eval_cell(ci);
		Eigen::VectorX<FloatType> ac_vec = h_phi_r*ac_fem_func.eval_cell(ci);
		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType detJxW = qweights(q)*detJ_r(q);
			for (int i = 0; i < h_mesh.dofs_per_cell(); i++) {
				rhs_vec(element_h(i)) += h_phi_r(q, i)*(
					h_vec(q) + dt*(uz_vec(q) - ux_vec(q)*dh_vec(q) + ac_vec(q))
				)*detJxW;
			}
		}
	}
}

void FreeSurfaceProblem::assemble_lhs_semi_implicit(
	FEMFunction1D &ux_fem_func, FloatType dt
) {
	lhs_coeffs.clear();
	lhs_coeffs.reserve((2*h_mesh.degree()+1)*h_mesh.nof_dofs());

	Eigen::MatrixX<FloatType> node_coords_h, node_coords_u, qpoints_x, h_phi_r, h_dphi_r,
		h_dphi_x, u_phi_r, u_dphi_r, u_dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi element_h, element_u;

	FEM1D::gauss_legendre_quadrature(
		gp_lhs, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		h_mesh.degree(), qpoints_r, h_phi_r, h_dphi_r
	);

	FEM1D::lagrange_basis(
		u_mesh.degree(), qpoints_r, u_phi_r, u_dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	h_dphi_x = Eigen::MatrixX<FloatType>::Zero(
		qpoints_r.rows(), h_mesh.dofs_per_cell()
	);
	u_dphi_x = Eigen::MatrixX<FloatType>::Zero(
		qpoints_r.rows(), u_mesh.dofs_per_cell()
	);

	Eigen::MatrixX<FloatType> A = Eigen::MatrixX<FloatType>::Zero(
		h_mesh.dofs_per_cell(), h_mesh.dofs_per_cell()
	);

	for (int ci = 0; ci < h_mesh.nof_cells(); ci++) {
		element_h = h_mesh.cmat(ci, Eigen::all);
		element_u = u_mesh.cmat(ci, Eigen::all);
		node_coords_h = h_mesh.pmat(element_h, Eigen::all);
		node_coords_u = u_mesh.pmat(element_u, Eigen::all);
		FEM1D::map_to_reference_cell(
			u_mesh.degree(), node_coords_u, qpoints_r,
			u_phi_r, u_dphi_r, detJ_r, qpoints_x, u_dphi_x
		);
		FEM1D::map_to_reference_cell(
			h_mesh.degree(), node_coords_h, qpoints_r,
			h_phi_r, h_dphi_r, detJ_r, qpoints_x, h_dphi_x
		);
		Eigen::VectorX<FloatType> ux_vec = u_phi_r*ux_fem_func.eval_cell(ci);

		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType detJxW = qweights(q)*detJ_r(q);
			for (int i = 0; i < h_mesh.dofs_per_cell(); i++) {
				for (int j = 0; j < h_mesh.dofs_per_cell(); j++) {
					A(i, j) += (
						h_phi_r(q, i)*h_phi_r(q, j) +
						dt*ux_vec(q)*h_phi_r(q, i)*h_dphi_x(q, j)
					)*detJxW;
				}
			}
		}

		for (int i = 0; i < h_mesh.dofs_per_cell(); i++) {
			for (int j = 0; j < h_mesh.dofs_per_cell(); j++) {
				lhs_coeffs.push_back(Eigen::Triplet<FloatType>(
					element_h(i),
					element_h(j),
					A(i, j)
				));
			}
		}
		A.setZero();
	}
}

void FreeSurfaceProblem::assemble_rhs_semi_implicit(
	FEMFunction1D &h_fem_func,
	FEMFunction1D &uz_fem_func,
	FloatType dt
) {
	rhs_vec.setZero();
	Eigen::MatrixX<FloatType> node_coords_h, node_coords_u, qpoints_x, h_phi_r, h_dphi_r,
		h_dphi_x, u_phi_r, u_dphi_r, u_dphi_x;
	Eigen::VectorX<FloatType> qweights, qpoints_r, detJ_r;
	Eigen::VectorXi element_h, element_u;

	FEM1D::gauss_legendre_quadrature(
		gp_rhs, qpoints_r, qweights
	);

	FEM1D::lagrange_basis(
		h_mesh.degree(), qpoints_r, h_phi_r, h_dphi_r
	);

	FEM1D::lagrange_basis(
		u_mesh.degree(), qpoints_r, u_phi_r, u_dphi_r
	);

	qpoints_x = Eigen::MatrixX<FloatType>::Zero(qpoints_r.rows(), 2);
	detJ_r = Eigen::VectorX<FloatType>::Zero(qpoints_r.rows());
	h_dphi_x = Eigen::MatrixX<FloatType>::Zero(
		qpoints_r.rows(), h_mesh.dofs_per_cell()
	);
	u_dphi_x = Eigen::MatrixX<FloatType>::Zero(
		qpoints_r.rows(), u_mesh.dofs_per_cell()
	);

	for (int ci = 0; ci < h_mesh.nof_cells(); ci++) {
		element_h = h_mesh.cmat(ci, Eigen::all);
		element_u = u_mesh.cmat(ci, Eigen::all);
		node_coords_h = h_mesh.pmat(element_h, Eigen::all);
		node_coords_u = u_mesh.pmat(element_u, Eigen::all);
		FEM1D::map_to_reference_cell(
			u_mesh.degree(), node_coords_u, qpoints_r,
			u_phi_r, u_dphi_r, detJ_r, qpoints_x, u_dphi_x
		);
		FEM1D::map_to_reference_cell(
			h_mesh.degree(), node_coords_h, qpoints_r,
			h_phi_r, h_dphi_r, detJ_r, qpoints_x, h_dphi_x
		);

		Eigen::VectorX<FloatType> h_vec = h_phi_r*h_fem_func.eval_cell(ci);
		Eigen::VectorX<FloatType> uz_vec = u_phi_r*uz_fem_func.eval_cell(ci);
		Eigen::VectorX<FloatType> ac_vec = h_phi_r*ac_fem_func.eval_cell(ci);
		for (int q = 0; q < qpoints_r.rows(); q++) {
			FloatType detJxW = qweights(q)*detJ_r(q);
			for (int i = 0; i < h_mesh.dofs_per_cell(); i++) {
				rhs_vec(element_h(i)) += h_phi_r(q, i)*(
					h_vec(q) + dt*(uz_vec(q) + ac_vec(q))
				)*detJxW;
			}
		}
	}
}

void FreeSurfaceProblem::commit_lhs()
{
	lhs_mat.setFromTriplets(lhs_coeffs.begin(), lhs_coeffs.end());
}

void FreeSurfaceProblem::solve_linear_system()
{
	if (linear_solver == SPARSE_LU) {
		Eigen::SparseLU<Eigen::SparseMatrix<FloatType>> solver;
		solver.analyzePattern(lhs_mat);
		solver.factorize(lhs_mat);
		zs_vec = solver.solve(rhs_vec);
	} else if (linear_solver == BICGSTAB) {
		Eigen::BiCGSTAB<Eigen::SparseMatrix<FloatType>> solver;
		solver.compute(lhs_mat);
		zs_vec = solver.solve(rhs_vec);
	} else {
		// TODO: throw exception
	}
}
