#pragma once
#include <interval_mesh.hpp>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_spmat_addons.hpp>
#include <Eigen/Sparse>

namespace FEM1D {
	void map_to_reference_cell(
		int degree, 
		Eigen::MatrixX<FloatType> &node_coords,
		Eigen::VectorX<FloatType> &qpoints_r,
		Eigen::MatrixX<FloatType> &phi_r,
		Eigen::MatrixX<FloatType> &grad_phi_r,
		Eigen::VectorX<FloatType> &detJ_r_ret,
		Eigen::MatrixX<FloatType> &qpoints_x_ret,
		Eigen::MatrixX<FloatType> &grad_phi_x_ret
	);

    Eigen::MatrixX<FloatType> reference_element_points_rs(
    	const int degree
    );

	void lagrange_basis(
		const int degree,
		Eigen::VectorX<FloatType> &qpoints_r,
		Eigen::MatrixX<FloatType> &phi_r_ret,
		Eigen::MatrixX<FloatType> &grad_phi_r_ret
	);

	void gauss_legendre_quadrature(
		const int precision,
		Eigen::VectorX<FloatType> &points,
		Eigen::VectorX<FloatType> &weights
	);

	Eigen::SparseMatrix<FloatType> assemble_mass_matrix(IntervalMesh &mesh, int gp);
}
