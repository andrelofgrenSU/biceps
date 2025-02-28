#pragma once

#include <float_type.hpp>
#include <Eigen/Dense>
#include <structured_mesh.hpp>
#include <fem_function_2d.hpp>

namespace Eigen {
	typedef Eigen::MatrixX<FloatType> MatrixXld;
	typedef Eigen::VectorX<FloatType> VectorXld;
}

namespace FEM2D {

	void gauss_legendre_quadrature(
		const int precision,
		const int cell_type,
		Eigen::MatrixX<FloatType> &points_ret,
		Eigen::VectorX<FloatType> &weights_ret
	);

	void gauss_legendre_quadrature_edge_quadrilateral(
		const int precision,
		const int edge_id,
		Eigen::MatrixX<FloatType> &points_ret,
		Eigen::VectorX<FloatType> &weights_ret
	);

	void map_rs_to_xz(
		Eigen::MatrixX<FloatType> &node_coords,
		Eigen::MatrixX<FloatType> &qpoints_rs,
		Eigen::MatrixX<FloatType> &qpoints_xz_ret
	);

	void map_xz_to_rs(
		Eigen::MatrixX<FloatType> &node_coords,
		Eigen::MatrixX<FloatType> &qpoints_xz,
		Eigen::MatrixX<FloatType> &qpoints_rs_ret
	);

	void lagrange_basis(
		const int degree,
		const int cell_type,
		Eigen::MatrixX<FloatType> &qpoints_rs,
		Eigen::MatrixX<FloatType> &phi_rs_ret,
		Eigen::MatrixX<FloatType> &grad_phi_rs_ret
	);

	void map_to_reference_cell(
		const int degree, 
		const int cell_type,
		const Eigen::MatrixX<FloatType> &node_coords,
		const Eigen::MatrixX<FloatType> &qpoints_rs,
		const Eigen::MatrixX<FloatType> &phi_rs,
		const Eigen::MatrixX<FloatType> &grad_phi_rs,
		Eigen::VectorX<FloatType> &detJ_rs_ret,
		Eigen::MatrixX<FloatType> &qpoints_xz_ret,
		Eigen::MatrixX<FloatType> &grad_phi_xz_ret
	);

	Eigen::MatrixX<FloatType> reference_element_points_rs(int cell_type, int degree);

	FloatType L2_norm(const FEMFunction2D &u, const Eigen::SparseMatrix<FloatType> &M_sp);
	FloatType L2_norm(const Eigen::VectorX<FloatType> &u_vec, const Eigen::SparseMatrix<FloatType> &M_sp);

	FloatType calculate_area(StructuredMesh &mesh, int gp);

	Eigen::SparseMatrix<FloatType> assemble_mass_matrix(
		StructuredMesh &mesh, int gp
	);
	Eigen::SparseMatrix<FloatType> assemble_stiffness_matrix(
		StructuredMesh &mesh, int gp
	);
	Eigen::SparseMatrix<FloatType> assemble_stiffness_xx_matrix(
		StructuredMesh &mesh, int gp
	);
	Eigen::SparseMatrix<FloatType> assemble_stiffness_zz_matrix(
		StructuredMesh &mesh, int gp
	);
	Eigen::SparseMatrix<FloatType> assemble_stiffness_x_matrix(
		StructuredMesh &mesh, int gp
	);
	Eigen::SparseMatrix<FloatType> assemble_stiffness_z_matrix(
		StructuredMesh &mesh, int gp
	);
	Eigen::MatrixX<FloatType> assemble_expansion_matrix(
		StructuredMesh &mesh,
		std::function<FloatType(FloatType, FloatType)> force,
		int gp
	);

	Eigen::SparseMatrix<FloatType> assemble_stress_matrix(
		StructuredMesh &mesh_u,
		StructuredMesh &mesh_p,
		FEMFunction2D ux,
		FEMFunction2D uz,
		FloatType A,
		FloatType n_i,
		int gp
	);

	Eigen::VectorX<FloatType> assemble_load_vector(
		StructuredMesh &mesh, std::function<FloatType(FloatType, FloatType)> f, int gp
	);
}
