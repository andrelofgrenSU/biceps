#pragma once

#include <fem_function_1d.hpp>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_spmat_addons.hpp>
#include <Eigen/Sparse>
#include <enums.hpp>

class FreeSurfaceProblem {
	private:
		std::vector<Eigen::Triplet<FloatType>> lhs_coeffs;
		Eigen::SparseMatrix<FloatType> lhs_mat;
		Eigen::VectorX<FloatType> rhs_vec;

	public:
		int gp_lhs = 5;
		int gp_rhs = 5;
		int linear_solver = SPARSE_LU;

		FEMFunction1D ac_fem_func;
		IntervalMesh &h_mesh;
		IntervalMesh &u_mesh;
		Eigen::VectorX<FloatType> zs_vec;

		FreeSurfaceProblem(IntervalMesh &h_mesh, IntervalMesh &u_mesh);

		void set_accumulation(const Eigen::VectorX<FloatType> &ac_vec);
		void assemble_lhs_explicit();
		void assemble_rhs_explicit(
			FEMFunction1D &h,
			FEMFunction1D &ux,
			FEMFunction1D &uz,
			FloatType dt
		);
		void assemble_lhs_semi_implicit(
			FEMFunction1D &ux,
			FloatType dt
		);
		void assemble_rhs_semi_implicit(
			FEMFunction1D &h,
			FEMFunction1D &uz,
			FloatType dt
		);
		void commit_lhs();
		void solve_linear_system();
};
