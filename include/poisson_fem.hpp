#pragma once

#include <vector>
#include <functional>
#include <fem_function_2d.hpp>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_spmat_addons.hpp>
#include <Eigen/Sparse>

class PoissonProblem {

	private:
		int n_dofs;

	public:
		std::vector<int> free_dofs;
		std::vector<int> fixed_dofs;

		PoissonProblem(int nx, int nz, int deg, int cell_type);
		
		StructuredMesh mesh;
		std::vector<Eigen::Triplet<FloatType>> lhs_coeffs;

		Eigen::SparseMatrix<FloatType> lhs_mat;
		Eigen::SparseMatrix<FloatType> lhs_mat_free;
		Eigen::SparseMatrix<FloatType> lhs_mat_fixed;
		Eigen::VectorX<FloatType> bc_vec;
		Eigen::VectorX<FloatType> rhs_vec;
		Eigen::VectorX<FloatType> rhs_vec_free;
		Eigen::VectorX<FloatType> sol_vec_free;
		Eigen::VectorX<FloatType> sol_vec;

		void reset();
		void assemble_stiffness_block(
			std::function<FloatType(FloatType, FloatType)> alpha,
			std::function<FloatType(FloatType, FloatType)> beta,
			int gp
		);
		void assemble_mass_block(
			std::function<FloatType(FloatType, FloatType)> beta,
			int gauss_precision
		);
		void assemble_FSSA_block(
			FEMFunction2D &f,
			FloatType theta,
			FloatType dt,
			int gp
		);
		void assemble_robin_block(
			std::function<FloatType(FloatType, FloatType)> a_robin,
			std::function<FloatType(FloatType, FloatType)> b_robin,
			std::function<FloatType(FloatType, FloatType)> g_robin,
			int boundary_id,
			int gauss_precision
		);
		void commit_lhs_mat();

		void assemble_force_rhs(
			std::function<FloatType(FloatType, FloatType)> f,
			int gauss_precision
		);

		void assemble_force_rhs(
			FEMFunction2D &f,
			int gauss_precision
		);

		void assemble_neumann_rhs(
			std::function<FloatType(FloatType, FloatType)> g_neumann,
			int boundary_id,
			int gauss_precision
		);

		void assemble_robin_rhs(
			std::function<FloatType(FloatType, FloatType)> b_robin,
			std::function<FloatType(FloatType, FloatType)> g_robin,
			int boundary_id,
			int gauss_precision
		);

		void apply_zero_dirichlet_bc(int boundary_id);
		void apply_dirichlet_bc(
			int boundary_id,
			std::function<FloatType(FloatType, FloatType)> g_dirichlet
		);
		void apply_robin_neumann_bc();

		FEMFunction2D solution();

		void solve_linear_system(std::string sp_solver_name);
		void reset_lhs();
};
