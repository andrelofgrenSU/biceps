#pragma once
#include <enums.hpp>
#include <logger.hpp>
#include <fem_function_2d.hpp>

class pStokesProblem {

	private:
		int nv_dofs;
		int np_dofs;
		int n_dofs;

		Eigen::VectorXi ux_constrained_dofs;
		Eigen::VectorXi uz_constrained_dofs;
		int boundary_ux;
		int boundary_uz;

	public:
		Logger logger;
		FloatType A, n_i, eps_reg_2, fssa_param = 0;
		int gp_stress = 5;
		int gp_incomp = 5;
		int gp_fssa = 5;
		int gp_rhs = 5;

		bool assemble_stress = true;
		bool assemble_incomp = true;
		int fssa_version = FSSA_NONE;
		int linear_solver = SPARSE_LU;
		int picard_max_iter = 100;
		FloatType picard_atol = 1e-10;
		FloatType prune_threshold = 1e-12;

		std::function<FloatType(FloatType, FloatType)> force_x;
		std::function<FloatType(FloatType, FloatType)> force_z;
		std::function<FloatType(FloatType, FloatType)> fssa_accum = [](FloatType x, FloatType z) {return 0.0;};

		std::vector<Eigen::Triplet<FloatType>> lhs_coeffs;
		int nof_pushed_elements = 0;

		StructuredMesh &u_mesh; 
		StructuredMesh &p_mesh;

		Eigen::VectorXi ux_v2d;
		Eigen::VectorXi uz_v2d;
		Eigen::VectorXi u_v2d;
		Eigen::VectorXi p_v2d;

		Eigen::VectorXi ux_d2v;
		Eigen::VectorXi uz_d2v;
		Eigen::VectorXi p_d2v;	
		Eigen::VectorXi w_d2v;

		Eigen::MatrixX<FloatType> pmat_p1;
		Eigen::MatrixXi cmat_p1;

		Eigen::SparseMatrix<FloatType> lhs_mat;
		Eigen::VectorX<FloatType> rhs_vec;
		Eigen::VectorX<FloatType> w_vec;

		pStokesProblem(
			FloatType rate_factor,
			FloatType glen_exp,
			FloatType eps_reg_2,
			std::function<FloatType(FloatType, FloatType)> force_x,
			std::function<FloatType(FloatType, FloatType)> force_z,
			StructuredMesh &u_mesh,
			StructuredMesh &p_mesh
		);

		void assemble_stress_block();
		void assemble_incomp_block();
		void assemble_fssa_vertical_block();
		void assemble_fssa_normal_block();
		void commit_lhs_mat(); 

		void assemble_rhs_vec();
		void assemble_fssa_vertical_rhs();

		void apply_standard_dirichlet_bc();
		void apply_dirichlet_bc(
			int boundary_part, int component, std::function<FloatType(FloatType, FloatType)> u_func
		);
		void prune_lhs();

		void solve_linear_system();
		void solve_picard(bool reset_linear_system = true);

		FEMFunction2D velocity_x();
		FEMFunction2D velocity_z();
		FEMFunction2D pressure();

		void reset_lhs();
		void reset_rhs();
		void reset_system();
};
