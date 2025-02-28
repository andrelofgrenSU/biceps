#pragma once
#include <pstokes_fem.hpp>
#include <free_surface_fem.hpp>

class TimeIntegrator {
	private:
		std::vector<int> svinds;

	public:
		pStokesProblem &psp;
		FreeSurfaceProblem &fsp;
		TimeIntegrator(pStokesProblem &psp, FreeSurfaceProblem &fsp);
		Eigen::VectorX<FloatType> step_explicit(FloatType dt, bool reset_pstokes = true);
		Eigen::VectorX<FloatType> step_semi_implicit(FloatType dt, bool reset_pstokes = true);
		Eigen::VectorX<FloatType> step_rk2(FloatType dt, bool reset_pstokes = true);
		Eigen::MatrixX<FloatType> transient_jacobian(int time_scheme, int jac_stencil, FloatType dt_eps);
		void extrude_mesh_z(const Eigen::VectorX<FloatType> &zs_vec);

// 		Eigen::VectorX<FloatType> step_rk4(
// 			FloatType dt,
// 			FloatType picard_ps_rtol,
// 			int picard_ps_max_iter,
// 			int gp_ps,
// 			int gp_fs
// 		);
// 
// 		Eigen::VectorX<FloatType> step_semi_implicit(
// 			FloatType dt,
// 			int fssa_version,
// 			FloatType fssa_param,
// 			FloatType picard_ps_rtol,
// 			int picard_ps_max_iter,
// 			int gp_ps,
// 			int gp_fs,
// 			bool reset_system
// 		);
// 
// 		Eigen::VectorX<FloatType> step_midpoint(
// 			FloatType dt,
// 			FloatType picard_ps_rtol,
// 			FloatType picard_fs_atol,
// 			int picard_ps_max_iter,
// 			int picard_fs_max_iter,
// 			int gp_ps,
// 			int gp_fs
// 		);
// 
// 		Eigen::VectorX<FloatType> step_crank_nicolson(
// 			FloatType dt,
// 			FloatType picard_ps_rtol,
// 			FloatType picard_fs_atol,
// 			int picard_ps_max_iter,
// 			int picard_fs_max_iter,
// 			int gp_ps,
// 			int gp_fs
// 		);
};
