#include <iostream>
#include <tuple>
#include <boost/format.hpp>
#include <poisson_test_suite.hpp>

BOOST_AUTO_TEST_CASE(test_manufactured_solution_case_1)
{
	populate_manufactured_sol_cases();	
	for (std::tuple<int, int, int, int, FloatType> test_case: test_manufactured_sol_cases) {
		int nx = std::get<0>(test_case);
		int nz = std::get<1>(test_case);
		int deg = std::get<2>(test_case);
		int cell_type = std::get<3>(test_case);
		FloatType error_tol = std::get<4>(test_case);
		PoissonProblem pp = PoissonProblem(nx, nz, deg, cell_type);
		pp.assemble_stiffness_block(&alpha_manufactured_case_2, &beta_manufactured_case_2, 3);
		pp.commit_lhs_mat();
		pp.assemble_force_rhs(&force_manufactured_case_1, 2);
		pp.apply_dirichlet_bc(BOUNDARY_ID, &bc_func_manufactured_case_1);
		pp.solve_linear_system("sparselu");
		Eigen::VectorX<FloatType> x_vec = pp.mesh.pmat(Eigen::all, 0);
		Eigen::VectorX<FloatType> z_vec = pp.mesh.pmat(Eigen::all, 1);

		Eigen::VectorX<FloatType> ue_vec = Eigen::VectorX<FloatType>::Zero(x_vec.size());
		for (int dof = 0; dof < x_vec.size(); dof++) {
			ue_vec(dof) = exact_sol_func_manufactured_case_1(
				x_vec(dof), z_vec(dof)
			);
		}
		Eigen::VectorX<FloatType> uh_vec = pp.sol_vec;
		FloatType error_inf_norm = (ue_vec - uh_vec).cwiseAbs().maxCoeff();
		BOOST_TEST_CHECK(
			error_inf_norm < error_tol,
			(boost::format{"Poisson solution error test with (alpha=1, beta=1, gamma=0) for StructuredMesh(nx=%d, nz=%d, degree=%d, cell_type=%s) failed: ||error_u||_inf = %g > %g"}
			% nx % nz % deg % cell_type %error_inf_norm %error_tol).str()
		);
	}
}

BOOST_AUTO_TEST_CASE(test_manufactured_solution_case_2)
{
	populate_manufactured_sol_cases();
	for (std::tuple<int, int, int, int, FloatType> test_case: test_manufactured_sol_cases) {
		int nx = std::get<0>(test_case);
		int nz = std::get<1>(test_case);
		int deg = std::get<2>(test_case);
		int cell_type = std::get<3>(test_case);
		FloatType error_tol = std::get<4>(test_case);
		PoissonProblem pp = PoissonProblem(nx, nz, deg, cell_type);
		pp.assemble_stiffness_block(&alpha_manufactured_case_2, &beta_manufactured_case_2, 2);
		pp.assemble_mass_block(&gamma_manufactured_case_2, 3);
		pp.commit_lhs_mat();
		pp.assemble_force_rhs(&force_manufactured_case_2, 2);
		pp.apply_dirichlet_bc(BOUNDARY_ID, &bc_func_manufactured_case_2);
		pp.solve_linear_system("sparselu");
		Eigen::VectorX<FloatType> x_vec = pp.mesh.pmat(Eigen::all, 0);
		Eigen::VectorX<FloatType> z_vec = pp.mesh.pmat(Eigen::all, 1);

		Eigen::VectorX<FloatType> ue_vec = Eigen::VectorX<FloatType>::Zero(x_vec.size());
		for (int dof = 0; dof < x_vec.size(); dof++) {
			ue_vec(dof) = exact_sol_func_manufactured_case_2(
				x_vec(dof), z_vec(dof)
			);
		}
		Eigen::VectorX<FloatType> uh_vec = pp.sol_vec;
		FloatType error_inf_norm = (ue_vec - uh_vec).cwiseAbs().maxCoeff();
		BOOST_TEST_CHECK(
			error_inf_norm < error_tol,
			(boost::format{"Poisson solution error test with (alpha=1, beta=1, gamma=1) for StructuredMesh(nx=%d, nz=%d, degree=%d, cell_type=%s) failed: ||error_u||_inf = %g > %g"}
			% nx % nz % deg % cell_type %error_inf_norm %error_tol).str()
		);
	}
}

BOOST_AUTO_TEST_CASE(test_manufactured_solution_case_3)
{
	populate_manufactured_sol_cases();
	for (std::tuple<int, int, int, int, FloatType> test_case: test_manufactured_sol_cases) {
		int nx = std::get<0>(test_case);
		int nz = std::get<1>(test_case);
		int deg = std::get<2>(test_case);
		int cell_type = std::get<3>(test_case);
		FloatType error_tol = std::get<4>(test_case);
		PoissonProblem pp = PoissonProblem(nx, nz, deg, cell_type);
		pp.assemble_stiffness_block(&alpha_manufactured_case_3, &beta_manufactured_case_3, 2);
		pp.assemble_mass_block(&gamma_manufactured_case_3, 3);
		pp.commit_lhs_mat();
		pp.assemble_force_rhs(&force_manufactured_case_3, 2);
		pp.apply_dirichlet_bc(BOUNDARY_ID, &bc_func_manufactured_case_3);
		pp.solve_linear_system("sparselu");
		Eigen::VectorX<FloatType> x_vec = pp.mesh.pmat(Eigen::all, 0);
		Eigen::VectorX<FloatType> z_vec = pp.mesh.pmat(Eigen::all, 1);

		Eigen::VectorX<FloatType> ue_vec = Eigen::VectorX<FloatType>::Zero(x_vec.size());
		for (int dof = 0; dof < x_vec.size(); dof++) {
			ue_vec(dof) = exact_sol_func_manufactured_case_3(
				x_vec(dof), z_vec(dof)
			);
		}
		Eigen::VectorX<FloatType> uh_vec = pp.sol_vec;
		FloatType error_inf_norm = (ue_vec - uh_vec).cwiseAbs().maxCoeff();
		BOOST_TEST_CHECK(
			error_inf_norm < error_tol,
			(boost::format{"Poisson solution error test (alpha=5, beta=3, gamma=-2) for StructuredMesh(nx=%d, nz=%d, degree=%d, cell_type=%s) failed: ||error_u||_inf = %g > %g"}
			% nx % nz % deg % cell_type %error_inf_norm %error_tol).str()
		);
	}
}

//BOOST_AUTO_TEST_CASE(test_lhs_mat_assembler, *boost::unit_test::tolerance(1e-12))
//{
//	populate_manufactured_sol_cases();	
//	for (std::tuple<int, int, int, int, FloatType> res: test_manufactured_sol_cases) {
//		int nx = std::get<0>(res);
//		int nz = std::get<1>(res);
//		int deg = std::get<2>(res);
//		int cell_type = std::get<3>(res);
//		FloatType error_ref = std::get<4>(res);
//		StructuredMesh mesh(nx, nz, deg, cell_type);
//		PoissonProblem pp = PoissonProblem(&mesh);
//		pp.assemble_lhs_mat(&alpha_manufactured, 3);
//		pp.assemble_rhs_vec(&force_manufactured, 2);
//		pp.apply_dirichlet_bc(BOUNDARY_ID, &bc_func_manufactured);
//		pp.solve_linear_system("umfpack");
//		Eigen::VectorX<FloatType> x_vec = mesh.pmat(Eigen::all, 0);
//		Eigen::VectorX<FloatType> z_vec = mesh.pmat(Eigen::all, 1);
//
//		Eigen::VectorX<FloatType> ue_vec = Eigen::VectorX<FloatType>::Zero(x_vec.size());
//		for (int dof = 0; dof < x_vec.size(); dof++) {
//			ue_vec(dof) = exact_sol_func_manufactured(
//				x_vec(dof), z_vec(dof)
//			);
//		}
//		Eigen::VectorX<FloatType> uh_vec = pp.sol_vec;
//		FloatType error_inf_norm = (ue_vec - uh_vec).cwiseAbs().maxCoeff();
//		BOOST_TEST_CHECK(
//			error_inf_norm < error_ref,
//			(boost::format{"Poisson solution error test for StructuredMesh(nx=%d, nz=%d, degree=%d, cell_type=%s) failed"}
//			% nx % nz % deg % cell_type).str()
//		);
//	}
//}
//
//BOOST_AUTO_TEST_CASE(test_rhs_vec_assembler, *boost::unit_test::tolerance(1e-12))
//{
//	populate_manufactured_sol_cases();	
//	for (std::tuple<int, int, int, int, FloatType> res: test_manufactured_sol_cases) {
//		int nx = std::get<0>(res);
//		int nz = std::get<1>(res);
//		int deg = std::get<2>(res);
//		int cell_type = std::get<3>(res);
//		FloatType error_ref = std::get<4>(res);
//		StructuredMesh mesh(nx, nz, deg, cell_type);
//		PoissonProblem pp = PoissonProblem(&mesh);
//		pp.assemble_lhs_mat(&alpha_manufactured, 3);
//		pp.assemble_rhs_vec(&force_manufactured, 2);
//		pp.apply_dirichlet_bc(BOUNDARY_ID, &bc_func_manufactured);
//		pp.solve_linear_system("umfpack");
//		Eigen::VectorX<FloatType> x_vec = mesh.pmat(Eigen::all, 0);
//		Eigen::VectorX<FloatType> z_vec = mesh.pmat(Eigen::all, 1);
//
//		Eigen::VectorX<FloatType> ue_vec = Eigen::VectorX<FloatType>::Zero(x_vec.size());
//		for (int dof = 0; dof < x_vec.size(); dof++) {
//			ue_vec(dof) = exact_sol_func_manufactured(
//				x_vec(dof), z_vec(dof)
//			);
//		}
//		Eigen::VectorX<FloatType> uh_vec = pp.sol_vec;
//		FloatType error_inf_norm = (ue_vec - uh_vec).cwiseAbs().maxCoeff();
//		BOOST_TEST_CHECK(
//			error_inf_norm < error_ref,
//			(boost::format{"Poisson solution error test for StructuredMesh(nx=%d, nz=%d, degree=%d, cell_type=%s) failed"}
//			% nx % nz % deg % cell_type).str()
//		);
//	}
//}
