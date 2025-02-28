//#include <Eigen/Dense>
//#include <Eigen/Sparse>
//#include <pstokes_fem.hpp>
//#include <free_surface_fem.hpp>
//#include <iostream>
//#include <cmath>
//
//enum {
//	DEFAULT,
//	BLACK = 30,
//	RED,
//	GREEN
//};
//
//template<typename T>
//void print_sparsity_pattern(Eigen::SparseMatrix<T> &sp_mat)
//{
//	for (int i = 0; i < sp_mat.rows(); i++) {
//		for (int j = 0; j < sp_mat.cols(); j++) {
//			T val = sp_mat.coeff(i, j);
//			if (fabs(val) < 1e-12) {
//				std::cout << "\x1b[" << RED << "m" << "0";
//			} else {
//				std::cout << "\x1b[" << GREEN << "m" << "1";
//			}
//		}
//		std::cout << "\n";
//	}
//	std::cout << "\x1b[" << DEFAULT << "m";
//}
//
//template<typename T>
//void get_sparse_nnz(Eigen::SparseMatrix<T> &sp_mat, std::vector<T> &nnz_ret)
//{
//	for (int i = 0; i < sp_mat.rows(); i++) {
//		for (int j = 0; j < sp_mat.cols(); j++) {
//			T val = sp_mat.coeff(i, j);
//			if (fabs(val) > 1e-12) {
//				nnz_ret.push_back(val);
//			}
//		}
//	}
//}
//
//int main(int argc, char *argv[])
//{
//	std::function<long double(long double, long double)> func = [](long double x, long double z) {
//		return 1.0;
//	};
//	std::function<long double(long double, long double)> a_func = [](long double x, long double z) {
//		return x;
//	};
//	StructuredMesh sm(2, 2, 2, TRIANGLE_RIGHT);
//	sm.corner_inds();
//	//FEMFunction h(&sm), ux(&sm), uz(&sm), as(&sm);
//	//h.assign(func);
//	//ux.assign(func);
//	//uz.assign(func);
//	//as.assign(a_func);
//	//FreeSurfaceProblem fp(h, ux, uz, as, &sm);
//	//fp.assemble_lhs(1.0, FORWARD_EULER, 3);
//	//fp.commit_lhs();
//	//fp.assemble_rhs(1.0, FORWARD_EULER, 3);	
//	//fp.solve("umfpack");	
//	//sm.extrude_z(fp.h_sol.vals);
//	//std::cout << sm.pmat << std::endl;
//	return 0;
//}
