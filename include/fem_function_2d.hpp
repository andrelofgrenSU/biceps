#pragma once
#include <structured_mesh.hpp>

class FEMFunction2D {
	private:

	public:
		explicit FEMFunction2D(StructuredMesh &mesh);

		StructuredMesh &mesh;
		Eigen::VectorX<FloatType> vals;

		void assign(std::function<FloatType(FloatType, FloatType)> func);
		void assign(const Eigen::VectorX<FloatType> &vec);
		void assign(FEMFunction2D &f);
		FloatType eval(FloatType x, FloatType z);
		Eigen::VectorX<FloatType> eval_cell(int cell_index);
		Eigen::VectorX<FloatType> eval_edge(int edge_index);
		Eigen::MatrixX<FloatType> extract_vertex_subvec(int domain_id);
		Eigen::MatrixX<FloatType> extract_dof_subvec(int domain_id);

		FEMFunction2D diff_x_interp();
		FEMFunction2D diff_z_interp();
		FEMFunction2D diff_x_proj(std::string sp_solver_name);
		FEMFunction2D diff_z_proj(std::string sp_solver_name);

		FloatType operator()(FloatType x, FloatType z);
		FEMFunction2D operator+(FloatType val);
		FEMFunction2D operator-(FloatType val);
		FEMFunction2D operator*(FloatType val);
		FEMFunction2D operator+(const FEMFunction2D &f);
		FEMFunction2D operator-(const FEMFunction2D &f);
		FEMFunction2D operator*(const FEMFunction2D &f);
};
