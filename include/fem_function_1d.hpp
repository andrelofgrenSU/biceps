#pragma once
#include <interval_mesh.hpp>

class FEMFunction1D {
	private:
		
	public:
		IntervalMesh &mesh;
		Eigen::VectorX<FloatType> vals;

		explicit FEMFunction1D(IntervalMesh &mesh);

		void assign(std::function<FloatType(FloatType, FloatType)> func);
		void assign(const Eigen::VectorX<FloatType> &vec);
		FloatType integrate(int gp);
		Eigen::VectorX<FloatType> eval_cell(int cell_index);
		Eigen::MatrixX<FloatType> extract_vertex_subvec(int domain_id);
		Eigen::MatrixX<FloatType> extract_dof_subvec(int domain_id);

		// TODO: Make operators work also for const FEMFunction1D
		FloatType operator()(FloatType x, FloatType z);
		FEMFunction1D operator+(FloatType b);
		FEMFunction1D operator-(FloatType b);
		FEMFunction1D operator*(FloatType b);
		FEMFunction1D operator+(const FEMFunction1D &b);
		FEMFunction1D operator-(const FEMFunction1D &b);
		FEMFunction1D operator*(const FEMFunction1D &b);
};
