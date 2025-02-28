#pragma once
#include <float_type.hpp>
#include <vector>
#include <Eigen/Dense>

namespace MESH1D {
	enum DOMAIN_IDS {
		EMPTY_ID = 0,
		INTERIOR_ID = 1,
		WEST_ID = 2,
		EAST_ID = 4,
		BOUNDARY_ID = WEST_ID | EAST_ID,
		DOMAIN_ID = BOUNDARY_ID | INTERIOR_ID
	};
}

class IntervalMesh {
	private:
		int _degree, _nof_cells, _nof_verts, _nof_dofs;
	public:
		Eigen::MatrixX<FloatType> pmat;
		Eigen::MatrixXi cmat;
		Eigen::VectorXi dimat;

		IntervalMesh(FloatType x0, FloatType x1, int n_cells, int degree);
		IntervalMesh(const Eigen::VectorX<FloatType> &xvec, int n_cells, int degree);
		std::vector<int> extract_dof_inds(int id);
		int nof_cells();
		int nof_verts();
		int nof_dofs();
		int dofs_per_cell();
		int degree();
};
