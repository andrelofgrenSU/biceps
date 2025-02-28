#include <interval_mesh.hpp>

IntervalMesh::IntervalMesh(FloatType x0, FloatType x1, int n_cells, int degree) :
	_degree(degree), _nof_cells(n_cells)
{
	pmat = Eigen::MatrixX<FloatType>::Zero(nof_dofs(), 2);
	cmat = Eigen::MatrixXi::Zero(nof_cells(), dofs_per_cell());
	dimat = Eigen::VectorXi::Zero(nof_dofs());

	for (int di = 0; di < nof_dofs(); di++) {
		pmat(di, 0) = x0 + (x1-x0)*di/(nof_dofs()-1);
		if (ABS_FUNC(pmat(di, 0) - x0) < 1e-10)
			dimat(di) = MESH1D::WEST_ID;
		else if (ABS_FUNC(pmat(di, 0) - x1) < 1e-10)
			dimat(di) = MESH1D::EAST_ID;
		else
			dimat(di) = MESH1D::INTERIOR_ID;
	}

	for (int ci = 0; ci < nof_cells(); ci++)
		for (int k = 0; k < dofs_per_cell(); k++)
			cmat(ci, k) = ci*degree+k;
}

IntervalMesh::IntervalMesh(
	const Eigen::VectorX<FloatType> &xvec, int n_cells, int degree
) :
	_degree(degree), _nof_cells(n_cells)
{
	pmat = Eigen::MatrixX<FloatType>::Zero(nof_dofs(), 2);
	cmat = Eigen::MatrixXi::Zero(nof_cells(), dofs_per_cell());
	dimat = Eigen::VectorXi::Zero(nof_dofs());
	if (xvec.size() != nof_dofs())
		throw std::length_error("Size of xvec does not match the number of dofs");
	FloatType x0 = xvec(0);
	FloatType x1 = xvec(Eigen::last);
	for (int di = 0; di < nof_dofs(); di++) {
		pmat(di, 0) = xvec(di);
		if (ABS_FUNC(pmat(di, 0) - x0) < 1e-10)
			dimat(di) = MESH1D::WEST_ID;
		else if (ABS_FUNC(pmat(di, 0) - x1) < 1e-10)
			dimat(di) = MESH1D::EAST_ID;
		else
			dimat(di) = MESH1D::INTERIOR_ID;
	}

	for (int ci = 0; ci < nof_cells(); ci++)
		for (int k = 0; k < dofs_per_cell(); k++)
			cmat(ci, k) = ci*degree+k;
}

std::vector<int> IntervalMesh::extract_dof_inds(int id)
{
	std::vector<int> dinds;
	for (int di = 0; di < dimat.rows(); di++)
		if (dimat(di) & id)
			dinds.push_back(di);
	return dinds;
}

int IntervalMesh::nof_cells()
{
	return _nof_cells;
}

int IntervalMesh::nof_verts()
{
	return _nof_cells + 1;
}

int IntervalMesh::nof_dofs()
{
	return degree()*nof_cells()+1;
}

int IntervalMesh::dofs_per_cell()
{
	return _degree + 1;
}

int IntervalMesh::degree()
{
	return _degree;
}
