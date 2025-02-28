#include <functional>
#include <tuple>
#include <set>
#include <boost/timer/timer.hpp>
#include <structured_mesh.hpp>

MeshGraph::MeshGraph(int nx, int nz, int degree)
{
	this->nx = nx;
	this->nz = nz;
	this->_degree = degree;

	hl = _degree*nx;
	vl = _degree*nz;
	_nof_verts = (nx+1)*(nz+1);
	_nof_dofs = (hl+1)*(vl+1);
	dof_adj_mat = Eigen::SparseMatrix<int>(_nof_dofs, _nof_dofs);

	v2d = Eigen::VectorXi(_nof_dofs);
	d2v = Eigen::VectorXi(_nof_dofs);

	assemble_dof_adj_mat();
}

void MeshGraph::assemble_dof_adj_mat()
{
	int nnz_adj_dof_mat = 2*(2*hl*vl + hl + vl); 
	std::vector<Eigen::Triplet<int>> adj_dof_coeffs;
	adj_dof_coeffs.reserve(nnz_adj_dof_mat);

	int dof_index = 0;
	for (int i = 0; i < vl+1; i++) {
		for (int j = 0; j < hl+1; j++) {
			if (i < vl)  {
				adj_dof_coeffs.push_back(
					Eigen::Triplet<int>(dof_index, (i+1)*(hl+1) + j, NORTH_ID)
				);
			}

			if (i < vl && j > 0)
				adj_dof_coeffs.push_back(
					Eigen::Triplet<int>(dof_index, (i+1)*(hl+1) + j-1, NORTH_WEST_ID)
				);
			
			if (j > 0)
				adj_dof_coeffs.push_back(
					Eigen::Triplet<int>(dof_index, i*(hl+1) + j-1, WEST_ID)
				);

			if (i > 0 && j > 0)
				adj_dof_coeffs.push_back(
					Eigen::Triplet<int>(dof_index, (i-1)*(hl+1) + j-1, SOUTH_WEST_ID)
				);
			
			if (i > 0)
				adj_dof_coeffs.push_back(
					Eigen::Triplet<int>(dof_index, (i-1)*(hl+1) + j, SOUTH_ID)
				);

			if (i > 0 && j < hl)
				adj_dof_coeffs.push_back(
					Eigen::Triplet<int>(dof_index, (i-1)*(hl+1) + j+1, SOUTH_EAST_ID)
				);
			
			if (j < hl)
				adj_dof_coeffs.push_back(
					Eigen::Triplet<int>(dof_index, i*(hl+1) + j+1, EAST_ID)
				);

			if (j < hl && i < vl) {
				adj_dof_coeffs.push_back(
					Eigen::Triplet<int>(dof_index, (i+1)*(hl+1) + j+1, NORTH_EAST_ID)
				);
			}

			d2v[dof_index] = dof_index;
			v2d[dof_index] = dof_index;
			dof_index++;	
		}
	}
	dof_adj_mat.setFromTriplets(
		adj_dof_coeffs.begin(), adj_dof_coeffs.end()
	);
}

int MeshGraph::vertex_degree(int vi)
{
	Eigen::SparseVector<int> vertex = dof_adj_mat.row(vi);
	return vertex.nonZeros();
}

int MeshGraph::nof_verts() {
	return _nof_verts;
}

int MeshGraph::nof_dofs() {
	return _nof_dofs;
}

int MeshGraph::degree() {
	return _degree;
}

StructuredMesh::StructuredMesh(int nx, int nz, int degree, int cell_type) : 
	MeshGraph(nx, nz, degree), _cell_type(cell_type)
{
	if (cell_type == TRIANGLE_LEFT || cell_type == TRIANGLE_RIGHT) {
		_dofs_per_cell = (_degree+1)*(_degree+2) >> 1;
		_nof_cells = 2*nx*nz;
		_nof_edges = 3*nx*nz + nx + nz;
		_dofs_per_edge = degree+1;
		_edges_per_cell = 3;
	} else if (cell_type == QUADRILATERAL) {
		_dofs_per_cell = (_degree+1)*(_degree+1);
		_nof_cells = nx*nz;
		_dofs_per_edge = degree+1;
		_nof_edges = 2*nx*nz + nx + nz;
		_edges_per_cell = 4;
	}

	pmat = Eigen::MatrixX<FloatType>(_nof_dofs, 2);
	pmat_unit_box = Eigen::MatrixX<FloatType>(_nof_dofs, 2);
	cmat = Eigen::MatrixXi(_nof_cells, _dofs_per_cell);
	dimat = Eigen::VectorXi(_nof_dofs);
	cimat = Eigen::VectorXi(_nof_cells);

	assemble_pmat();
	assemble_cmat();
	assemble_emat();
	assemble_cimat();
	assemble_dimat();
	assemble_eimat();
	compute_vertex_dof_inds();
	compute_normals_and_tangents();

	v2s_map = Eigen::VectorXi::Zero(_nof_dofs);
	v2b_map = Eigen::VectorXi::Zero(_nof_dofs);
	for (int vi = 0; vi < _nof_dofs; vi++) {
		v2b_map[vi] = vi % (hl+1);
		v2s_map[vi] = v2b_map[vi] + vl*(hl+1);
	}
}

void StructuredMesh::assemble_pmat()
{
	Eigen::VectorX<FloatType> xrange = Eigen::VectorX<FloatType>::LinSpaced(hl+1, 0.0, 1.0);
	Eigen::VectorX<FloatType> zrange = Eigen::VectorX<FloatType>::LinSpaced(vl+1, 0.0, 1.0);

	int dof = 0;
	for (FloatType z : zrange) {
		for (FloatType x : xrange) {
			pmat(dof, 0) = x;
			pmat(dof, 1) = z;
			pmat_unit_box(dof, 0) = x;
			pmat_unit_box(dof, 1) = z;
			dof++;
		}
	}
}

inline void StructuredMesh::triangle_right_connectivity(int ci, int i, int j)
{
	if (_degree == 1) {
		cmat(ci, 0) = (hl+1)*i + j;
		cmat(ci, 1) = (hl+1)*i + j+1;
		cmat(ci, 2) = (hl+1)*(i+1) + j;
	
		cmat(ci+1, 0) = (hl+1)*(i+1) + j+1;
		cmat(ci+1, 1) = (hl+1)*(i+1) + j;
		cmat(ci+1, 2) = (hl+1)*i + (j+1);
	} else if (_degree == 2) {
		cmat(ci, 0) = (hl+1)*i + j;
		cmat(ci, 1) = (hl+1)*i + j+1;
		cmat(ci, 2) = (hl+1)*i + j+2;	
		cmat(ci, 3) = (hl+1)*(i+1) + j+1;
		cmat(ci, 4) = (hl+1)*(i+2) + j;
		cmat(ci, 5) = (hl+1)*(i+1) + j;
	
		cmat(ci+1, 0) = (hl+1)*(i+2) + j+2;
		cmat(ci+1, 1) = (hl+1)*(i+2) + j+1;
		cmat(ci+1, 2) = (hl+1)*(i+2) + j;	
		cmat(ci+1, 3) = (hl+1)*(i+1) + j+1;
		cmat(ci+1, 4) = (hl+1)*i + j+2;
		cmat(ci+1, 5) = (hl+1)*(i+1) + j+2;
	}
}

inline void StructuredMesh::triangle_left_connectivity(int ci, int i, int j)
{
	if (_degree == 1) {
		cmat(ci, 0) = (hl+1)*(i+1) + j;
		cmat(ci, 1) = (hl+1)*i + j;
		cmat(ci, 2) = (hl+1)*(i+1) + j+1;
	
		cmat(ci+1, 0) = (hl+1)*i + j+1;
		cmat(ci+1, 1) = (hl+1)*(i+1) + j+1;
		cmat(ci+1, 2) = (hl+1)*i + j;
	} else if (_degree == 2) {
		cmat(ci, 0) = (hl+1)*(i+2) + j;
		cmat(ci, 1) = (hl+1)*(i+1) + j;
		cmat(ci, 2) = (hl+1)*i + j;	
		cmat(ci, 3) = (hl+1)*(i+1) + j+1;
		cmat(ci, 4) = (hl+1)*(i+2) + j+2;
		cmat(ci, 5) = (hl+1)*(i+2) + j+1;
	
		cmat(ci+1, 0) = (hl+1)*i + j+2;
		cmat(ci+1, 1) = (hl+1)*(i+1) + j+2;
		cmat(ci+1, 2) = (hl+1)*(i+2) + j+2;	
		cmat(ci+1, 3) = (hl+1)*(i+1) + j+1;
		cmat(ci+1, 4) = (hl+1)*i + j;
		cmat(ci+1, 5) = (hl+1)*i + j+1;
	}
}

inline void StructuredMesh::quadrilateral_connectivity(int ci, int i, int j)
{
	if (_degree == 1) {
		cmat(ci, 0) = (hl+1)*i + j;
		cmat(ci, 1) = (hl+1)*i + j+1;
		cmat(ci, 2) = (hl+1)*(i+1) + j+1;
		cmat(ci, 3) = (hl+1)*(i+1) + j;
	} else if (_degree == 2) {
		cmat(ci, 0) = (hl+1)*i + j;
		cmat(ci, 1) = (hl+1)*i + j+1;
		cmat(ci, 2) = (hl+1)*i + j+2;
		cmat(ci, 3) = (hl+1)*(i+1) + j+2;
		cmat(ci, 4) = (hl+1)*(i+2) + j+2;
		cmat(ci, 5) = (hl+1)*(i+2) + j+1;
		cmat(ci, 6) = (hl+1)*(i+2) + j;
		cmat(ci, 7) = (hl+1)*(i+1) + j;
		cmat(ci, 8) = (hl+1)*(i+1) + j+1;
	}
}

void StructuredMesh::assemble_cmat()
{
	int ci = 0;
	for (int i = 0; i < vl; i+=_degree)
		for (int j = 0; j < hl; j+=_degree) {
			if (_cell_type == TRIANGLE_RIGHT) {
				triangle_right_connectivity(ci, i, j);
				ci += 2;
			} else if (_cell_type == TRIANGLE_LEFT) {
				triangle_left_connectivity(ci, i, j);
				ci += 2;
			} else if (_cell_type == QUADRILATERAL) {
				quadrilateral_connectivity(ci, i, j);
				ci += 1;
			}
		}	
}

void StructuredMesh::assemble_emat()
{
	emat = Eigen::MatrixXi(_nof_edges, _dofs_per_edge);

	std::set<std::tuple<int, int>> visited_edges;
	std::tuple<int, int> edge;

	int first_ind, last_ind, ei = 0;
	for (int ci = 0; ci < cmat.rows(); ci++) {
		Eigen::VectorXi cell = cmat.row(ci);
		for (int i = 0; i < _edges_per_cell; i++) {
			first_ind = i*_degree;
			if (i != _edges_per_cell-1)
				last_ind = (i+1)*_degree;
			else
				last_ind = 0;

			if (cell(first_ind) < cell(last_ind))
				edge = std::tuple<int, int>(cell(first_ind), cell(last_ind));
			else
				edge = std::tuple<int, int>(cell(last_ind), cell(first_ind));

			if (visited_edges.find(edge) == visited_edges.end()) {
				visited_edges.insert(edge);
				for (int j = 0; j < _dofs_per_edge-1; j++) {
					emat(ei, j) = cell(first_ind + j); 
				}
				emat(ei, _dofs_per_edge-1) = cell(last_ind);
				ei++;
			}
		}
	}
}

void StructuredMesh::assemble_cimat() 
{
	if (_cell_type == QUADRILATERAL)
		assemble_cimat_quadrilateral();
	else if (_cell_type == TRIANGLE_LEFT)
		assemble_cimat_triangle_left();
	else if (_cell_type == TRIANGLE_RIGHT)
		assemble_cimat_triangle_right();
}

void StructuredMesh::assemble_cimat_triangle_left() 
{
}

void StructuredMesh::assemble_cimat_triangle_right() 
{
}

void StructuredMesh::assemble_cimat_quadrilateral()
{
	int i, j, ci;

	cimat(0) = SOUTH_WEST_ID;
	cimat(nx-1) = SOUTH_EAST_ID;
	cimat(_nof_cells-nx) = NORTH_WEST_ID;
	cimat(_nof_cells-1) = NORTH_EAST_ID;

	i = nz-1;
	for (j = 1; j < nx-1; j++) {
		ci = i*nx + j;
		cimat(ci) = NORTH_ID;
	}

	j = 0;
	for (i = 1; i < nz-1; i++) {
		ci = i*nx + j;
		cimat(ci) = WEST_ID; 
	}

	j = nx-1;
	for (i = 1; i < nz-1; i++) {
		ci = i*nx + j;
		cimat(ci) = EAST_ID; 
	}

	i = 0;
	for (j = 1; j < nx-1; j++) {
		ci = i*(nx-1) + j;
		cimat(ci) = SOUTH_ID;
	}

	for (i = 1; i < nz-1; i++) {
		for (j = 1; j < nx-1; j++) {
			ci = i*nx + j;
			cimat(ci) = INTERIOR_ID; 
		}
	}

}

void StructuredMesh::assemble_dimat()
{
	FloatType domain_width = pmat(Eigen::last, 0) - pmat(0, 0);
	FloatType xtol = domain_width*1e-3/hl;
	FloatType ztol = 1e-3/vl;
	for (int i = 0; i < _nof_dofs; i++) {
		if (pmat(i, 1) > 1.0 - ztol && pmat(i, 0) < xtol)
			dimat(i) = NORTH_WEST_ID;
		else if (pmat(i, 1) < ztol && pmat(i, 0) < xtol)
			dimat(i) = SOUTH_WEST_ID;
		else if (pmat(i, 1) < ztol && pmat(i, 0) > 1.0 - xtol)
			dimat(i) = SOUTH_EAST_ID;
		else if (pmat(i, 1) > 1.0 - ztol && pmat(i, 0) > 1.0 - xtol)
			dimat(i) = NORTH_EAST_ID;
		else if (pmat(i, 1) > 1.0 - ztol)
			dimat(i) = NORTH_ID;
		else if (pmat(i, 1) < ztol)
			dimat(i) = SOUTH_ID;
		else if (pmat(i, 0) < xtol)
			dimat(i) = WEST_ID;
		else if (pmat(i, 0) > 1.0 - xtol)
			dimat(i) = EAST_ID;
		else
			dimat(i) = INTERIOR_ID;
	}
}

void StructuredMesh::assemble_eimat()
{
	int dpe = _dofs_per_edge;
	eimat = Eigen::MatrixXi::Zero(_nof_edges, dpe+1);
	FloatType domain_width = pmat(Eigen::last, 0) - pmat(0, 0);

	FloatType xtol = domain_width*1e-3/hl;
	FloatType ztol = 1e-3/vl;
	for (int ei = 0; ei < emat.rows(); ei++) {
		FloatType x1 = pmat(emat(ei, 0), 0);
		FloatType z1 = pmat(emat(ei, 0), 1);
		FloatType x2 = pmat(emat(ei, 1), 0);
		FloatType z2 = pmat(emat(ei, 1), 1);
		if (z1 > 1.0 - ztol && z2 > 1.0 - ztol)
			eimat(ei, 0) = NORTH_ID;
		else if (z1 < ztol && z2 < ztol)
			eimat(ei, 0) = SOUTH_ID;
		else if (x1 < xtol && x2 < xtol)
			eimat(ei, 0) = WEST_ID;
		else if (x1 > 1.0 - xtol && x2 > 1.0 - xtol)
			eimat(ei, 0) = EAST_ID;
		else
			eimat(ei, 0) = INTERIOR_ID;
	}

	std::set<std::tuple<int, int>> visited_edges;
	std::tuple<int, int> edge;

	int v1, v2;
	int ei = 0;
	for (int ci = 0; ci < cmat.rows(); ci++) {
		Eigen::VectorXi cell = cmat.row(ci);
		for (int i = 0; i < _edges_per_cell; i++) {
			v1 = cell(i*_degree);
			if (i != _edges_per_cell-1)
				v2 = cell((i+1)*_degree);
			else
				v2 = cell(0);

			if (v1 < v2)
				edge = std::tuple<int, int>(v1, v2);
			else
				edge = std::tuple<int, int>(v2, v1);

			if (visited_edges.find(edge) == visited_edges.end()) {
				eimat(ei, 1) = v1;
				eimat(ei, dpe) = v2;
				if (_degree == 2)
					eimat(ei, 2) = cell(i*_degree + 1);
				visited_edges.insert(edge);	
				ei++;
			}
		}
	}
}

void StructuredMesh::extrude_x(const Eigen::VectorX<FloatType> &xvec_p1)
{
	Eigen::VectorX<FloatType> x_vec;
	if (_degree == 1)
		x_vec = xvec_p1;
	else if (_degree == 2) {
		x_vec = Eigen::VectorX<FloatType>::Zero(hl+1);
		x_vec(Eigen::seq(0, Eigen::last, 2)) = xvec_p1;
		for (int i = 1; i < x_vec.size()-1; i+=2)
			x_vec(i) = 0.5*(x_vec(i-1) + x_vec(i+1)); 
	}

	int dof = 0;
	for (int i = 0; i < vl+1; i++) {
		for (int j = 0; j < hl+1; j++) {
			pmat(dof, 0) = x_vec(j);
			dof++;
		}
	}
	// recompute boundary normals and tangents
	compute_normals_and_tangents();
}

void StructuredMesh::compute_edge_normals_and_tangents()
{
	edge_tangents = Eigen::MatrixX<FloatType>::Zero(_nof_edges, 2);
	edge_normals = Eigen::MatrixX<FloatType>::Zero(_nof_edges, 2);

	std::vector<int> boundary_einds = extract_edge_inds(BOUNDARY_ID);
	for (int ei: boundary_einds) {
		int v1 = eimat(ei, 1);
		int v2 = eimat(ei, _dofs_per_edge);
		edge_tangents(ei, Eigen::all) = pmat.row(v2) - pmat.row(v1);
		edge_tangents(ei, Eigen::all) /= edge_tangents.row(ei).norm();
		edge_normals(ei, 0) = edge_tangents(ei, 1);
		edge_normals(ei, 1) = -edge_tangents(ei, 0);
	}
}

void StructuredMesh::compute_vertex_normals_and_tangents()
{
	vertex_tangents = dof_tangents(vertex_dof_inds, Eigen::all);
	vertex_normals = dof_normals(vertex_dof_inds, Eigen::all);
}

const int inline StructuredMesh::corner_dof_SW()
{
	return 0;
}

const int inline StructuredMesh::corner_dof_SE()
{
	return hl;
}

const int inline StructuredMesh::corner_dof_NE()
{
	return (hl+1)*(vl+1) - 1;
}

const int inline StructuredMesh::corner_dof_NW()
{
	return (hl+1)*vl;
}

bool inline StructuredMesh::is_corner_dof(int vi)
{
	return (
		vi == corner_dof_SW() || vi == corner_dof_SE() ||
		vi == corner_dof_NE() || vi == corner_dof_NW()
	);
}

void StructuredMesh::compute_dof_normals_and_tangents()
{
	dof_tangents = Eigen::MatrixX<FloatType>::Zero(_nof_dofs, 2);
	dof_normals = Eigen::MatrixX<FloatType>::Zero(_nof_dofs, 2);

	// compute normals/tangents for vertex dofs
	std::vector<int> boundary_einds = extract_edge_inds(BOUNDARY_ID);
	for (int i = 0; i < boundary_einds.size(); i++) {
		int ei = boundary_einds[i];
		int v1 = eimat(ei, 1);
		int v2 = eimat(ei, _dofs_per_edge);
		// exclude corners
		dof_tangents(v1, Eigen::all) += edge_tangents.row(ei);
		dof_normals(v1, Eigen::all) += edge_normals.row(ei);
		dof_tangents(v2, Eigen::all) += edge_tangents.row(ei);
		dof_normals(v2, Eigen::all) += edge_normals.row(ei);
	}

	// normalize vertex dofs normals/tangents
	std::vector<int> boundary_vertex_dofs = extract_vertex_dof_inds(BOUNDARY_ID);
	for (int di: boundary_vertex_dofs) {
		dof_tangents(di, Eigen::all) /= dof_tangents.row(di).norm();
		dof_normals(di, Eigen::all) /= dof_normals.row(di).norm();
	}

	// compute normals/tangents for edge dofs
	if (_degree == 2)
		for (int ei: boundary_einds) {
			int v1 = eimat(ei, 1);
			int ve = eimat(ei, 2);
			int v2 = eimat(ei, _dofs_per_edge);
			dof_tangents(ve, Eigen::all) = edge_tangents.row(ei);
			dof_normals(ve, Eigen::all) = edge_normals.row(ei);
		}
	
	// set normal/tangent at SW corner to edge normal/tangent
	dof_tangents(
		corner_dof_SW(), Eigen::all
	) = dof_tangents(corner_dof_SW()+1, Eigen::all);
	dof_normals(
		corner_dof_SW(), Eigen::all
	) = dof_normals(corner_dof_SW()+1, Eigen::all);

	// set normal/tangent at SE corner to edge normal/tangent
	dof_tangents(
		corner_dof_SE(), Eigen::all
	) = dof_tangents(corner_dof_SE()-1, Eigen::all);
	dof_normals(
		corner_dof_SE(), Eigen::all
	) = dof_normals(corner_dof_SE()-1, Eigen::all);

	// set normal/tangent at NE corner to edge normal/tangent
	dof_tangents(
		corner_dof_NE(), Eigen::all
	) = dof_tangents(corner_dof_NE()-1, Eigen::all);
	dof_normals(
		corner_dof_NE(), Eigen::all
	) = dof_normals(corner_dof_NE()-1, Eigen::all);

	// set normal/tangent at NW corner to edge normal/tangent
	dof_tangents(
		corner_dof_NW(), Eigen::all
	) = dof_tangents(corner_dof_NW()+1, Eigen::all);
	dof_normals(
		corner_dof_NW(), Eigen::all
	) = dof_normals(corner_dof_NW()+1, Eigen::all);
}

void StructuredMesh::compute_normals_and_tangents()
{
	compute_edge_normals_and_tangents();
	compute_dof_normals_and_tangents();
	compute_vertex_normals_and_tangents();
}

void StructuredMesh::extrude_x(FloatType x0, FloatType x1)
{
	Eigen::VectorX<FloatType> ones = Eigen::VectorX<FloatType>::Ones(pmat.rows());
	pmat(Eigen::all, 0) = (ones-pmat(Eigen::all, 0))*x0 + pmat(Eigen::all, 0)*x1;
	// recompute boundary normals and tangents
	compute_normals_and_tangents();
}

std::vector<int> StructuredMesh::extract_edge_inds(int id)
{
	std::vector<int> einds;
	for (int ei = 0; ei < eimat.rows(); ei++)
		if (eimat(ei, 0) & id)
			einds.push_back(ei);
	return einds;
}

std::vector<int> StructuredMesh::extract_cell_inds(int id)
{
	std::vector<int> cell_inds;
	for (int ci = 0; ci < cimat.rows(); ci++)
		if (cimat(ci) & id)
			cell_inds.push_back(ci);
	return cell_inds;
}

std::vector<int> StructuredMesh::extract_vertex_dof_inds(int id)
{
	std::vector<int> dinds;
	for (int vi: vertex_dof_inds)
		if (dimat(vi) & id)
			dinds.push_back(vi);
	return dinds;
}

std::vector<int> StructuredMesh::extract_dof_inds(int id)
{
	std::vector<int> dinds;
	for (int di = 0; di < dimat.rows(); di++)
		if (dimat(di) & id)
			dinds.push_back(di);
	return dinds;
}

/**
 * 
 * Determine interior node position by linearly interpolating between bottom z=b(x), and top z=s(x): h(x) = (1-x)*b(x) + x*s(x)
 * @param values bedrock vector b and surface vector s
 *
 */
void StructuredMesh::extrude_z(
	std::function<FloatType (FloatType)> zb_expr,
	std::function<FloatType (FloatType)> zs_expr
)
{
	std::vector<int> surf_vert_inds = extract_vertex_dof_inds(SURFACE_ID);	
	Eigen::VectorX<FloatType> xvec_p1 = pmat(surf_vert_inds, 0);
	Eigen::VectorX<FloatType> zb_vec_p1 = xvec_p1.unaryExpr(zb_expr);
	Eigen::VectorX<FloatType> zs_vec_p1 = xvec_p1.unaryExpr(zs_expr);
	extrude_z(zb_vec_p1, zs_vec_p1);
}

void StructuredMesh::extrude_z(const Eigen::VectorX<FloatType> &zs_vec_p1)
{	
	std::vector<int> bed_vert_inds = extract_vertex_dof_inds(BED_ID);
	Eigen::VectorX<FloatType> zb_vec_p1 = pmat(bed_vert_inds, 1);
	extrude_z(zb_vec_p1, zs_vec_p1);
}

/** 
 *
 * Determine interior node position by linearly interpolating between bottom z=b(x), and top z=s(x): h(x) = (1-x)*b(x) + x*s(x)
 * @param values bedrock b and surface s
 *
 */
void StructuredMesh::extrude_z(
	const Eigen::VectorX<FloatType> &zb_vec_p1, const Eigen::VectorX<FloatType> &zs_vec_p1
)
{
	Eigen::VectorX<FloatType> zb_vec_p2 = Eigen::VectorX<FloatType>::Zero(hl+1);
	Eigen::VectorX<FloatType> zs_vec_p2 = Eigen::VectorX<FloatType>::Zero(hl+1);

	// precalculate displacement of surface and bedrock
	int k = 0;
	for (int i = 0; i < hl+1; i+=_degree) {
		zb_vec_p2[i] = zb_vec_p1[k];
		zs_vec_p2[i] = zs_vec_p1[k];
		k++;
	}
	// set the correct z position of all cell corner nodes
	for (int i = 0; i < vl+1; i+=_degree) {
		for (int j = 0; j < hl+1; j+=_degree) {
			int vi = i*(hl+1) + j;
			FloatType z = pmat_unit_box(vi, 1);
			pmat(vi, 1) = (1.0-z)*zb_vec_p2(v2b_map(vi)) + z*zs_vec_p2(v2b_map(vi));
		}	
	}	

	if (_degree == 2) {
		// set the correct z position of all cell interior nodes
		for (int i = 1; i < vl; i+=_degree) {
			for (int j = 1; j < hl; j+=_degree) {
				int vi_NW = (i+1)*(hl+1) + j-1;
				int vi_NE = (i+1)*(hl+1) + j+1;
				int vi_SW = (i-1)*(hl+1) + j-1;
				int vi_SE = (i-1)*(hl+1) + j+1;
				int vi = i*(hl+1) + j;
				if (_cell_type == TRIANGLE_LEFT)
					pmat(vi, 1) = 0.5*(pmat(vi_NE, 1) + pmat(vi_SW, 1));
				else if (_cell_type == TRIANGLE_RIGHT)
					pmat(vi, 1) = 0.5*(pmat(vi_NW, 1) + pmat(vi_SE, 1));
				else if (_cell_type == QUADRILATERAL)
					pmat(vi, 1) = 0.25*(pmat(vi_NW, 1) + pmat(vi_NE, 1) + pmat(vi_SW, 1) + pmat(vi_SE, 1));
			}	
		}	

		// set the correct z position of all horizontal edge nodes
		for (int i = 0; i < vl+1; i+=_degree) {
			for (int j = 1; j < hl; j+=_degree) {
				int vi_W = i*(hl+1) + j-1;
				int vi_E = i*(hl+1) + j+1;
				int vi = i*(hl+1) + j;	
				pmat(vi, 1) = 0.5*(pmat(vi_W, 1) + pmat(vi_E, 1));
			}	
		}	

		// set the correct z position of all vertical edge nodes
		for (int i = 1; i < vl; i+=_degree) {
			for (int j = 0; j < hl+1; j+=_degree) {
				int vi_N = (i+1)*(hl+1) + j;
				int vi_S = (i-1)*(hl+1) + j;
				int vi = i*(hl+1) + j;
				pmat(vi, 1) = 0.5*(pmat(vi_N, 1) + pmat(vi_S, 1));
			}	
		}	
	}

	// recompute boundary normals and tangents
	compute_normals_and_tangents();
}

void StructuredMesh::compute_vertex_dof_inds() {
	int n_verts = (nx+1)*(nz+1);
	vertex_dof_inds.reserve(n_verts);

	for (int i = 0; i < vl+1; i+=_degree) {
		for (int j = 0; j < hl+1; j+=_degree) {
			vertex_dof_inds.push_back(i*(hl+1) + j);
		}
	}
}

int StructuredMesh::nof_cells()
{
	return _nof_cells;
}

int StructuredMesh::dofs_per_cell()
{
	return _dofs_per_cell;
}

int StructuredMesh::dofs_per_edge()
{
	return _dofs_per_edge;
}

int StructuredMesh::cell_type()
{
	return _cell_type;
}

int StructuredMesh::lhs_nnz()
{
	return dof_adj_mat.nonZeros();
}
