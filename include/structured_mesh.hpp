#pragma once

#include <float_type.hpp>
#include <vector>
#include <Eigen/Dense>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_spmat_addons.hpp>
#include <Eigen/Sparse>

enum CELL_TYPE {
	QUADRILATERAL,
	TRIANGLE_LEFT,
	TRIANGLE_RIGHT
};

enum DOMAIN_IDS {
	EMPTY_ID = 0,
	INTERIOR_ID = 1,
	NORTH_ID = 2,
	NORTH_WEST_ID = 4,
	WEST_ID = 8,
	SOUTH_WEST_ID = 16,
	SOUTH_ID = 32,
	SOUTH_EAST_ID = 64,
	EAST_ID = 128,
	NORTH_EAST_ID = 256,
	SURFACE_ID = NORTH_WEST_ID | NORTH_ID | NORTH_EAST_ID,
	BED_ID = SOUTH_WEST_ID | SOUTH_ID | SOUTH_EAST_ID,
	BOUNDARY_ID = NORTH_ID | NORTH_WEST_ID | WEST_ID | SOUTH_WEST_ID |
		SOUTH_ID | SOUTH_EAST_ID | EAST_ID | NORTH_EAST_ID,
	DOMAIN_ID = BOUNDARY_ID | INTERIOR_ID
};

class MeshGraph {
	private:

	protected:
		int _degree;
		int _nof_dofs;
		int _nof_verts;
		Eigen::SparseMatrix<int> cell_adj_mat;
		Eigen::SparseMatrix<int> dof_adj_mat;

		void assemble_dof_adj_mat();

	public:
		int nx, nz;
		int hl, vl;
		Eigen::VectorXi d2v;
		Eigen::VectorXi v2d;

		MeshGraph() {};
		MeshGraph(int nx, int nz, int degree);

		int vertex_degree(int vi);
		int degree();
		int nof_verts();
		int nof_dofs();
};

// TODO: implement a generic mesh class and make StructuredMesh into a subclass
class StructuredMesh : public MeshGraph {

	private:
		int _nof_cells;
		int _nof_edges;
		int _dofs_per_cell;
		int _edges_per_cell;
		int _dofs_per_edge;
		int _cell_type;

		void assemble_pmat();
		void assemble_cmat();
		void assemble_emat();
		void assemble_cimat();
		void assemble_eimat();
		void assemble_dimat();

		void assemble_cimat_quadrilateral();
		void assemble_cimat_triangle_left();
		void assemble_cimat_triangle_right();

		void compute_vertex_dof_inds();
		void compute_edge_normals_and_tangents();
		void compute_vertex_normals_and_tangents();
		void compute_dof_normals_and_tangents();
		void compute_normals_and_tangents();

		void triangle_left_connectivity(int ci, int i, int j);
		void triangle_right_connectivity(int ci, int i, int j);
		void quadrilateral_connectivity(int ci, int i, int j);

		const int corner_dof_SW();
		const int corner_dof_SE();
		const int corner_dof_NE();
		const int corner_dof_NW();
		bool is_corner_dof(int vi);

		std::vector<int> vertex_dof_inds;

	public:
		Eigen::MatrixX<FloatType> pmat_unit_box;
		Eigen::MatrixX<FloatType> pmat;
		Eigen::MatrixXi cmat;
		Eigen::MatrixXi emat;
		Eigen::MatrixXi eimat;
		Eigen::VectorXi dimat;
		Eigen::VectorXi cimat;

		Eigen::MatrixXi edge_list;
		Eigen::MatrixX<FloatType> edge_tangents;
		Eigen::MatrixX<FloatType> edge_normals;
		Eigen::MatrixX<FloatType> vertex_tangents;
		Eigen::MatrixX<FloatType> vertex_normals;
		Eigen::MatrixX<FloatType> dof_tangents;
		Eigen::MatrixX<FloatType> dof_normals;
		Eigen::VectorXi v2s_map;
		Eigen::VectorXi v2b_map;

		StructuredMesh(int nx, int nz, int degree, int cell_type);

		void extrude_x(const Eigen::VectorX<FloatType> &xvec_p1);
		void extrude_x(FloatType x0, FloatType x1);
		void extrude_z(
			std::function<FloatType (FloatType)> zb_expr,
			std::function<FloatType (FloatType)> zs_expr
		);
		void extrude_z(const Eigen::VectorX<FloatType> &zb_vec_p1, const Eigen::VectorX<FloatType> &zs_vec_p1);
		void extrude_z(const Eigen::VectorX<FloatType> &zs_vec_p1);
		std::vector<int> extract_cell_inds(int id);
		std::vector<int> extract_edge_inds(int id);
		std::vector<int> extract_dof_inds(int id);
		std::vector<int> extract_vertex_dof_inds(int id);
		int nof_cells();
		int nof_edges();
		int dofs_per_cell();
		int dofs_per_edge();
		int cell_type();
		int lhs_nnz();
};
