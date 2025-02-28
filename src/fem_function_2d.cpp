#include <iostream>
// #include <Eigen/UmfPackSupport>
#include <fem_2d.hpp>

FEMFunction2D::FEMFunction2D(StructuredMesh &mesh) : mesh(mesh)
{
	vals = Eigen::VectorX<FloatType>::Zero(mesh.nof_dofs());
}

void FEMFunction2D::assign(std::function<FloatType(FloatType, FloatType)> func)
{
	for (int i = 0; i < mesh.nof_dofs(); i++) {
		FloatType x = mesh.pmat(i, 0);
		FloatType z = mesh.pmat(i, 1);
		vals(i) = func(x, z);
	}
}

void FEMFunction2D::assign(const Eigen::VectorX<FloatType> &vec)
{
	if (vec.size() != mesh.nof_dofs()) {
		std::cerr << "error:fem_function_2d:assign(): vec.size() != mesh.nof_dofs()" << std::endl;
		exit(-1);
	}
	for (int dof = 0; dof < vec.size(); dof++) {
		vals(dof) = vec(dof);
	}
}

FEMFunction2D FEMFunction2D::diff_x_interp()
{
	Eigen::MatrixX<FloatType> node_coords, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> detJ_rs;
	Eigen::VectorXi element;
	FEMFunction2D u_diff(mesh);

	qpoints_rs = FEM2D::reference_element_points_rs(mesh.cell_type(), mesh.degree());
	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), mesh.dofs_per_cell());

	FEM2D::lagrange_basis(
		mesh.degree(),
		mesh.cell_type(),
		qpoints_rs,
		phi_rs, dphi_rs
	);

	Eigen::VectorX<FloatType> diff_vals = Eigen::VectorX<FloatType>::Zero(mesh.nof_dofs());
	Eigen::VectorXi duplicate_dof_count = Eigen::VectorXi::Zero(mesh.nof_dofs());
	for (int ci = 0; ci < mesh.cmat.rows(); ci++) {
		element = mesh.cmat.row(ci);
		node_coords = mesh.pmat(element, Eigen::all);
		FEM2D::map_to_reference_cell(
			mesh.degree(), mesh.cell_type(), node_coords, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);
		for (int q = 0; q < qpoints_rs.rows(); q++) {
			int w = 2*q;
			duplicate_dof_count(element(q)) += 1;
			diff_vals(element(q)) += vals(element).dot(dphi_xz(w, Eigen::all));
		}
	}
	for (int dof = 0; dof < mesh.nof_dofs(); dof++) {
		u_diff.vals(dof) = diff_vals(dof)/duplicate_dof_count(dof);
	}
	return u_diff;
}

FEMFunction2D FEMFunction2D::diff_z_interp()
{
	Eigen::MatrixX<FloatType> node_coords, qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> detJ_rs;
	Eigen::VectorXi element;
	FEMFunction2D u_diff(mesh);

	qpoints_rs = FEM2D::reference_element_points_rs(mesh.cell_type(), mesh.degree());
	detJ_rs = Eigen::VectorX<FloatType>::Zero(qpoints_rs.rows());
	qpoints_xz = Eigen::MatrixX<FloatType>::Zero(qpoints_rs.rows(), 2);
	dphi_xz = Eigen::MatrixX<FloatType>::Zero(2*qpoints_rs.rows(), mesh.dofs_per_cell());

	FEM2D::lagrange_basis(
		mesh.degree(),
		mesh.cell_type(),
		qpoints_rs,
		phi_rs, dphi_rs
	);

	Eigen::VectorX<FloatType> diff_vals = Eigen::VectorX<FloatType>::Zero(mesh.nof_dofs());
	Eigen::VectorXi duplicate_dof_count = Eigen::VectorXi::Zero(mesh.nof_dofs());
	for (int ci = 0; ci < mesh.cmat.rows(); ci++) {
		element = mesh.cmat.row(ci);
		node_coords = mesh.pmat(element, Eigen::all);
		FEM2D::map_to_reference_cell(
			mesh.degree(), mesh.cell_type(), node_coords, qpoints_rs,
			phi_rs, dphi_rs, detJ_rs, qpoints_xz, dphi_xz
		);
		for (int q = 0; q < qpoints_rs.rows(); q++) {
			int w = 2*q;
			duplicate_dof_count(element(q)) += 1;
			diff_vals(element(q)) += vals(element).dot(dphi_xz(w+1, Eigen::all));
		}
	}
	for (int dof = 0; dof < mesh.nof_dofs(); dof++) {
		u_diff.vals(dof) = diff_vals(dof)/duplicate_dof_count(dof);
	}
	return u_diff;
}

FEMFunction2D FEMFunction2D::diff_x_proj(std::string sp_solver_name)
{
	FEMFunction2D u_diff(mesh);

	int deg = mesh.degree();
	Eigen::SparseMatrix<FloatType> M = FEM2D::assemble_mass_matrix(mesh, 2*deg);
	Eigen::SparseMatrix<FloatType> Ax = FEM2D::assemble_stiffness_x_matrix(mesh, deg+1);
	Eigen::VectorX<FloatType> rhs_vec = Ax*vals;

	if (sp_solver_name.compare("sparselu") == 0) {
		Eigen::SparseLU<Eigen::SparseMatrix<FloatType>> sp_solver;
		sp_solver.analyzePattern(M);
		sp_solver.factorize(M);
		u_diff.vals = sp_solver.solve(rhs_vec);
	} 
	// else if (sp_solver_name.compare("umfpack") == 0) {
	// 	Eigen::UmfPackLU<Eigen::SparseMatrix<FloatType>> sp_solver;
	// 	sp_solver.analyzePattern(M);
	// 	sp_solver.factorize(M);
	// 	u_diff.vals = sp_solver.solve(rhs_vec);
	// }

	return u_diff;
}

FEMFunction2D FEMFunction2D::diff_z_proj(std::string sp_solver_name)
{
	FEMFunction2D u_diff(mesh);

	int deg = mesh.degree();
	Eigen::SparseMatrix<FloatType> M = FEM2D::assemble_mass_matrix(mesh, 2*deg);
	Eigen::SparseMatrix<FloatType> Az = FEM2D::assemble_stiffness_z_matrix(mesh, deg+1);
	Eigen::VectorX<FloatType> rhs_vec = Az*vals;

	if (sp_solver_name.compare("sparselu") == 0) {
		Eigen::SparseLU<Eigen::SparseMatrix<FloatType>> sp_solver;
		sp_solver.analyzePattern(M);
		sp_solver.factorize(M);
		u_diff.vals = sp_solver.solve(rhs_vec);
	} 
	// else if (sp_solver_name.compare("umfpack") == 0) {
	// 	Eigen::UmfPackLU<Eigen::SparseMatrix<FloatType>> sp_solver;
	// 	sp_solver.analyzePattern(M);
	// 	sp_solver.factorize(M);
	// 	u_diff.vals = sp_solver.solve(rhs_vec);
	// }

	return u_diff;
}

Eigen::VectorX<FloatType> FEMFunction2D::eval_cell(int ci)
{
	return vals(mesh.cmat.row(ci), Eigen::all);	
}

Eigen::VectorX<FloatType> FEMFunction2D::eval_edge(int ei)
{
	return vals(mesh.emat.row(ei), Eigen::all);	
}

Eigen::MatrixX<FloatType> FEMFunction2D::extract_vertex_subvec(int domain_id)
{
	std::vector<int> subdomain_inds = mesh.extract_vertex_dof_inds(domain_id);
	Eigen::MatrixX<FloatType> vector_xzf(subdomain_inds.size(), 3);
	vector_xzf(Eigen::all, Eigen::seq(0, 1)) = mesh.pmat(subdomain_inds, Eigen::all);
	vector_xzf(Eigen::all, 2) = vals(subdomain_inds);
	return vector_xzf;
}

Eigen::MatrixX<FloatType> FEMFunction2D::extract_dof_subvec(int domain_id)
{
	std::vector<int> subdomain_inds = mesh.extract_dof_inds(domain_id);
	Eigen::MatrixX<FloatType> vector_xzf(subdomain_inds.size(), 3);
	vector_xzf(Eigen::all, Eigen::seq(0, 1)) = mesh.pmat(subdomain_inds, Eigen::all);
	vector_xzf(Eigen::all, 2) = vals(subdomain_inds);
	return vector_xzf;
}

// TODO: implement this!
FloatType FEMFunction2D::operator()(FloatType x, FloatType z)
{
	Eigen::MatrixX<FloatType> qpoints_rs, qpoints_xz,
		phi_rs, dphi_rs, dphi_xz;
	Eigen::VectorX<FloatType> qweights, detJ_rs;
	Eigen::VectorXi element_v;
	return 0.0;
}

FEMFunction2D FEMFunction2D::operator+(const FEMFunction2D &f)
{
	FEMFunction2D func(mesh);
	func.vals = this->vals + f.vals;
	return func;
}

FEMFunction2D FEMFunction2D::operator-(const FEMFunction2D &f)
{
	FEMFunction2D func(mesh);
	func.vals = this->vals - f.vals;
	return func;
}

FEMFunction2D FEMFunction2D::operator*(const FEMFunction2D &f)
{
	FEMFunction2D func(mesh);
	func.vals = this->vals * f.vals;
	return func;
}

FEMFunction2D FEMFunction2D::operator+(FloatType val)
{
	FEMFunction2D func(mesh);
	func.vals = this->vals.array() + val;
	return func;
}

FEMFunction2D FEMFunction2D::operator-(FloatType val)
{
	FEMFunction2D func(mesh);
	func.vals = this->vals.array() - val;
	return func;
}

FEMFunction2D FEMFunction2D::operator*(FloatType val)
{
	FEMFunction2D func(mesh);
	func.vals = this->vals * val;
	return func;
}
