#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <io_handler.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

IOHandler::IOHandler(std::string output_path)
{
	xdmf_dir = (boost::format{"%s/xdmf/"} %output_path).str();
	csv_dir = (boost::format{"%s/csv/"} %output_path).str();
}

Eigen::VectorX<long double> IOHandler::csv_read_data(std::string filename)
{
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << 
			boost::format("failed to open file \"%s\"") %filename
		<< std::endl;
		exit(-1);
	}

	std::string line_str;
	std::vector<long double> data_vec;
	while (std::getline(file, line_str)) {
		data_vec.push_back(atof(line_str.c_str()));
	}
	file.close();
	Eigen::VectorX<long double> vec_ret = Eigen::VectorX<long double>::Map(data_vec.data(), data_vec.size());
	return vec_ret;
}

//Eigen::VectorX<long double> IOHandler::csv_write_data(std::string filename)
//{
//	std::ifstream file(filename);
//	if (!file.is_open()) {
//		std::cerr << 
//			boost::format("failed to open file \"%s\"") %filename
//		<< std::endl;
//		exit(-1);
//	}
//
//	std::string line_str;
//	std::vector<long double> data_vec;
//	while (std::getline(file, line_str)) {
//		data_vec.push_back(atof(line_str.c_str()));
//	}
//	file.close();
//	Eigen::VectorX<long double> vec_ret = Eigen::VectorX<long double>::Map(data_vec.data(), data_vec.size());
//	return vec_ret;
//}

void IOHandler::xdmf_init_time_series(std::string filename)
{
	boost::filesystem::create_directories(xdmf_dir);
	xdmf_file_path = (boost::format("%s/%s") %xdmf_dir %filename).str();
	shared_ptr<XdmfDomain> domain = XdmfDomain::New();
	shared_ptr<XdmfGridCollection> grid_collection = XdmfGridCollection::New();
	grid_collection->setType(XdmfGridCollectionType::Temporal());
	grid_collection->setName("TimeSeries");
	domain->insert(grid_collection);

	domain->getGridCollection("TimeSeries");
	shared_ptr<XdmfWriter> writer = XdmfWriter::New(xdmf_file_path);
	domain->accept(writer);
}

void IOHandler::xdmf_prepare_data_insertion(StructuredMesh &mesh, long double time)
{
	shared_ptr<XdmfReader> reader = XdmfReader::New();
	xdmf_domain = shared_dynamic_cast<XdmfDomain>(reader->read(xdmf_file_path));
	xdmf_grid_collection = xdmf_domain->getGridCollection("TimeSeries");
	shared_ptr<XdmfTime> xdmf_time = XdmfTime::New(time);
	xdmf_grid = XdmfUnstructuredGrid::New();
	shared_ptr<XdmfTopology> topology = XdmfTopology::New();
	shared_ptr<XdmfGeometry> geometry = XdmfGeometry::New();
	xdmf_writer = XdmfWriter::New(xdmf_file_path);

	xdmf_grid->setName("Mesh");
	xdmf_grid->setTime(xdmf_time);

	// The geometry type specifies the mesh node coordinates
	shared_ptr<const XdmfGeometryType> geo_type = XdmfGeometryType::XY();
	geometry->setType(geo_type);

	// The topology type specifies the connectivity between nodes (triangles, quads, etc.)
	shared_ptr<const XdmfTopologyType> topo_type;
	if (
		mesh.cell_type() == TRIANGLE_LEFT ||
		mesh.cell_type() == TRIANGLE_RIGHT
	)
		topo_type = XdmfTopologyType::Triangle();
	else if (mesh.cell_type() == QUADRILATERAL)
		topo_type = XdmfTopologyType::Quadrilateral();
	else
		;//throw exception
	
	// Eigen matrices are stored column major; need to transpose
	// in order to get the correct memory layout
	Eigen::MatrixX<long double> pmat = mesh.pmat.transpose();
	Eigen::MatrixXi cmat = mesh.cmat.transpose();

	geometry->insert(0, &pmat(0, 0), pmat.size());
	topology->insert(0, &cmat(0, 0), cmat.size());

	topology->setType(topo_type);

	xdmf_grid->setGeometry(geometry);
	xdmf_grid->setTopology(topology);
}

void IOHandler::xdmf_insert_data(FEMFunction2D &fem_function, std::string attr_name)
{
	// insert fem values
	std::vector<int> vinds = fem_function.mesh.extract_vertex_dof_inds(DOMAIN_ID);
	Eigen::VectorX<long double> vals = fem_function.vals(vinds);
	shared_ptr<XdmfAttribute> val_attr = XdmfAttribute::New();
	val_attr->setName(attr_name);
	val_attr->setCenter(XdmfAttributeCenter::Node());
	val_attr->setType(XdmfAttributeType::Scalar());
	val_attr->initialize(XdmfArrayType::Float64());
	val_attr->insert(0, &vals(0), vals.size());
	xdmf_grid->insert(val_attr);
}

void IOHandler::xdmf_insert_marker_data(StructuredMesh &mesh)
{
	std::vector<int> vinds = mesh.extract_vertex_dof_inds(DOMAIN_ID);
	Eigen::VectorXi vm_vec = mesh.dimat(vinds);
	shared_ptr<XdmfAttribute> vm_attr = XdmfAttribute::New();
	vm_attr->setName("Vertex marker");
	vm_attr->setCenter(XdmfAttributeCenter::Node());
	vm_attr->setType(XdmfAttributeType::Scalar());
	vm_attr->initialize(XdmfArrayType::Int32());
	vm_attr->insert(0, &vm_vec(0), vm_vec.size());
	xdmf_grid->insert(vm_attr);
}

void IOHandler::xdmf_finalize_data_insertion()
{
	xdmf_grid_collection->insert(xdmf_grid);
	xdmf_domain->accept(xdmf_writer);
}
