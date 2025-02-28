#pragma once
#include <string>
#include <XdmfReader.hpp>
#include <XdmfWriter.hpp>
#include <XdmfDomain.hpp>
#include <XdmfGrid.hpp>
#include <XdmfGeometry.hpp>
#include <XdmfGeometryType.hpp>
#include <XdmfTopology.hpp>
#include <XdmfTopologyType.hpp>
#include <XdmfGridCollection.hpp>
#include <XdmfGridCollectionType.hpp>
#include <XdmfTime.hpp>
#include <fem_function_2d.hpp>

class IOHandler {

	private:
		std::string csv_dir;
		std::string xdmf_dir;
		std::string xdmf_file_path;

		shared_ptr<XdmfUnstructuredGrid> xdmf_grid;
		shared_ptr<XdmfWriter> xdmf_writer;
		shared_ptr<XdmfGridCollection> xdmf_grid_collection;
		shared_ptr<XdmfDomain> xdmf_domain;
	
	public:
		explicit IOHandler(std::string output_dir);

		Eigen::VectorX<long double> csv_read_data(std::string filename);
		void csv_write_pmat(const Eigen::MatrixX<long double> &pmat);
		void csv_write_cmat(const Eigen::MatrixXi &cmat);
		void csv_write_data(const Eigen::VectorX<long double> data_vec);
		
		void xdmf_init_time_series(std::string filename);
		void xdmf_prepare_data_insertion(StructuredMesh &mesh, long double time);
		void xdmf_insert_data(
			FEMFunction2D &fem_function, std::string attr_name
		);
		void xdmf_insert_marker_data(StructuredMesh &mesh);
		void xdmf_finalize_data_insertion();
};
