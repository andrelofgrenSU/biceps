cmake_minimum_required(VERSION 3.10)
project(biceps VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(BUILD_SHARED_LIBS False)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")

set(USE_LONG_DOUBLE False)

find_package(UMFPACK)
find_package(HDF5 REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(XDMF REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(Boost 1.78 COMPONENTS unit_test_framework timer filesystem REQUIRED)

add_subdirectory(src)
add_subdirectory(unit_tests)

install(
	FILES 
		${CMAKE_BINARY_DIR}/float_type.hpp
		include/eigen_spmat_addons.hpp
		include/enums.hpp
		include/logger.hpp
		include/interval_mesh.hpp
		include/structured_mesh.hpp
		include/fem_1d.hpp
		include/fem_2d.hpp
		include/fem_function_1d.hpp
		include/fem_function_2d.hpp
		include/io_handler.hpp
		include/pstokes_fem.hpp
		include/poisson_fem.hpp
		include/free_surface_fem.hpp
		include/time_integrator.hpp
	DESTINATION
		include
)
