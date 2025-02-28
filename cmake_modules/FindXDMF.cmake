include(FindPackageHandleStandardArgs)

find_path(
	XDMF_INCLUDE
	NAMES
	Xdmf.hpp
	PATHS
	$ENV{XDMF_DIR}
	${INCLUDE_INSTALL_DIR}
)

find_library(
	XDMF_LIBRARY
	Xdmf
	PATHS
	${INCLUDE_INSTALL_DIR}
	$ENV{XDMF_DIR}
)

find_library(
	XDMF_CORE_LIBRARY
	XdmfCore
	PATHS
	${INCLUDE_INSTALL_DIR}
	$ENV{XDMF_DIR}
)

find_package_handle_standard_args(
	XDMF
	DEFAULT_MSG
	XDMF_INCLUDE
	XDMF_LIBRARY
	XDMF_CORE_LIBRARY
)

mark_as_advanced(XDMF_LIBRARY XDMF_CORE_LIBRARY)
