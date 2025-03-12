cmake_minimum_required(VERSION 3.10)
project(biceps VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BUILD_SHARED_LIBS ON)
set(CMAKE_BUILD_TYPE Release)
set(USE_LONG_DOUBLE OFF)
set(DOXYGEN_HTML ON)
set(DOXYGEN_LATEX ON)


find_package(Doxygen QUIET)
if(ENABLE_DOCS AND DOXYGEN_FOUND)
    message(STATUS "Doxygen found, enabling documentation target.")

    set(DOXYGEN_OUTPUT_DIR ${CMAKE_BINARY_DIR}/docs)
    set(DOXYGEN_CONFIG_FILE ${CMAKE_BINARY_DIR}/Doxyfile)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${DOXYGEN_CONFIG_FILE} @ONLY)

    add_custom_target(
        doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONFIG_FILE}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Generate documentation using Doxygen"
        VERBATIM
    )
else()
    if(ENABLE_DOCS)
        message(WARNING "Doxygen not found, but ENABLE_DOCS is ON. Documentation will not be built.")
    else()
        message(STATUS "Documentation build is disabled.")
    endif()
endif()

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(Python3 COMPONENTS Development NumPy REQUIRED)
find_package(eigenpy REQUIRED)
find_package(Boost 1.74 COMPONENTS python unit_test_framework REQUIRED)

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
        include/pstokes_fem.hpp
        include/poisson_fem.hpp
        include/free_surface_fem.hpp
    DESTINATION
        include
)
