#
# Copyright (C) 2025 André Löfgren
#
# This file is part of Biceps.
#
# Biceps is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Biceps is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Biceps. If not, see <https://www.gnu.org/licenses/>.
#
include_directories(
    ${CMAKE_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
)

add_executable(mesh_test_suite ${CMAKE_CURRENT_SOURCE_DIR}/mesh_test_suite.cpp)
target_link_libraries(mesh_test_suite biceps Eigen3::Eigen ${Boost_LIBRARIES})

add_executable(fem_2d_test_suite ${CMAKE_CURRENT_SOURCE_DIR}/fem_2d_test_suite.cpp)
target_link_libraries(fem_2d_test_suite Eigen3::Eigen biceps ${Boost_LIBRARIES} ${UMFPACK_LIBRARIES})

add_executable(poisson_test_suite ${CMAKE_CURRENT_SOURCE_DIR}/poisson_test_suite.cpp)
target_link_libraries(poisson_test_suite Eigen3::Eigen biceps ${Boost_LIBRARIES} ${UMFPACK_LIBRARIES})

add_executable(stokes_test_suite ${CMAKE_CURRENT_SOURCE_DIR}/stokes_test_suite.cpp)
target_link_libraries(stokes_test_suite Eigen3::Eigen biceps ${Boost_LIBRARIES} ${UMFPACK_LIBRARIES})

add_custom_target(
    test_mesh
    DEPENDS
    biceps
    COMMAND
    ./mesh_test_suite 
)

add_custom_target(
    test_fem_2d
    DEPENDS
    biceps
    COMMAND
    ./fem_2d_test_suite
)

add_custom_target(
    test_poisson
    DEPENDS
    biceps
    COMMAND
    ./poisson_test_suite 
)

add_custom_target(
    test_stokes
    DEPENDS
    biceps
    COMMAND
    ./stokes_test_suite
)

add_custom_target(
    tests
    DEPENDS
    test_mesh
    test_fem_2d
    test_poisson
    test_stokes
    COMMAND
    ./mesh_test_suite && ./fem2d_test_suite && ./poisson_test_suite.cpp && ./stokes_test_suite
)
