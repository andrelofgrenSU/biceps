/*
 * Copyright (C) 2025 André Löfgren
 *
 * This file is part of Biceps.
 *
 * Biceps is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Biceps is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Biceps. If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_spmat_addons.hpp>

#include <interval_mesh.hpp>
#include <Eigen/Sparse>

/**
 * @namespace FEM1D
 * @brief Namespace for 1D Finite Element Method (FEM) functions and operations.
 */
namespace FEM1D {
    /**
     * @brief Maps quadrature points and shape functions from reference to physical element.
     * 
     * @param[in] degree Polynomial degree of the basis functions.
     * @param[in] node_coords Coordinates of the element nodes.
     * @param[in] qpoints_r Quadrature points in the reference element.
     * @param[in] phi_r Shape function values in the reference element.
     * @param[in] grad_phi_r Shape function gradients in the reference element.
     * @param[out] detJ_r_ret Determinants of the Jacobian for mapping.
     * @param[out] qpoints_x_ret Quadrature points in the physical element.
     * @param[out] grad_phi_x_ret Shape function gradients in the physical element.
     */
    void map_to_reference_cell(
        int degree, 
        Eigen::MatrixXd &node_coords,
        Eigen::VectorXd &qpoints_r,
        Eigen::MatrixXd &phi_r,
        Eigen::MatrixXd &grad_phi_r,
        Eigen::VectorXd &detJ_r_ret,
        Eigen::MatrixXd &qpoints_x_ret,
        Eigen::MatrixXd &grad_phi_x_ret
    );

    /**
     * @brief Computes Lagrange basis functions and their gradients.
     * 
     * @param[in] degree Polynomial degree of the basis functions.
     * @param[in] qpoints_r Quadrature points in the reference element.
     * @param[out] phi_r_ret Shape function values at quadrature points.
     * @param[out] grad_phi_r_ret Shape function gradients at quadrature points.
     */
    void lagrange_basis(
        const int degree,
        Eigen::VectorXd &qpoints_r,
        Eigen::MatrixXd &phi_r_ret,
        Eigen::MatrixXd &grad_phi_r_ret
    );

    /**
     * @brief Returns the reference element points for a given degree.
     * 
     * @param[in] degree Polynomial degree of the basis functions.
     * @return Eigen::MatrixXd Matrix of reference element points.
     */
    Eigen::MatrixXd reference_element_points_rs(const int degree);

    /**
     * @brief Computes Gauss-Legendre quadrature points and weights.
     * 
     * @param[in] precision Number of quadrature points.
     * @param[out] points Computed quadrature points.
     * @param[out] weights Corresponding quadrature weights.
     */
    void gauss_legendre_quadrature(
        const int precision,
        Eigen::VectorXd &points,
        Eigen::VectorXd &weights
    );

    /**
     * @brief Assembles the mass matrix using Gaussian quadrature.
     * 
     * @param[in] mesh Interval mesh for the finite element method.
     * @param[in] gp Number of quadrature points.
     * @return Eigen::SparseMatrix<double> Sparse mass matrix.
     */
    Eigen::SparseMatrix<double> assemble_mass_matrix(IntervalMesh &mesh, int gp);

    /**
     * @brief Computes the inner product u^T * M * v.
     * 
     * @param[in] u_vec First vector.
     * @param[in] M Sparse mass matrix.
     * @param[in] v_vec Second vector.
     * @return double Inner product result.
     */
    double inner(
        const Eigen::VectorXd &u_vec,
        const Eigen::SparseMatrix<double> &M,
        const Eigen::VectorXd &v_vec
    );
};
