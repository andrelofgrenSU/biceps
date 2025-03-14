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
#include <logger.hpp>
#include <fem_function_2d.hpp>

/**
 * @class pStokesProblem
 * @brief A class to solve the 2D Stokes problem using the finite element method.
 * 
 * This class provides functionality for assembling the system of equations, 
 * applying boundary conditions, solving linear and nonlinear systems, 
 * and storing the results for velocity and pressure fields.
 */
class pStokesProblem {

private:
    int nv_dofs; ///< Number of velocity degrees of freedom
    int np_dofs; ///< Number of pressure degrees of freedom
    int n_dofs;  ///< Total number of degrees of freedom
    Logger logger; ///< Logger object for logging messages

public:
    std::vector<Eigen::Triplet<FloatType>> lhs_coeffs; ///< Left-hand side coefficients (stiffness matrix)
    int nof_pushed_elements = 0; ///< Counter for the number of elements pushed during assembly

    StructuredMesh &u_mesh; ///< Mesh for velocity field (horizontal and vertical)
    StructuredMesh &p_mesh; ///< Mesh for pressure field

    Eigen::VectorXi ux_v2d; ///< Velocity vector (horizontal)
    Eigen::VectorXi uz_v2d; ///< Velocity vector (vertical)
    Eigen::VectorXi u_v2d; ///< Combined velocity vector
    Eigen::VectorXi p_v2d; ///< Pressure vector

    Eigen::VectorXi ux_d2v; ///< Velocity degrees of freedom mapping (horizontal)
    Eigen::VectorXi uz_d2v; ///< Velocity degrees of freedom mapping (vertical)
    Eigen::VectorXi p_d2v; ///< Pressure degrees of freedom mapping
    Eigen::VectorXi w_d2v; ///< Mapping for the solution vector (unknowns)

    Eigen::SparseMatrix<FloatType> lhs_mat; ///< Left-hand side system matrix
    Eigen::VectorX<FloatType> rhs_vec; ///< Right-hand side force vector
    Eigen::VectorX<FloatType> w_vec; ///< Solution vector

    /**
     * @brief Constructor for initializing the pStokesProblem with meshes for velocity and pressure.
     * 
     * @param[in] u_mesh Reference to the mesh for velocity field (both horizontal and vertical).
     * @param[in] p_mesh Reference to the mesh for pressure field.
     */
    pStokesProblem(StructuredMesh &u_mesh, StructuredMesh &p_mesh);

    /**
     * @brief Assembles the stress block of the left-hand side matrix.
     * 
     * @param[in] A The rate factor.
     * @param[in] n_i The flow law exponent.
     * @param[in] eps_reg_2 The regularization parameter (avoids infinite viscosity).
     * @param[in] gauss_precision The precision for Gaussian quadrature used in integration.
     */
    void assemble_stress_block(
        FloatType A, FloatType n_i, FloatType eps_reg_2, int gauss_precision
    );

    /**
     * @brief Assembles the incompressibility block of the left-hand side matrix.
     * 
     * @param[in] gauss_precision The precision for Gaussian quadrature used in integration.
     */
    void assemble_incomp_block(int gauss_precision);

    /**
     * @brief Assembles the FSSA vertical block of the left-hand side matrix.
     * 
     * @param[in] stab_param The stabilization parameter for the FSSA method.
     * @param[in] force_x The force function in the x-direction (as a function of coordinates).
     * @param[in] force_z The force function in the z-direction (as a function of coordinates).
     * @param[in] gauss_precision The precision for Gaussian quadrature used in integration.
     */
    void assemble_fssa_vertical_block(
        FloatType stab_param,
        std::function<double(double, double)> force_x,
        std::function<double(double, double)> force_z,
        int gauss_precision
    );

    /**
     * @brief Assembles the FSSA normal block of the left-hand side matrix.
     * 
     * @param[in] stab_param The stabilization parameter for the FSSA method.
     * @param[in] force_z The force function in the z-direction (as a function of coordinates).
     * @param[in] gauss_precision The precision for Gaussian quadrature used in integration.
     */
    void assemble_fssa_normal_block(
        FloatType stab_param,
        std::function<double(double, double)> force_z,
        int gauss_precision
    );

    /**
     * @brief Commits the assembled left-hand side matrix to the system.
     */
    void commit_lhs_mat(); 

    /**
     * @brief Assembles the right-hand side vector based on external forces.
     * 
     * @param[in] force_x The external force function in the x-direction.
     * @param[in] force_z The external force function in the z-direction.
     * @param[in] gauss_precision The precision for Gaussian quadrature used in integration.
     */
    void assemble_rhs_vec(
        std::function<FloatType(FloatType, FloatType)> force_x,
        std::function<FloatType(FloatType, FloatType)> force_z,
        int gauss_precision
    );

    /**
     * @brief Assembles the FSSA vertical right-hand side vector with external forces.
     * 
     * @param[in] stab_param The stabilization parameter for the FSSA method.
     * @param[in] force_x The external force function in the x-direction.
     * @param[in] force_z The external force function in the z-direction.
     * @param[in] accum The accumulation function used in the FSSA method.
     * @param[in] gauss_precision The precision for Gaussian quadrature used in integration.
     */
    void assemble_fssa_vertical_rhs(
        FloatType stab_param,
        std::function<FloatType(FloatType, FloatType)> force_x,
        std::function<FloatType(FloatType, FloatType)> force_z,
        std::function<FloatType(FloatType, FloatType)> accum,
        int gauss_precision
    );

    /**
     * @brief Applies zero Dirichlet boundary conditions on the velocity components.
     * 
     * @param[in] ux_boundary_id The ID for the boundary where horizontal velocity is zero.
     * @param[in] uz_boundary_id The ID for the boundary where vertical velocity is zero.
     */
    void apply_zero_dirichlet_bc(int ux_boundary_id, int uz_boundary_id);

    /**
     * @brief Applies Dirichlet boundary conditions for a specific velocity component.
     * 
     * @param[in] boundary_part The ID for the boundary part.
     * @param[in] velocity_component The velocity component to which the boundary condition is applied.
     * @param[in] u_func The function representing the Dirichlet boundary condition.
     */
    void apply_dirichlet_bc(
        const int boundary_part,
        const int velocity_component,
        std::function<FloatType(FloatType, FloatType)> u_func
    );

    /**
     * @brief Prunes small coefficients from the left-hand side matrix to improve sparsity.
     * 
     * @param[in] threshold The coefficient magnitude below which entries are pruned.
     */
    void prune_lhs(FloatType threshold);

    /**
     * @brief Solves the linear system of equations.
     */
    void solve_linear_system();

    /**
     * @brief Solves the nonlinear system using the Picard method (iterative approach).
     * 
     * @param[in] A The matrix coefficient for stress.
     * @param[in] n_i The index for the stress term.
     * @param[in] eps_reg_2 The regularization parameter for the stress block.
     * @param[in] fssa_version The version of the FSSA method to use.
     * @param[in] fssa_param The FSSA parameter to use.
     * @param[in] force_x The external force function in the x-direction.
     * @param[in] force_z The external force function in the z-direction.
     * @param[in] ux_boundary_id The boundary ID for the horizontal velocity.
     * @param[in] uz_boundary_id The boundary ID for the vertical velocity.
     * @param[in] max_iter The maximum number of iterations allowed for the Picard method.
     * @param[in] stol The solution tolerance for the Picard method.
     * @param[in] gauss_precision The precision for Gaussian quadrature used in integration.
     */
    void solve_nonlinear_system_picard(
        FloatType A, FloatType n_i, FloatType eps_reg_2,
        int fssa_version, FloatType fssa_param,
        std::function<FloatType(FloatType, FloatType)> force_x,
        std::function<FloatType(FloatType, FloatType)> force_z,
        int ux_boundary_id, int uz_boundary_id,
        int max_iter, FloatType stol, int gauss_precision
    );

    /**
     * @brief Retrieves the horizontal velocity as a FEM function.
     * 
     * @return A FEMFunction2D object representing the horizontal velocity field.
     */
    FEMFunction2D velocity_x();

    /**
     * @brief Retrieves the veloctiy velocity as a FEM function.
     * @return FEMFunction2D representing the solution.
     */
    FEMFunction2D velocity_z();

    /**
     * @brief Retrieves the pressure as a FEM function.
     * @return FEMFunction2D representing the solution.
     */
    FEMFunction2D pressure();

    /**
     * @brief Retrieves the current log level of the logger.
     * 
     * @return The current log level.
     */
    int log_level();

    /**
     * @brief Sets the log level for the logger.
     * 
     * @param[in] value The log level to set.
     */
    void set_log_level(int value);

    /**
     * @brief Resets the left-hand side matrix to its initial state
     */
    void reset_lhs();

    /**
     * @brief Resets the right-hand side vector to its initial state.
     */
    void reset_rhs();

    /**
     * @brief Resets the system matrices and vectors.
     */
    void reset_system();
};

