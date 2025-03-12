#pragma once
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_spmat_addons.hpp>

#include <vector>
#include <fem_function_1d.hpp>
#include <Eigen/Sparse>
#include <enums.hpp>

/**
 * @class FreeSurfaceProblem
 * @brief Solves a free surface problem with FEM for height and velocity fields.
 *
 * This class manages the solution of the free surface problem, including the assembly of the left-hand side (LHS)
 * and right-hand side (RHS) matrices, as well as the solution of the linear system.
 */
class FreeSurfaceProblem {

private:
    /**
     * @brief Stores the coefficients for the left-hand side (LHS) matrix in sparse format.
     */
    std::vector<Eigen::Triplet<FloatType>> lhs_coeffs; ///< Coefficients for LHS matrix.

    /**
     * @brief Sparse matrix representing the left-hand side (LHS) system matrix.
     */
    Eigen::SparseMatrix<FloatType> lhs_mat; ///< The sparse LHS matrix of the system.

    /**
     * @brief Right-hand side (RHS) vector of the system.
     */
    Eigen::VectorX<FloatType> rhs_vec; ///< The RHS vector.

public:
    /**
     * @brief Mesh for the height field in the free surface problem.
     */
    IntervalMesh &h_mesh; ///< The height mesh for the free surface problem.

    /**
     * @brief Mesh for the velocity field in the free surface problem.
     */
    IntervalMesh &u_mesh; ///< The velocity mesh for the free surface problem.

    /**
     * @brief Vector storing the height values of the free surface.
     */
    Eigen::VectorX<FloatType> zs_vec; ///< Free surface height values.

    /**
     * @brief Number of elements pushed for the LHS matrix assembly.
     */
    int nof_pushed_elements = 0; ///< Number of elements pushed into the LHS matrix.

    /**
     * @brief Constructor to initialize the free surface problem.
     * 
     * @param[in] h_mesh The mesh for the height function.
     * @param[in] u_mesh The mesh for the velocity function.
     */
    FreeSurfaceProblem(IntervalMesh &h_mesh, IntervalMesh &u_mesh);

    /**
     * @brief Assembles the left-hand side (LHS) matrix for the explicit scheme.
     * 
     * @param[in] gauss_precision The number of Gauss points for integration.
     */
    void assemble_lhs_explicit(int gauss_precision);

    /**
     * @brief Assembles the right-hand side (RHS) vector for the explicit scheme.
     * 
     * @param[in] h0_fem_func The FEM function for the initial height field.
     * @param[in] ux_fem_func The FEM function for the velocity in the x-direction.
     * @param[in] uz_fem_func The FEM function for the velocity in the z-direction.
     * @param[in] ac_fem_func The FEM function for the acceleration.
     * @param[in] dt The time step for the simulation.
     * @param[in] gauss_precision The number of Gauss points for integration.
     */
    void assemble_rhs_explicit(
        FEMFunction1D &h0_fem_func,
        FEMFunction1D &ux_fem_func,
        FEMFunction1D &uz_fem_func,
        FEMFunction1D &ac_fem_func,
        FloatType dt,
        int gauss_precision
    );

    /**
     * @brief Assembles the left-hand side (LHS) matrix for the implicit scheme.
     * 
     * @param[in] ux_fem_func The FEM function for the velocity in the x-direction.
     * @param[in] dt The time step for the simulation.
     * @param[in] gauss_precision The number of Gauss points for integration.
     */
    void assemble_lhs_simplicit(
        FEMFunction1D &ux_fem_func,
        FloatType dt,
        int gauss_precision
    );

    /**
     * @brief Assembles the right-hand side (RHS) vector for the implicit scheme.
     * 
     * @param[in] h0_fem_func The FEM function for the initial height field.
     * @param[in] uz_fem_func The FEM function for the velocity in the z-direction.
     * @param[in] ac_fem_func The FEM function for the acceleration.
     * @param[in] dt The time step for the simulation.
     * @param[in] gauss_precision The number of Gauss points for integration.
     */
    void assemble_rhs_simplicit(
        FEMFunction1D &h0_fem_func,
        FEMFunction1D &uz_fem_func,
        FEMFunction1D &ac_fem_func,
        FloatType dt,
        int gauss_precision
    );

    /**
     * @brief Commits the LHS matrix by converting it from triplet form to sparse matrix.
     */
    void commit_lhs();

    /**
     * @brief Solves the linear system of equations using LU decomposition.
     */
    void solve_linear_system();

    /**
     * @brief Retrieves the height function as an FEM function.
     * 
     * @return FEMFunction1D The computed height function for the free surface.
     */
    FEMFunction1D height();

    /**
     * @brief Resets the LHS matrix coefficients.
     */
    void reset_lhs();

    /**
     * @brief Resets the RHS vector.
     */
    void reset_rhs();

    /**
     * @brief Resets the entire system, including both LHS and RHS.
     */
    void reset_system();
};

