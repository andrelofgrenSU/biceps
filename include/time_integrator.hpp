#pragma once
#include <pstokes_fem.hpp>
#include <free_surface_fem.hpp>
#include <Eigen/Dense>

/**
 * @class TimeIntegrator
 * @brief A class for coupling the pStokesProblem and the Free-surface problem.
 *
 * This class provides predefined time-stepping methods for the pStokes-coupled free-surface problem.
*/
class TimeIntegrator {
private:

public:
    pStokesProblem &psp;  ///< pStokesProblem
    FreeSurfaceProblem &fsp;  ///< pStokesProblem
                              ///
    /**
     * @brief Constructor for initializing TimeIntegrator
     * 
     * @param[in] psp pStokesProblem
     * @param[in] fsp FreeSurfaceProblem
     */
    TimeIntegrator(pStokesProblem &psp, FreeSurfaceProblem &fsp);
    Eigen::VectorX<FloatType> step_explicit(FloatType dt);
    void extrude_mesh_z(const Eigen::VectorX<FloatType> &zs_vec);
};
