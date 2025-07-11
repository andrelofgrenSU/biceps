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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Stokes

#include <boost/test/unit_test.hpp>
#include <pstokes_fem.hpp>
#include <enums.hpp>
#include <math.h>

using namespace MESH2D;

inline double ux_func(double x, double z) {
    return sin(4.0*M_PI*x)*cos(4.0*M_PI*z);
}

inline double uz_func(double x, double z) {
    return -cos(4.0*M_PI*x)*sin(4.0*M_PI*z);
}

inline double uz_func_bed(double x, double z) {
    return cos(4.0*M_PI*x)*sin(4.0*M_PI*z);
}

inline double p_func(double x, double z) {
    return M_PI*cos(4.0*M_PI*x)*cos(4.0*M_PI*z);
}

inline double force_x(double x, double z) {
    return 28.0*M_PI*M_PI*sin(4.0*M_PI*x)*cos(4.0*M_PI*z);
}

inline double force_z(double x, double z) {
    return -36.0*M_PI*M_PI*cos(4.0*M_PI*x)*sin(4.0*M_PI*z);
}

static std::vector<std::tuple<int, int, int, double, double>> test_manufactured_sol_cases;

inline void populate_manufactured_sol_cases() {
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        8, 8, TRIANGLE_RIGHT, 0.05, 0.08
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        16, 8, TRIANGLE_RIGHT, 0.04, 0.03
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        8, 16, TRIANGLE_RIGHT, 0.015, 0.05
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        16, 16, TRIANGLE_RIGHT, 0.075, 0.075
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        32, 32, TRIANGLE_RIGHT, 0.005, 0.005
    ));

    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        8, 8, TRIANGLE_LEFT, 0.05, 0.08
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        16, 8, TRIANGLE_LEFT, 0.04, 0.03
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        8, 16, TRIANGLE_LEFT, 0.015, 0.05
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        16, 16, TRIANGLE_LEFT, 0.075, 0.075
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        32, 32, TRIANGLE_LEFT, 0.005, 0.005
    ));

    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        8, 8, QUADRILATERAL, 0.05, 0.08
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        16, 8, QUADRILATERAL, 0.04, 0.03
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        8, 16, QUADRILATERAL, 0.015, 0.05
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        16, 16, QUADRILATERAL, 0.075, 0.075
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, double, double>(
        32, 32, QUADRILATERAL, 0.005, 0.005
    ));
}
