#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Poisson

#include <boost/test/unit_test.hpp>
#include <enums.hpp>
#include <poisson_fem.hpp>

using namespace MESH2D;

static std::vector<std::tuple<int, int, int, int, FloatType>> test_manufactured_sol_cases;

inline void populate_manufactured_sol_cases()
{
    test_manufactured_sol_cases.clear();

    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 8, 1, TRIANGLE_RIGHT, 0.015
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 8, 1, TRIANGLE_RIGHT, 0.01
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 16, 1, TRIANGLE_RIGHT, 0.01
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 16, 1, TRIANGLE_RIGHT, 0.005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        32, 32, 1, TRIANGLE_RIGHT, 0.001
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        64, 64, 1, TRIANGLE_RIGHT, 0.0005
    ));

    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 8, 1, TRIANGLE_LEFT, 0.015
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 8, 1, TRIANGLE_LEFT, 0.01
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 16, 1, TRIANGLE_LEFT, 0.01
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 16, 1, TRIANGLE_LEFT, 0.005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        32, 32, 1, TRIANGLE_LEFT, 0.001
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        64, 64, 1, TRIANGLE_LEFT, 0.0005
    ));

    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 8, 1, QUADRILATERAL, 0.015
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 8, 1, QUADRILATERAL, 0.01
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 16, 1, QUADRILATERAL, 0.01
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 16, 1, QUADRILATERAL, 0.005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        32, 32, 1, QUADRILATERAL, 0.001
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        64, 64, 1, QUADRILATERAL, 0.0005
    ));

    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 8, 2, TRIANGLE_RIGHT, 0.0005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 8, 2, TRIANGLE_RIGHT, 0.00025
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 16, 2, TRIANGLE_RIGHT, 0.00025
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 16, 2, TRIANGLE_RIGHT, 0.00005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        32, 32, 2, TRIANGLE_RIGHT, 0.000005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        64, 64, 2, TRIANGLE_RIGHT, 0.0000001
    ));

    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 8, 2, TRIANGLE_LEFT, 0.0005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 8, 2, TRIANGLE_LEFT, 0.00025
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 16, 2, TRIANGLE_LEFT, 0.00025
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 16, 2, TRIANGLE_LEFT, 0.00005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        32, 32, 2, TRIANGLE_LEFT, 0.000005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        64, 64, 2, TRIANGLE_LEFT, 0.0000001
    ));

    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 8, 2, QUADRILATERAL, 0.0005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 8, 2, QUADRILATERAL, 0.00025
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        8, 16, 2, QUADRILATERAL, 0.00025
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        16, 16, 2, QUADRILATERAL, 0.00005
    ));
    test_manufactured_sol_cases.push_back(std::tuple<int, int, int, int, FloatType>(
        32, 32, 2, QUADRILATERAL, 0.000005
    ));
}

inline FloatType alpha_manufactured_case_1(FloatType x, FloatType z)
{
    return 1.0;
}

inline FloatType beta_manufactured_case_1(FloatType x, FloatType z)
{
    return 1.0;
}

inline FloatType force_manufactured_case_1(FloatType x, FloatType z)
{
    return 2.0*PI_CONST*PI_CONST*SIN_FUNC(PI_CONST*x)*SIN_FUNC(PI_CONST*z);
}

inline FloatType bc_func_manufactured_case_1(FloatType x, FloatType z)
{
    return 0.0;
}

inline FloatType exact_sol_func_manufactured_case_1(FloatType x, FloatType z)
{
    return SIN_FUNC(PI_CONST*x)*SIN_FUNC(PI_CONST*z);
}

inline FloatType alpha_manufactured_case_2(FloatType x, FloatType z)
{
    return 1.0;
}

inline FloatType beta_manufactured_case_2(FloatType x, FloatType z)
{
    return 1.0;
}

inline FloatType gamma_manufactured_case_2(FloatType x, FloatType z)
{
    return 1.0;
}

inline FloatType force_manufactured_case_2(FloatType x, FloatType z)
{
    return (1.0 + 2.0*PI_CONST*PI_CONST)*SIN_FUNC(PI_CONST*x)*SIN_FUNC(PI_CONST*z);
}

inline FloatType bc_func_manufactured_case_2(FloatType x, FloatType z)
{
    return 0.0;
}

inline FloatType exact_sol_func_manufactured_case_2(FloatType x, FloatType z)
{
    return SIN_FUNC(PI_CONST*x)*SIN_FUNC(PI_CONST*z);
}

inline FloatType alpha_manufactured_case_3(FloatType x, FloatType z)
{
    return 3.0;
}

inline FloatType beta_manufactured_case_3(FloatType x, FloatType z)
{
    return 2.0;
}

inline FloatType gamma_manufactured_case_3(FloatType x, FloatType z)
{
    return -2.0;
}

inline FloatType force_manufactured_case_3(FloatType x, FloatType z)
{
    return (5.0*PI_CONST*PI_CONST - 2.0)*SIN_FUNC(PI_CONST*x)*SIN_FUNC(PI_CONST*z);
}

inline FloatType bc_func_manufactured_case_3(FloatType x, FloatType z)
{
    return 0.0;
}

inline FloatType exact_sol_func_manufactured_case_3(FloatType x, FloatType z)
{
    return SIN_FUNC(PI_CONST*x)*SIN_FUNC(PI_CONST*z);
}
