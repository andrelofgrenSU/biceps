#/bin/env python3
import numpy as np
import argparse


def leggauss128(precision):
    # Compute the roots of the Legendre polynomial with high precision
    coeff = np.zeros(precision + 1, dtype=np.float128)
    coeff[-1] = 1  # Highest order term is 1
    roots = np.polynomial.legendre.Legendre(coeff).roots().astype(np.float128)  # Compute roots (nodes)

    # Compute weights
    Pn_derivative = np.polynomial.legendre.Legendre(coeff).deriv()
    weights = 2 / ((1 - roots**2) * (Pn_derivative(roots) ** 2))

    return roots, weights


# Setup commandline parser
arg_parser = argparse.ArgumentParser(
    description="Prints Gauss-Legendre points for the 2d reference element [0, 1]x[0, 1] up to specified precision"
)
arg_parser.add_argument(
    "--precision-max", help="", type=int, required=True
)
arg_parser.add_argument(
    "--format", help="", type=str, required=True
)
cmd_args = arg_parser.parse_args()
precision_max = cmd_args.precision_max
fmt = cmd_args.format

# Print points
print("\tif (cell_type == QUADRILATERAL) {")

for precision in range(1, precision_max+1):
    points_1d, weights_1d = leggauss128(precision)

    weights_2d = np.zeros(weights_1d.size**2)
    points_2d = np.zeros([weights_1d.size**2, 2])

    k = 0
    for zp, wz in zip(points_1d, weights_1d):
        for xp, wx in zip(points_1d, weights_1d):
            points_2d[k, :] = 0.5*(np.array([xp, zp]) + 1)
            weights_2d[k] = 0.25*wx*wz
            k += 1

    # Print points
    if (precision == 1):
        print("\t\tif (precision == 1) {")
    else:
        print(f"\t\t}} else if (precision == {precision:d}) {{")

    print("\t\t\tpoints <<")
    for i, point in enumerate(points_2d):
        if i < points_2d.shape[0]-1:
            print(f"\t\t\t\t{point[0]:{fmt}}, {point[1]:{fmt}},")
        else:
            print(f"\t\t\t\t{point[0]:{fmt}}, {point[1]:{fmt}};\n")

    # Print weights
    print("\t\t\tweights <<")
    for i, weight in enumerate(weights_2d):
        if i < weights_2d.shape[0]-1:
            print(f"\t\t\t\t{weight:{fmt}},")
        else:
            print(f"\t\t\t\t{weight:{fmt}};")
print("\t\t} else {")
print("\t\t}")
print("\t} else if (cell_type == TRIANGLE_LEFT || cell_type == TRIANGLE_RIGHT) {")
