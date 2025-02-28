import numpy as np
import matplotlib.pyplot as plt


def draw_quad(corners, *plt_args, **plt_kwargs):
    for i in range(corners.shape[0]):
        x1, y1 = corners[i]
        x2, y2 = corners[(i+1) % corners.shape[0]]
        plt.plot([x1, x2], [y1, y2], *plt_args, **plt_kwargs)


def rs_to_xz_map(element, r, s):
    return np.matmul(np.array([(1-r)*(1-s), r*(1-s), r*s, s*(1-r)]), element)


def jac_F(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return np.array([
        [(s-1)*x0 + (1-s)*x1 + s*x2 - s*x3, (r-1)*x0 - r*x1 + r*x2 + (1-r)*x3],
        [(s-1)*z0 + (1-s)*z1 + s*z2 - s*z3, (r-1)*z0 - r*z1 + r*z2 + (1-r)*z3]
    ])


def xz_to_rs_map(element, x, z):
    x0, z0 = element[0]
    x1, z1 = element[1]
    x2, z2 = element[2]
    x3, z3 = element[3]
    detF = (x2-x1)*(z3-z1) - (x3-x1)*(z2-z1)
    F = np.array([
        [z3-z1, x1-x3],
        [z1-z2, x2-x1]
    ])
    return np.matmul(F/detF, np.array([x-x1, z-z1]))


def calculate_area(element, gp):
    qp, wq = np.polynomial.legendre.leggauss(gp)
    area_sum = 0.0
    det_J = lambda r, s: abs(np.linalg.det(jac_F(element, r, s)))
    for i in range(len(qp)):
        for j in range(len(qp)):
            area_sum += 0.25*wq[i]*wq[j]*det_J(0.5*(1 + qp[i]), 0.5*(1 + qp[j]))
    return area_sum


if __name__ == "__main__":
    ref_quad = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    element = np.array([[0.3, 0.3], [0.5, -0.1], [0.5, 0.5], [0.3, 0.7]]) + \
        np.array([0.3, 0.3])
    element = np.array([[0.0, 0.0], [2.0, -1/np.sqrt(2)], [2.0, 1/np.sqrt(2)],
                        [0.0, 1/np.sqrt(3)]])
    rp_0, sp_0 = 0.0, 0.0
    rp_1, sp_1 = 1.0, 0.0
    rp_2, sp_2 = 1.0, 1.0
    rp_3, sp_3 = 0.0, 1.0

    rp_mid, sp_mid = 0.5, 0.5

    xp_0, zp_0 = rs_to_xz_map(element, rp_0, sp_0)
    xp_1, zp_1 = rs_to_xz_map(element, rp_1, sp_1)
    xp_2, zp_2 = rs_to_xz_map(element, rp_2, sp_2)
    xp_3, zp_3 = rs_to_xz_map(element, rp_3, sp_3)

    xp_mid, zp_mid = rs_to_xz_map(element, rp_mid, sp_mid)

    element = np.array([[0.0, 0.0], [2.0, -1/np.sqrt(2)], [2.0, 1/np.sqrt(2)],
                        [0.0, 1/np.sqrt(3)]])
    print(np.linalg.det(jac_F(element, 0.25, 0.3)))
    print(np.linalg.det(jac_F(element, 0.7, 0.5)))
    A = calculate_area(element, 1)

    draw_quad(ref_quad, color="blue")
    draw_quad(element, color="red")
    plt.plot(rp_0, sp_0, "o", color="blue")
    plt.plot(xp_0, zp_0, "o", color="red")
    plt.plot(rp_1, sp_1, "s", color="blue")
    plt.plot(xp_1, zp_1, "s", color="red")
    plt.plot(rp_2, sp_2, "^", color="blue")
    plt.plot(xp_2, zp_2, "^", color="red")
    plt.plot(rp_3, sp_3, "*", color="blue")
    plt.plot(xp_3, zp_3, "*", color="red")
    plt.plot(rp_mid, sp_mid, "P", color="blue")
    plt.plot(xp_mid, zp_mid, "P", color="red")
    plt.show()
