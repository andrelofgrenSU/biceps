import numpy as np
import matplotlib.pyplot as plt


def draw_triangle(corners, *plt_args, **plt_kwargs):
    for i in range(corners.shape[0]):
        x1, y1 = corners[i]
        x2, y2 = corners[(i+1) % corners.shape[0]]
        plt.plot([x1, x2], [y1, y2], *plt_args, **plt_kwargs)


def rs_to_xz_map(element, r, s):
    return np.matmul(np.array([1-r-s, r, s]), element)


def xz_to_rs_map(element, x, z):
    x1, z1 = element[0]
    x2, z2 = element[1]
    x3, z3 = element[2]
    detF = (x2-x1)*(z3-z1) - (x3-x1)*(z2-z1)
    F = np.array([
        [z3-z1, x1-x3],
        [z1-z2, x2-x1]
    ])
    return np.matmul(F/detF, np.array([x-x1, z-z1]))


if __name__ == "__main__":
    ref_tri = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    element = np.array([[1.0, 0], [1.0, 1.0], [0.0, 0.0]])
    rp, sp = 0.0, 1.0
    xp, zp = rs_to_xz_map(element, rp, sp)
    print(xz_to_rs_map(element, xp, zp))

    draw_triangle(ref_tri, color="blue")
    draw_triangle(element, color="red")
    plt.plot(rp, sp, "o", color="blue")
    plt.plot(xp, zp, "o", color="red")
    plt.show()
