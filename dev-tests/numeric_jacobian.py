import sympy as sy
import numpy as np

def dFxdr(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return (x0 - x1 + x2 - x3)*s + x1 - x0

def dFzdr(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return (z0 - z1 + z2 - z3)*s + z1 - z0

def dFxds(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return (x0 - x1 + x2 - x3)*r + x3 - x0

def dFzds(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return (z0 - z1 + z2 - z3)*r + z3 - z0

def T1(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return dFxdr(element, r, s)*dFzds(element, r, s)

def T2(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return dFxds(element, r, s)*dFzdr(element, r, s)

def det_J(element, r, s):
    return T1(element, r, s) - T2(element, r, s)

def T1_symbolic(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return r*s*(
        (x0-x1)*(z0-z1) + (x0-x1)*(z2-z3) + (x2-x3)*(z0-z1) + (x2-x3)*(z2-z3)
    ) + r*((x1-x0)*(z0-z1) + (x1-x0)*(z2-z3)) + \
        s*((x0-x1)*(z3-z0) + (x2-x3)*(z3-z0)) + \
        (z3-z0)*(x1-x0)

def T2_symbolic(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return r*s*(
        (x0-x1)*(z0-z1) + (x2-x3)*(z0-z1) + (x0-x1)*(z2-z3) + (x2-x3)*(z2-z3)
    ) + r*((x0-x1)*(z1-z0) + (x2-x3)*(z1-z0)) + \
        s*((x3-x0)*(z0-z1) + (x3-x0)*(z2-z3)) + \
        (x3-x0)*(z1-z0)

#def det_J_symbolic(element, r, s):
#    x0, z0 = element[0, :]
#    x1, z1 = element[1, :]
#    x2, z2 = element[2, :]
#    x3, z3 = element[3, :]
#    return r*((x1-x0)*(z2-z3) - (z1-z0)*(x2-x3)) + \
#        s*((z3-z0)*(x1-x0) - (x3-x0)*(z1-z0)) + \
#        (x1-x0)*(z3-z0) - (x3-x0)*(z1-z0)

def det_J_symbolic(element, r, s):
    x0, z0 = element[0, :]
    x1, z1 = element[1, :]
    x2, z2 = element[2, :]
    x3, z3 = element[3, :]
    return (x1-x0)*((1-r+s)*(z3-z0) + (r-s)*(z2-z1))


def area(element, gp):
    qp, wq = np.polynomial.legendre.leggauss(gp)
    area_sum = 0.0
    det_Jrs = lambda r, s: abs(det_J(element, r, s))
    for i in range(len(qp)):
        qi = 0.5*(1 + qp[i])
        for j in range(len(qp)):
            qj = 0.5*(1 + qp[j])
            area_sum += 0.25*wq[i]*wq[j]*det_Jrs(qi, qj)
    return area_sum

#element = np.array([
#    [0.0, 0.0], [2.0, -1/np.sqrt(2)],
#    [2.0, 1/np.sqrt(2)], [0.0, 1/np.sqrt(3)]
#])

element = np.array([
    [0.0, 0.0], [0.5, 0.25],
    [0.5, 0.4], [0.0, 0.8]
])

for r in np.arange(0, 1, 101):
    for s in np.arange(0, 1, 101):
        assert(abs(T1(element, r, s) - T1_symbolic(element, r, s)) < 1e-12)
        assert(abs(T2(element, r, s) - T2_symbolic(element, r, s)) < 1e-12)
        assert(abs(det_J(element, r, s) - det_J_symbolic(element, r, s)) < 1e-12)
