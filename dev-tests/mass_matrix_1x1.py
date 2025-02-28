import numpy as np
import sympy as sy
sy.init_printing(use_unicode=True)

x0, x1, x2, x3, z0, z1, z2, z3, r, s, Lx, h0, h1, dt = sy.symbols(
    "x0 x1 x2 x3 z0 z1 z2 z3 r s Lx h0 h1, dt"
)

I = sy.Matrix([[1, 0], [0, 1]])

phi_0 = (1-r)*(1-s)
phi_1 = r*(1-s)
phi_2 = r*s
phi_3 = (1-r)*s

Fz = h0*(1-r)*s + h1*r*s
f = 1-Fz**2

u0 = f.subs([(s, 0), (r, 0)])
u1 = f.subs([(s, 0), (r, 1)])

M20 = Lx*sy.integrate((1-r)*phi_3*phi_0*h0 + r*phi_3*phi_0*h1, (r, 0, 1), (s, 0, 1))
M21 = Lx*sy.integrate((1-r)*phi_3*phi_1*h0 + r*phi_3*phi_1*h1, (r, 0, 1), (s, 0, 1))
M22 = Lx*sy.integrate((1-r)*phi_3*phi_3*h0 + r*phi_3*phi_3*h1, (r, 0, 1), (s, 0, 1))
M23 = Lx*sy.integrate((1-r)*phi_3*phi_2*h0 + r*phi_3*phi_2*h1, (r, 0, 1), (s, 0, 1))
M30 = Lx*sy.integrate((1-r)*phi_2*phi_0*h0 + r*phi_2*phi_0*h1, (r, 0, 1), (s, 0, 1))
M31 = Lx*sy.integrate((1-r)*phi_2*phi_1*h0 + r*phi_2*phi_1*h1, (r, 0, 1), (s, 0, 1))
M32 = Lx*sy.integrate((1-r)*phi_2*phi_3*h0 + r*phi_2*phi_3*h1, (r, 0, 1), (s, 0, 1))
M33 = Lx*sy.integrate((1-r)*phi_2*phi_2*h0 + r*phi_2*phi_2*h1, (r, 0, 1), (s, 0, 1))
M = sy.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [M20, M21, M22, M23], [M30, M31, M32, M33]])

b2 = Lx*sy.integrate((1-r)*f*phi_3*h0 + r*f*phi_3*h1, (r, 0, 1), (s, 0, 1))
b3 = Lx*sy.integrate((1-r)*f*phi_2*h0 + r*f*phi_2*h1, (r, 0, 1), (s, 0, 1))
b = sy.Matrix([[1], [1], [b2], [b3]])

M_2x2 = M[2:4, 2:4]
b_2x1 = sy.Matrix([
    b2 - M20*u0 - M21*u1,
    b3 - M30*u0 - M31*u1
])

P00 = Lx*sy.integrate(f*phi_3*(1-r) - u0*(1-r)*phi_3*phi_0 - u1*(1-r)*phi_3*phi_1, (r, 0, 1), (s, 0, 1))
P01 = Lx*sy.integrate(f*phi_3*r - u0*r*phi_3*phi_0 - u1*r*phi_3*phi_1, (r, 0, 1), (s, 0, 1))
P10 = Lx*sy.integrate(f*phi_2*(1-r) - u0*(1-r)*phi_2*phi_0 - u1*(1-r)*phi_2*phi_1, (r, 0, 1), (s, 0, 1))
P11 = Lx*sy.integrate(f*phi_2*r - u0*r*phi_2*phi_0 - u1*r*phi_2*phi_1, (r, 0, 1), (s, 0, 1))
P = sy.Matrix([[P00, P01], [P10, P11]])
Q = M_2x2.inv()*P
Q.simplify()
dQ_dh0 = Q.diff(h0)
dQ_dh0.simplify()
dQ_dh1 = Q.diff(h1)
dQ_dh1.simplify()

h = sy.Matrix([[h0], [h1]])
F = Q*h
F.simplify()

dQxdh0xh = dQ_dh0*h
dQxdh1xh = dQ_dh1*h
dQ_dhxh = sy.Matrix([[dQxdh0xh, dQxdh1xh]])
dF = Q.transpose() + dQ_dhxh.transpose()
dF.simplify()

eig_vals = list(dF.eigenvals())
import IPython; IPython.embed()
