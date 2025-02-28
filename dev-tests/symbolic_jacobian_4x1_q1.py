import numpy as np
import sympy as sy
sy.init_printing(use_unicode=True)

x0, x1, x2, x3, z0, z1, z2, z3, r, s, Lx, h0, h1 = sy.symbols(
    "x0 x1 x2 x3 z0 z1 z2 z3 r s Lx h0 h1"
)

nx, nz = sy.symbols("nx nz", integer=False)

phi_0 = (1-r)*(1-s)
phi_1 = r*(1-s)
phi_2 = r*s
phi_3 = (1-r)*s

phi_list = [phi_0, phi_1, phi_2, phi_3]

Fx = phi_0*x0 + phi_1*x1 + phi_2*x2 + phi_3*x3
Fz = phi_0*z0 + phi_1*z1 + phi_2*z2 + phi_3*z3

det_J = Fx.diff(r)*Fz.diff(s) - Fx.diff(s)*Fz.diff(r)
det_J = det_J.subs([
    (x3, x0), (x2, x1), (x1, x0 + Lx/nx), (x2, x0 + Lx/nx),
    (z3, z0 + h0/nz), (z2, z1 + h1/nz)
])

det_J = det_J.expand().collect(h0*Lx/(nx*nz)).subs([(nx, 4), (nz, 1), (Lx, 1)])

I = np.zeros((10, 5))
# cell 0
I[0, 0] += sy.integrate(phi_0*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[0, 1] += sy.integrate(phi_0*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[1, 0] += sy.integrate(phi_1*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[1, 1] += sy.integrate(phi_1*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[6, 0] += sy.integrate(phi_2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[6, 1] += sy.integrate(phi_2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[5, 0] += sy.integrate(phi_3*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[5, 1] += sy.integrate(phi_3*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

# cell 1
I[1, 1] += sy.integrate(phi_0*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[1, 2] += sy.integrate(phi_0*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[2, 1] += sy.integrate(phi_1*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[2, 2] += sy.integrate(phi_1*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[7, 1] += sy.integrate(phi_2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[7, 2] += sy.integrate(phi_2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[6, 1] += sy.integrate(phi_3*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[6, 2] += sy.integrate(phi_3*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

# cell 2
I[2, 2] += sy.integrate(phi_0*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[2, 3] += sy.integrate(phi_0*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[3, 2] += sy.integrate(phi_1*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[3, 3] += sy.integrate(phi_1*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[8, 2] += sy.integrate(phi_2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[8, 3] += sy.integrate(phi_2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[7, 2] += sy.integrate(phi_3*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[7, 3] += sy.integrate(phi_3*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

# cell 3
I[3, 3] += sy.integrate(phi_0*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[3, 4] += sy.integrate(phi_0*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[4, 3] += sy.integrate(phi_1*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[4, 4] += sy.integrate(phi_1*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[9, 3] += sy.integrate(phi_2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[9, 4] += sy.integrate(phi_2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

I[8, 3] += sy.integrate(phi_3*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1)).evalf()
I[8, 4] += sy.integrate(phi_3*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1)).evalf()

import IPython; IPython.embed()
