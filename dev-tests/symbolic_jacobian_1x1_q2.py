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

Fx = phi_0*x0 + phi_1*x1 + phi_2*x2 + phi_3*x3
Fz = phi_0*z0 + phi_1*z1 + phi_2*z2 + phi_3*z3

det_J = Fx.diff(r)*Fz.diff(s) - Fx.diff(s)*Fz.diff(r)
det_J = det_J.subs([
    (x3, x0), (x2, x1), (x1, x0 + Lx/nx), (x2, x0 + Lx/nx),
    (z3, z0 + h0/nz), (z2, z1 + h1/nz)
])

phi_0_q2 = 4*(0.5-r)*(1-r)*(1-s)*(0.5-s)
phi_1_q2 = 8*r*(1-r)*(1-s)*(0.5-s)
phi_2_q2 = 4*r*(r-0.5)*(0.5-s)*(1-s)
phi_3_q2 = 8*r*(r-0.5)*s*(1-s)
phi_4_q2 = 4*r*(r-0.5)*s*(s-0.5)
phi_5_q2 = 8*r*(1-r)*s*(s-0.5)
phi_6_q2 = 4*(0.5-r)*(1-r)*s*(s-0.5)
phi_7_q2 = 8*(0.5-r)*(1-r)*s*(1-s)
phi_8_q2 = 16*r*(1-r)*s*(1-s)

det_J = det_J.expand().collect(h0*Lx/(nx*nz))
I00 = sy.integrate(phi_0_q2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1))
I01 = sy.integrate(phi_0_q2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1))

I10 = sy.integrate(phi_1_q2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1))
I11 = sy.integrate(phi_1_q2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1))

I20 = sy.integrate(phi_2_q2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1))
I21 = sy.integrate(phi_2_q2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1))

I30 = sy.integrate(phi_3_q2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1))
I31 = sy.integrate(phi_3_q2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1))

I40 = sy.integrate(phi_4_q2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1))
I41 = sy.integrate(phi_4_q2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1))

I50 = sy.integrate(phi_5_q2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1))
I51 = sy.integrate(phi_5_q2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1))

I60 = sy.integrate(phi_6_q2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1))
I61 = sy.integrate(phi_6_q2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1))

I70 = sy.integrate(phi_7_q2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1))
I71 = sy.integrate(phi_7_q2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1))

I80 = sy.integrate(phi_8_q2*det_J.coeff(h0, 1), (r, 0, 1), (s, 0, 1))
I81 = sy.integrate(phi_8_q2*det_J.coeff(h1, 1), (r, 0, 1), (s, 0, 1))
import IPython; IPython.embed()
