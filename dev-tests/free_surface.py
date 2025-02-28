import fenics as fe

mesh = fe.IntervalMesh(2, 0, 1)
V = fe.FunctionSpace(mesh, "Lagrange", 2)

v = fe.TestFunction(V)
h = fe.TrialFunction(V)
h0 = fe.interpolate(fe.Constant(1.0), V)
ux = fe.interpolate(fe.Constant(1.0), V)
uz = fe.interpolate(fe.Constant(1.0), V)
acc = fe.interpolate(fe.Expression("x[0]", degree=1), V)

dt = 1.0
lhs_form = v*h*fe.dx
rhs_form = v*(h0 + dt*(uz + acc))*fe.dx

lhs_mat = fe.assemble(lhs_form)
rhs_vec = fe.assemble(rhs_form)

print(rhs_vec[:])
h = fe.Function(V)
fe.solve(lhs_form == rhs_form, h, [])
print(h.vector()[:])
