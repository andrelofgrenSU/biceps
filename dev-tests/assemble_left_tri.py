import numpy as np
import fenics as fe

mesh = fe.UnitSquareMesh(1, 1)

V = fe.FunctionSpace(mesh, "Lagrange", 2)
Q = fe.FunctionSpace(mesh, "Lagrange", 1)
W = fe.FunctionSpace(
    mesh,
    fe.VectorElement("Lagrange", "triangle", 2) *
    fe.FiniteElement("Lagrange", "triangle", 1)
)
y = fe.TrialFunction(W)
z = fe.TestFunction(W)
u, p = fe.split(y)
v, q = fe.split(z)
vx, vz = fe.split(v)
ux, uz = fe.split(u)

eta = 1.0
a = 2.0*eta*fe.inner(fe.sym(fe.grad(u)), fe.grad(v))*fe.dx
b = (-fe.div(v)*p - q*fe.div(u))*fe.dx

rho = 910
g = 9.8
f = fe.Constant(((0.0),(-1e-3*rho*g)))
l = fe.dot(f, v)*fe.dx
lhs_array = fe.assemble(a+b).array().flatten()
rhs_array = fe.assemble(l).get_local()

nnz_inds = np.abs(lhs_array) > 1e-12
nnz_mat_sorted = np.sort(lhs_array[nnz_inds])

print("LHS MAT")
for nnz in nnz_mat_sorted:
    print(f"{nnz:.3g}")

nnz_inds = np.abs(rhs_array) > 1e-12
nnz_vec_sorted = np.sort(rhs_array[nnz_inds])
print("RHS VEC")
for nnz in nnz_vec_sorted:
    print(f"{nnz:.3g}")
