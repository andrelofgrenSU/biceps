# About
Biceps is a prognostic two-dimensional full-Stokes ice-sheet solver. It is mainly intended for theoretical investigation of numerical stability. Furthermore, it comes bundled with its own FEM libraries implemented using the numerical linear algebra library [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page).

# Build instruction
In the following instructions, note that commands requiring elevated privileges are prepended with a '#', and commands that can be run as regular user with a '$'.
## Linux (Ubuntu/Debian)

### C++
This project has rather few dependencies, a minimal C++ installation requires only a working C++ compiler and tool chain (e.g., [gcc](https://gcc.gnu.org/)), [CMake](https://cmake.org/), [boost](https://www.boost.org/), and [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page):

```# apt install gcc build-essentials cmake libboost-dev libeigen3-dev```

### Python

Building the python interface requires [EigenPy](https://github.com/stack-of-tasks/eigenpy) and the python module of [boost](https://www.boost.org/). [EigenPy](https://github.com/stack-of-tasks/eigenpy) in turn depends on [NumPy](https://numpy.org/) and [SciPy](https://scipy.org/):

```# apt install python3-numpy python3-scipy libboost-python-dev```

After installing dependencies for [EigenPy](https://github.com/stack-of-tasks/eigenpy), grab the latest release from [here](https://github.com/stack-of-tasks/eigenpy/archive/refs/tags/v3.10.3.tar.gz) and compile it:

```$ tar -xvzf v3.10.3.tar.gz && cd v3.10.3 && mkdir -p .build && cd .build && cmake .. -DCMAKE_BUILD_TYPE=release```

Then to install run:

```# make install```

CMake flag: -DWITH_PYTHON:boolean

### Documentation

Generating documentation requires [Doxygen](https://www.doxygen.nl/) and [Graphviz](https://graphviz.org/):

```# apt install doxygen graphviz```

CMake flag: -DBUILD_DOCS:boolean

### Example build

The steps for installing this repository is:
1. Clone repository: ```git clone https://github.com/andrelofgrenSU/biceps.git```

2. Compile: ```$ cd biceps && mkdir -p .build && cd .build && cmake .. -DCMAKE_BUILD_TYPE=release -DWITH_PYTHON=ON -DBUILD_DOCS=OFF -DUSE_LONG_DOUBLE=OFF -DTESTS=OFF```

3. Run tests (optional): ```# make install```

4. Install: ```# make install```

# Governing equations

## The pStokes equations

### Strong formulation
The velocity- and pressure distribution, denoted $\mathbf{u}$ and $p$ respectively, inside the ice sheet is governed by the Stokes equation (a simplification of the Navier-Stokes equation, valid only for viscous dominated flows), which consists of the momentum balance and an incompressibility condition:

$$
\begin{aligned}
    \nabla \cdot (2 \eta(\mathbf{u}) \dot{\varepsilon}(\mathbf{u})) - \nabla p + \mathbf{f} &= \mathbf{0}, \quad \mathbf{x} \in \Omega, \\
    \nabla \cdot \mathbf{u} &= 0, \quad \mathbf{x} \in \Omega.
\end{aligned}
$$

Here $\dot{\varepsilon}(\mathbf{u}) = \frac{1}{2} \left ( \nabla \mathbf{u} + \nabla{\mathbf{u}}^T \right )$ is the strain-rate tensor, $\mathbf{f}$ is the volumetric external body force acting on each fluid element (gravity in case of ice). Furthermore, ice is a shear thinning fluid with the viscosity function $\eta$ following a power-law rheology known as Glen's flow law


$$\eta(\mathbf{u}) = A^{\frac{1}{n}} \left (\dot{\varepsilon}^2_e (\mathbf{u}) + \dot{\varepsilon}^2_0 \right)^{\frac{1-n}{2n}}.$$

Here $A$ is the so-called rate factor or ice softness parameter, $n \approx 3$ is the glen exponent, and $\dot{\varepsilon}_e$ is the effective strain rate

$$\dot{\varepsilon}^2_e = \frac{1}{2} \left (\text{tr} \left (\dot{\varepsilon}^2 \right ) - \text{tr}^2 \left (\dot{\varepsilon} \right ) \right)$$

In addition, a small regularization term $\varepsilon_0$ is included to prevent infinite viscosity at zero effective strain rate.

This type of power-law viscosity coupled Stokes equation is referred to as the pStokes equation.

### Boundary conditions
For wellposedness of the pStokes equation, boundary conditions needs to be specified on all parts of the boundary. Specifically, for grounded ice sheets the following are typically imposed:

1. Stress free condition on the ice-atmosphere interface $\Gamma_s$: $\sigma \hat{\mathbf{n}} = \mathbf{0}$
2. No-slip on the bedrock $\Gamma_b$: $\mathbf{u} \cdot \hat{\mathbf{n}} = 0$
3. Impenetrability on the lateral boundaries $\Gamma_E$ and $\Gamma_W$: $\mathbf{u} \cdot \hat{\mathbf{n}} = 0$

These boundary conditions can, however, be used interchangeably on all boundary parts.

### Weak formulation
The pStokes equations are solved using the finite element method (FEM), which discretizes the weak formulation. To state the weak form, the momentum equation and the incompressibility condition are multiplied by test functions $\mathbf{v}$ and $q$, respectively, and then integrating over the domain $\Omega$, resulting in the following weak formulation:

Find $\mathbf{u} \in U$ and $p \in q$, such that

```math
\left (\dot{\varepsilon}(\mathbf{v}), 2 \eta(\mathbf{u}) \dot{\varepsilon}(\mathbf{u}) \right )_{\Omega} - (\nabla \cdot \mathbf{v}, p) - (q, \nabla \cdot \mathbf{u})_{\Omega} = (\mathbf{v}, \mathbf{f})_{\Omega}
```

for all $\mathbf{v} \in V$ and all $q \in Q$. Here $U, V$ and $Q$ are appropriate Sobolev spaces, in particular the discretized trial spaces $U$ and $Q$ should be chosen so that they satisfy a discrete *inf-sup* stability condition. A common choice is so-called Taylor-Hood element, using quadratic basis functions for $U$ and linear basis for $Q$.

### Nonlinear iterations
To resolve the nonlinearity a Picard iteration scheme is employed where $2 \eta(\mathbf{u}) \dot{\varepsilon}(\mathbf{u}) \approx 2 \eta(\mathbf{u}_0) \dot{\varepsilon}(\mathbf{u})$, with $\mathbf{u}_0$ being some known approximation of $\mathbf{u}$. The following problem is then solved to obtain an improved guess $\mathbf{u}^{m+1}$ from the known guess $\mathbf{u}^m$:

Find $\mathbf{u}^{m+1} \in U$ and $p \in q$, such that

```math
\begin{equation}
    \left (\dot{\varepsilon}(\mathbf{v}), 2 \eta(\mathbf{u}^m) \dot{\varepsilon}(\mathbf{u}^{m+1}) \right )_{\Omega} - (q^{m+1}, \nabla \cdot \mathbf{u}^{m+1})_{\Omega} - (\nabla \cdot \mathbf{v}, p^{m+1})_{\Omega} = (\mathbf{v}, \mathbf{f})_{\Omega}
\end{equation}
```

for all $\mathbf{v} \in V$ and all $q \in Q$. This is then iterated upon until a user defined step tolerance $\epsilon_s$ is reached:

```math
\begin{equation}
    || \mathbf{u}^{m+1} - \mathbf{u}^{m}||_{L^2(\Omega)} < \epsilon_s || \mathbf{u}^{m+1}||_{L^2(\Omega)}
\end{equation}
```

## The free-surface equation
The interface between the ice and atmosphere is modeled as freely moving boundary, this interface moves either due to ice particles being transported by the ice flow across the boundary, or due snow accumulating or ablating on top of it. Tracking the free-surface height $h(x, t)$, its evolution until time $T$ is described by the so-called free-surface equation

```
\begin{equation}
\frac{\partial h}{\partial t} + u_x^s(h(x, t)) \frac{\partial h}{\partial x} = u_z(h(x, t))^s + a_s(x, t), \quad (x, t) \in \Gamma^{\perp}_s \times [0, T],
\end{equation}
```

where $u^s_x$ and $u^s_z$ are the respective horizontal and vertical ice surface velocities, a_s is the surface mass balance, and \Gamma_s^{\perp} is the projection of surface onto the line $z = 0$.

### Weak formulation
Similarly, in this case the weak formulation is derived by multiplying by a test function $w$ and integrating over $\Gamma_s^{\perp}$, resulting in the following problem:

Find $h \in Z$ such that

```
\left (w, \frac{\partial h}{\partial t} \right)_{\Gamma_s^{\perp}} + \left (w, u_x^s \frac{\partial h}{\partial x} \right)_{\Gamma_s^{\perp}} = \left (w, u^s_z + a_s\right )_{\Gamma_s^{\perp}}
```

for all $w \in Z$, where $Z$ is an appropriate Sobolev space.

### Time discretization
This equation is then numerically integrated in time using either explicit- or semi-implicit time stepping. The weak formulation in each case is as follows:

##### Explicit
Find $h^{k+1} \in Z$ such that

```math
\left (w,  h^{k+1}\right)_{\Gamma_s^{\perp}} = \left (w, h^k\right)_{\Gamma_s^{\perp}} - \Delta t \left (w, u^s_x \frac{\partial h^k}{\partial x})_{\Gamma_s^{\perp}} + \Delta t (w, u^s_z + a_s \right)_{\Gamma_s^{\perp}}
```

for all $w \in Z$.

##### Semi implicit
Find $h^{k+1} \in Z$ such that

```math
\left (w,  h^{k+1}\right)_{\Gamma_s^{\perp}} + \Delta t \left (w, u^s_x \frac{\partial h^{k+1}}{\partial x} \right)_{\Gamma_s^{\perp}} = \Delta t \left (w, u^s_z\right )_{\Gamma_s^{\perp}} + \Delta t\left (w, a_s \right)_{\Gamma_s^{\perp}}
```

for all $w \in Z$.

### Free-surface stabilization algorithm (FSSA)
The pStokes equations coupled to the free-surface equation can be viewed as infinite dimensional dynamical system.

The modules for solving these equations are the pStokesProblem and the FreeSurfaceProblem. It is then up to the user to couple these equations. Full working examples on how this is done is provided in the demos.

# DEMOS
Below are commented examples using the C++ and Python interfaces. The demos can be found under demos/

## C++
```C++

```

## Python

Running this demo also requires [Matplotlib](https://matplotlib.org). To install on debian-based system, run:
```# apt install python3-matplotlib```

```python
import numpy as np
import matplotlib.pyplot as plt
import biceps as bp

GRAVITY = 9.8  # Gravitational acceleration (m/s^2)
ICE_DENSITY = 910  $ Ice density (kg/m^3)
FORCE_Z = -1e-3*ICE_DENSITY*GRAVITY  $ Force in z direction converted to units in MPa, yr and km

# Define the expressions for the bottom and surface elevations
def zb_expr(x): return 0.0  # Flat bottom surface
def zs_expr(x): return H + z0*np.cos(np.pi*x/L)  # Undulating surface elevation


# Define external force functions
def force_x(x, z): return 0.0  # No horizontal force
def force_z(x, z): return FORCE_Z  # Gravity-driven force in the z-direction


# Define domain and grid parameters
x0 = 0.0  # Left end of domain (km)
x1 = 100.0  # Right end of domain (km)
L = x1 - x0  # Length of the domain
H = 1  # Mean height of the domain (km)
z0 = 0.5  # Amplitude of surface undulation (km)

A = 100.0  # Rate factor
n_i = 3.0  # Glen exponent
fssa_version = bp.FSSA_VERSION.FSSA_NONE  # No FSSA stabilization
fssa_param = 0  # No additional parameter for FSSA
eps_reg_2 = 1e-10  # Regularization parameter

nx = 50  # Number of elements in x-direction
nz = 5  # Number of elements in z-direction
nt = 10  # Number of time steps
dt = 1  # Time step size

deg_u = 2  # Polynomial degree for velocity field
deg_p = 1  # Polynomial degree for pressure field
deg_h = 1  # Polynomial degree for height field
gp = 5  # Number of Gauss points per element

max_iter = 100  # Maximum number of iterations for solver
stol = 1e-10  # Convergence tolerance for solver

cell_type = bp.CELL_TYPE_2D.TRIANGLE_LEFT  # Type of triangular element

# Create structured meshes for velocity, pressure, and height fields
u_mesh_2d = bp.StructuredMesh(nx, nz, deg_u, cell_type)
p_mesh_2d = bp.StructuredMesh(nx, nz, deg_p, cell_type)

# Extrude the meshes in the x and z directions
u_mesh_2d.extrude_x(x0, x1)
u_mesh_2d.extrude_z(zb_expr, zs_expr)
p_mesh_2d.extrude_x(x0, x1)
p_mesh_2d.extrude_z(zb_expr, zs_expr)

# Extract degrees of freedom for surface nodes
sdofs_u = u_mesh_2d.extract_dof_inds(bp.DOMAIN_IDS_2D.SURFACE_ID)
sdofs_h = u_mesh_2d.extract_vertex_dof_inds(bp.DOMAIN_IDS_2D.SURFACE_ID)

# Extract surface coordinates
spmat_u = u_mesh_2d.pmat[sdofs_u, :].copy()
spmat_h = u_mesh_2d.pmat[sdofs_h, :].copy()
xs_vec = spmat_h[:, 0].copy()
zs_vec = spmat_h[:, 1].copy()
spmat_u[:, 1] = 0
spmat_h[:, 1] = 0

# Create 1D meshes for velocity and height on the surface
u_mesh_1d = bp.IntervalMesh(spmat_u, deg_u)
h_mesh_1d = bp.IntervalMesh(spmat_h, deg_h)

# Apply boundary conditions
ux_boundary_id = (
    bp.DOMAIN_IDS_2D.NORTH_WEST_ID |
    bp.DOMAIN_IDS_2D.WEST_ID |
    bp.DOMAIN_IDS_2D.BED_ID |
    bp.DOMAIN_IDS_2D.EAST_ID |
    bp.DOMAIN_IDS_2D.NORTH_EAST_ID
)
uz_boundary_id = bp.DOMAIN_IDS_2D.BED_ID

# Initialize FEM functions for velocity and height
ux_func = bp.FEMFunction1D(u_mesh_1d)
uz_func = bp.FEMFunction1D(u_mesh_1d)
h0_func = bp.FEMFunction1D(h_mesh_1d)
ac_func = bp.FEMFunction1D(h_mesh_1d)

# Initialize Stokes and Free Surface Problems
psp = bp.pStokesProblem(u_mesh_2d, p_mesh_2d)
fsp = bp.FreeSurfaceProblem(h_mesh_1d, u_mesh_1d)

# Plot initial surface profile
plt.plot(xs_vec, zs_vec)

# Time-stepping loop
for k in range(nt):
    # Assemble and solve the nonlinear pStokes problem
    psp.solve_nonlinear_system_picard(
        A, n_i, eps_reg_2, fssa_version, fssa_param, force_x, force_z,
        ux_boundary_id, uz_boundary_id, max_iter, stol, gp
    )

    # Extract velocity field solutions
    ux_vec = psp.velocity_x().vals
    uz_vec = psp.velocity_z().vals
    ux_func.vals = ux_vec[sdofs_u]
    uz_func.vals = uz_vec[sdofs_u]
    h0_func.vals = zs_vec  # Update initial height values

    # Print surface energy and domain area
    print(f"||E|| = {h0_func.L2_norm()}")
    print(f"A = {u_mesh_2d.area()}")

    # Solve the free surface problem using explicit time stepping
    fsp.assemble_lhs_explicit(gp)
    fsp.commit_lhs()
    fsp.assemble_rhs_explicit(h0_func, ux_func, uz_func, ac_func, dt, gp)
    fsp.solve_linear_system()

    # Update surface elevation
    zs_vec = fsp.zs_vec

    # Update mesh with new surface elevation
    u_mesh_2d.extrude_z(zs_vec)
    p_mesh_2d.extrude_z(zs_vec)

    # Plot updated surface profile
    plt.plot(xs_vec, zs_vec)

# Show final plot
plt.show()
```
