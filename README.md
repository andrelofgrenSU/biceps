# BICEPS
## Overview
Biceps is a prognostic two-dimensional full-Stokes ice-sheet solver. It is mainly intended for theoretical investigation of numerical stability. Furthermore, it comes bundled with its own FEM libraries implemented using the numerical linear algebra library [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page).

## Governing equations
Modeling ice essentially boils down to solving a highly viscous gravity-driven free-surface problem. The dynamics of the ice are described by Stokes equation (a simplified version of the Navier-Stokes equation, valid only for viscous flows)

The interface between the ice and atmosphere is modeled as freely moving boundary, and is described by the so-called free-surface equation

The two main modules are the pStokesProblem and the FreeSurfaceProblem

### Boundary conditions


## BUILD INSTRUCTIONS
In the instructions that follow, commands requiring elevated privileges are prepended with a '#' and commands that can be run as regular user with a '$'.
### Linux (Ubuntu/Debian)
This project has rather few dependencies, to build it you only need a working C++ compiler (e.g., gcc or clang), CMake and Eigen3:

```# apt install gcc cmake libeigen3-dev```

The python interface in addition requires [eigenpy](https://github.com/stack-of-tasks/eigenpy) and boost-python. Eigenpy in turn depends on numpy and scipy:

```# apt install python3-numpy python3-scipy```

Once eigenpy dependencies are installed, grab the latest release from [here](https://github.com/stack-of-tasks/eigenpy/archive/refs/tags/v3.10.3.tar.gz) and compile it:

```$ wget https://github.com/stack-of-tasks/eigenpy/archive/refs/tags/v3.10.3.tar.gz && tar -xvzf v3.10.3.tar.gz && cd v3.10.3 && mkdir .build && cmake .. -DCMAKE_BUILD_TYPE=release```

Then install it:

```# make install```

## DEMOS
Below are commented examples using C++ and Python interface. The demos are available under demos/
### C++

### Python
