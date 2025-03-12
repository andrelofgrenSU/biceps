# biceps
2D full-Stokes ice-sheet solver

## BUILD INSTRUCTIONS
To build this project you need a C++ compiler (e.g., gcc or clang) and CMake:

```sh sudo apt install clang cmake```

Furthermore, the python interface requires [eigenpy](https://github.com/stack-of-tasks/eigenpy) and boost-python. Eigenpy in turn depends on numpy and scipy:

```sh sudo apt install python3-numpy python3-scipy```

Once eigenpy dependencies are installed, grab the latest release from [here] (https://github.com/stack-of-tasks/eigenpy/archive/refs/tags/v3.10.3.tar.gz) and compile it

```sh wget https://github.com/stack-of-tasks/eigenpy/archive/refs/tags/v3.10.3.tar.gz && tar -xvzf v3.10.3.tar.gz && cd v3.10.3 && mkdir .build && cmake .. -DCMAKE_BUILD_TYPE=release```

The install by running

```sh sudo make install```


WIP
