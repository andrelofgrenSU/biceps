#!/bin/sh

BUILD_DIR=.build
mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release
make -j8 && make test_mesh && make test_fem_2d && make test_poisson && make test_stokes && sudo make install
