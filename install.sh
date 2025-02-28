#!/bin/sh

BUILD_DIR=.build
ENV_PATH=$HOME/opt/biceps-env
mkdir -p $BUILD_DIR
cd $BUILD_DIR
export CPATH=$ENV_PATH/include/libxml2
cmake .. -DCMAKE_PREFIX_PATH=$ENV_PATH -DCMAKE_INSTALL_PREFIX=$ENV_PATH -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
make preprocess_float_type_header &&  make -j8 && make test_mesh && make test_fem_2d && make test_poisson && make test_stokes && make install
