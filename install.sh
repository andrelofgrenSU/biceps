#!/bin/sh

NOF_CORES=$(cat /proc/cpuinfo | grep "core id" | sort | uniq | wc -l)
BUILD_DIR=.build
mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DENABLE_PYTHON=ON -DENABLE_TESTS=ON -DENABLE_DOCS=OFF
# make -j${NOF_CORES} && make docs && make test && sudo make install
make -j${NOF_CORES} && make test && sudo make install
cd ..
