#!/bin/sh

NOF_CORES=$(cat /proc/cpuinfo | grep "core id" | sort | uniq | wc -l)
BUILD_DIR=.build
mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DENABLE_PYTHON=ON -DENABLE_DOCS=ON -DENABLE_TESTS=OFF
make -j${NOF_CORES} && make docs && make test && sudo make install
