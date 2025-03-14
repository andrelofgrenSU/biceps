#!/bin/sh

BUILD_DIR=.build
mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release -DENABLE_PYTHON=ON -DENABLE_DOCS=ON -DENABLE_TESTS=ON
make -j8 && make tests && sudo make install
