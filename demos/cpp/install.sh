#!/bin/sh

BUILD_DIR=.build
mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=true
make && make install
cd ..
