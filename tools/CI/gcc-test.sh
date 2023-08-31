#!/usr/bin/env bash

folder="test-gcc"
rm -rf $folder && mkdir $folder && pushd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=ON 
make
ctest --output-on-failure
