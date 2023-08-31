#!/usr/bin/env bash

folder="test-gcc-debug"
rm -rf $folder && mkdir $folder && pushd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=ON 
make 
ctest --output-on-failure
