#!/usr/bin/env bash

folder="test-gnu"
rm -rf $folder && mkdir $folder && cd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=ON 
make
ctest --output-on-failure
