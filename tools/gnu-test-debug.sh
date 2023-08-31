#!/usr/bin/env bash

folder="test-gnu-debug"
rm -rf $folder && mkdir $folder && cd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=TRUE 
make
ctest --output-on-failure  
