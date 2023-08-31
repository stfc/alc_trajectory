#!/usr/bin/env bash

folder="test-intel-debug"
rm -rf $folder && mkdir $folder && cd $folder
FC=ifort cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=ON 
make
ctest --output-on-failure  
