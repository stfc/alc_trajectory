#!/usr/bin/env bash

folder="test-ifx-debug"
rm -rf $folder && mkdir $folder && cd $folder
FC=ifx cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=ON 
make
ctest --output-on-failure  
