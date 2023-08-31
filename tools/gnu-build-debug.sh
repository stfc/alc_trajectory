#!/usr/bin/env bash

folder="build-gnu-debug"
rm -rf $folder && mkdir $folder && cd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=Off
make
