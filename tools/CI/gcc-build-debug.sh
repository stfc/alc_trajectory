#!/usr/bin/env bash

folder="build-gcc-debug"
rm -rf $folder && mkdir $folder && pushd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=Off
make
