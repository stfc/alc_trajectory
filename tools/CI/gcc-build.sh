#!/usr/bin/env bash

folder="build-gcc"
rm -rf $folder && mkdir $folder && pushd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=Off
make
