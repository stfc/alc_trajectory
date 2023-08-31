#!/usr/bin/env bash

folder="build-gnu"
rm -rf $folder && mkdir $folder && cd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=Off
make
