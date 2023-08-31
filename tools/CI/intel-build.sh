#!/usr/bin/env bash

source /opt/intel/oneapi/setvars.sh > /dev/null

folder="build-intel"
rm -rf $folder && mkdir $folder && pushd $folder
FC=ifort cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=Off
make
