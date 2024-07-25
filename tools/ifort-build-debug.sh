#!/usr/bin/env bash

folder="build-ifort-debug"
rm -rf $folder && mkdir $folder && cd $folder
FC=ifort cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=Off
make
