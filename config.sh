#!/bin/bash
path=`pwd`
#cd external/scons
#python bootstrap.py
export PATH=$path/external/scons/build/test-zip/bin:$PATH
export PATH=$path/src/io:$PATH
export PATH=$path/build/examples:$PATH

#cd $path
#cd external/armadillo/
#cmake . -DCMAKE_INSTALL_PREFIX:PATH=./build
#make
#make install
export ARMA_DIR=$path/external/armadillo/build
export LD_LIBRARY_PATH=$path/external/armadillo/build/lib64:$path/lib:$LD_LIBRARY_PATH
#cd $path
