#!/bin/bash
path=`pwd`
#cd external/scons
#python bootstrap.py
#export PATH=$path/external/scons/build/test-zip/bin:$PATH
#cd $path
#cd external/armadillo/
#cmake . -DCMAKE_INSTALL_PREFIX:PATH=./build
#make
#make install
export ARMA_DIR=$path/external/armadillo/build
export LD_LIBRARY_PATH=$ARMA_DIR/lib:$path/lib:$LD_LIBRARY_PATH
#cd $path

# tempfix
export LD_LIBRARY_PATH=/home/anandps/local/anaconda2/lib:$LD_LIBRARY_PATH
