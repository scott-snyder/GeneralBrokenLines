#!/bin/bash
cd /Package
mkdir build
cd build
source ../.travis-ci.d/init_x86_64.sh
source $ROOTSYS/bin/thisroot.sh
cmake -GNinja -DROOT_DIR=$ROOTSYS -DEigen3_DIR=$Eigen3_DIR -DCMAKE_CXX_FLAGS="-fdiagnostics-color=always" .. && \
ninja && \
ninja doc && \
ninja install
