#!/bin/bash
cd /Package
mkdir build
cd build
source init_x86_64.sh
source $ROOTSYS/bin/thisroot.sh
cmake -GNinja -DROOT_DIR=$ROOTSYS -DEigen3_DIR=$Eigen3_DIR .. && \
ninja && \
ninja install
