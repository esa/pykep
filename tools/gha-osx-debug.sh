#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Core deps.
sudo apt-get install wget

# Install conda+deps.
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh -O mambaforge.sh
export deps_dir=$HOME/local
export PATH="$HOME/mambaforge/bin:$PATH"
bash mambaforge.sh -b -p $HOME/mambaforge
mamba create -y -q -p $deps_dir c-compiler cxx-compiler cmake boost boost-cpp python=3.11 scipy matplotlib
source activate $deps_dir

# Create the build dir and cd into it.
mkdir build

# Install keplerian_toolbox
cd build
cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=Debug \
    -DPYKEP_BUILD_KEP_TOOLBOX=yes \
    -DPYKEP_BUILD_PYKEP=no \
    -DPYKEP_BUILD_SPICE=yes \
    -DPYKEP_BUILD_TESTS=yes \
    ..
make VERBOSE=1 install
ctest -j4 -V
cd ..

# Install pykep 
mkdir build_pykep
cd build_pykep
cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=Debug \
    -DPYKEP_BUILD_KEP_TOOLBOX=no \
    -DPYKEP_BUILD_PYKEP=yes \
    ..

make VERBOSE=1 install
python -c "import pykep.test; pykep.test.run_test_suite()"
set +e
set +x