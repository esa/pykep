#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Core deps.
sudo apt-get install wget

# Install conda+deps.
MINIFORGE_VERSION=24.11.3-1
wget "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-Linux-x86_64.sh" -O miniforge3.sh
export deps_dir=$HOME/local
export PATH="$HOME/miniforge3/bin:$PATH"
bash miniforge3.sh -b -p $HOME/miniforge3
conda env create --file=kep3_devel.yml -q -p $deps_dir
source activate $deps_dir
conda list

# Install additional packages for Python compiling and docs building (19/11/2023 sphinx 7 not working)
conda install numpy "sphinx<9" sphinx-book-theme sphinxcontrib-bibtex myst-nb matplotlib pybind11 sgp4 spiceypy 

# We build and install pykep (and the kep3 library)
mkdir build
cd build
cmake -G "Ninja" ../ -DCMAKE_INSTALL_PREFIX=$deps_dir -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_BUILD_TYPE=Release -Dkep3_BUILD_TESTS=no -Dkep3_BUILD_BENCHMARKS=no -Dkep3_BUILD_PYTHON_BINDINGS=yes
cmake --build . --target=install --config=Release -- -j 2

# We get out of build as to test the global installation
cd /
python -c "import pykep.test; pykep.test.run_test_suite()"

# Build the documentation.
cd ${GITHUB_WORKSPACE}/doc
make html linkcheck

set +e
set +x
