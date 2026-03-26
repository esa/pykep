#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Install conda+deps with a pinned Miniforge release and host-arch installer.
MINIFORGE_VERSION=24.11.3-1
case "$(uname -m)" in
  arm64)
    MINIFORGE_ARCH=arm64
    ;;
  x86_64)
    MINIFORGE_ARCH=x86_64
    ;;
  *)
    echo "Unsupported macOS architecture: $(uname -m)" >&2
    exit 1
    ;;
esac

wget "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-MacOSX-${MINIFORGE_ARCH}.sh" -O miniforge3.sh
export deps_dir=$HOME/local
export PATH="$HOME/miniforge3/bin:$PATH"
bash miniforge3.sh -b -p "$HOME/miniforge3"
conda env create --file=kep3_devel.yml -q -p "$deps_dir"
source activate "$deps_dir"
conda list

# Install additional packages for Python compiling and docs building.
conda install numpy "sphinx<9" sphinx-book-theme sphinxcontrib-bibtex myst-nb matplotlib pybind11 sgp4 spiceypy -y

# Build and install pykep (and the kep3 library).
mkdir build
cd build
cmake -G "Ninja" ../ -DCMAKE_INSTALL_PREFIX="$deps_dir" -DCMAKE_PREFIX_PATH="$deps_dir" -DCMAKE_BUILD_TYPE=Release -Dkep3_BUILD_TESTS=no -Dkep3_BUILD_BENCHMARKS=no -Dkep3_BUILD_PYTHON_BINDINGS=yes
cmake --build . --target=install --config=Release -- -j 2

# Test the global installation.
cd /
python -c "import pykep.test; pykep.test.run_test_suite()"

set +e
set +x
