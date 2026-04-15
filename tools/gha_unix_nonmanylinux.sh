#!/usr/bin/env bash

# Echo each command.
set -x
# Exit on first error.
set -e

# Require the Python version input from the workflow.
: "${KEP3_PYTHON_VERSION:?KEP3_PYTHON_VERSION must be set}"

# Ensure conda-forge is used consistently across all non-manylinux Unix jobs.
# Add conda-forge to the channel list.
conda config --add channels conda-forge
# Enforce strict channel priority.
conda config --set channel_priority strict
# Install build and runtime dependencies in the active env.
conda install -y -q \
    c-compiler cxx-compiler ninja "cmake>=3.28,<3.31" \
    "python=${KEP3_PYTHON_VERSION}" \
    "libboost>=1.73" "fmt>=10" "heyoka>=7" "heyoka.py>=7" spdlog \
    "xtensor>=0.26" xtensor-blas pagmo-devel pygmo \
    pybind11 sgp4 spiceypy matplotlib scipy \
    "eigen>=3,<5" "nlopt<2.10.1"

# Cache the active conda prefix for build/install paths.
deps_dir="${CONDA_PREFIX}"

# On macOS arm64 force the conda clang toolchain to avoid AppleClang/header incompatibilities.
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
    # Ensure the conda clang executable exists.
    if [[ ! -x "${deps_dir}/bin/clang" || ! -x "${deps_dir}/bin/clang++" ]]; then
        # Report missing toolchain and stop.
        echo "Expected conda clang toolchain not found in ${deps_dir}/bin"
        exit 1
    fi
    # Force C compiler to conda clang.
    export CC="${deps_dir}/bin/clang"
    # Force C++ compiler to conda clang++.
    export CXX="${deps_dir}/bin/clang++"
fi

# Build and install kep3 + pykep.
# Remove any previous build directory.
rm -rf build
# Create a fresh build directory.
mkdir build
# Enter the build directory.
cd build
# Configure the project with CMake.
cmake ../ \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${deps_dir}" \
    -DCMAKE_PREFIX_PATH="${deps_dir}" \
    -Dkep3_BUILD_TESTS=OFF \
    -Dkep3_BUILD_BENCHMARKS=OFF \
    -Dkep3_BUILD_PYTHON_BINDINGS=ON
# Build and install targets.
cmake --build . --target install --parallel 2
# Return to repository root.
cd ..

# Run the Python test suite from outside the source tree so the installed package is imported.
(
    # Run tests from home so imports resolve to installed package.
    cd "${HOME}"
    # Execute the pykep Python test suite.
    python -c "import pykep; pykep.test.run_test_suite()"
)

# Disable exit-on-error at script end.
set +e
# Disable command echoing at script end.
set +x