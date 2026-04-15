#!/usr/bin/env bash

set -x
set -e

: "${KEP3_PYTHON_VERSION:?KEP3_PYTHON_VERSION must be set}"

# Ensure conda-forge is used consistently across all non-manylinux Unix jobs.
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -y -q \
    c-compiler cxx-compiler ninja "cmake>=3.28,<3.31" \
    "python=${KEP3_PYTHON_VERSION}" \
    numpy \
    "libboost>=1.73" "fmt>=10" "heyoka>=7" "heyoka.py>=7" spdlog \
    "xtensor>=0.26" xtensor-blas pagmo-devel \
    pybind11 pygmo sgp4 spiceypy matplotlib scipy \
    "eigen>=3,<5" "nlopt<2.10.1"

# Force CMake to use the exact compiler wrappers selected by the activated conda
# toolchain. This avoids accidental fallback to system compilers (e.g. Xcode
# clang on macOS) which can be ABI- or feature-incompatible with conda-forge
# dependencies. The prefix guard makes such mismatches fail fast.
: "${CC:?CC must be set by the active conda toolchain}"
: "${CXX:?CXX must be set by the active conda toolchain}"
if [[ "${CC}" != "${CONDA_PREFIX}/bin/"* ]] || [[ "${CXX}" != "${CONDA_PREFIX}/bin/"* ]]; then
    echo "ERROR: expected CC/CXX from ${CONDA_PREFIX}/bin, got CC=${CC} CXX=${CXX}" >&2
    exit 1
fi

cmake_compiler_args=("-DCMAKE_C_COMPILER=${CC}" "-DCMAKE_CXX_COMPILER=${CXX}")

deps_dir="${CONDA_PREFIX}"

# Build and install kep3 + pykep.
rm -rf build
mkdir build
cd build
cmake ../ \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${deps_dir}" \
    -DCMAKE_PREFIX_PATH="${deps_dir}" \
    "${cmake_compiler_args[@]}" \
    -Dkep3_BUILD_TESTS=OFF \
    -Dkep3_BUILD_BENCHMARKS=OFF \
    -Dkep3_BUILD_PYTHON_BINDINGS=ON
cmake --build . --target install --parallel 2
cd ..

# Run the Python test suite from outside the source tree so the installed package is imported.
(
    cd "${HOME}"
    python -c "import pykep; pykep.test.run_test_suite()"
)

set +e
set +x