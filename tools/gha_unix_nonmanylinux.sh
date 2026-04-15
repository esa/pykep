#!/usr/bin/env bash

set -x
set -e

: "${KEP3_PYTHON_VERSION:?KEP3_PYTHON_VERSION must be set}"

heyoka_spec="heyoka>=7"
heyoka_py_spec="heyoka.py>=7"
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
    # On macOS arm64, newer heyoka packages may trigger template-constraint
    # failures in callable/event headers with the current toolchain.
    heyoka_spec="heyoka>=7,<8"
    heyoka_py_spec="heyoka.py>=7,<8"
fi

# Ensure conda-forge is used consistently across all non-manylinux Unix jobs.
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -y -q \
    c-compiler cxx-compiler ninja "cmake>=3.28,<3.31" \
    "python=${KEP3_PYTHON_VERSION}" \
    numpy \
    "libboost>=1.73" "fmt>=10" "${heyoka_spec}" "${heyoka_py_spec}" spdlog \
    "xtensor>=0.26" xtensor-blas pagmo-devel \
    pybind11 pygmo sgp4 spiceypy matplotlib scipy \
    "eigen>=3,<5" "nlopt<2.10.1"

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