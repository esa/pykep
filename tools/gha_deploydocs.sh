#!/usr/bin/env bash

# Echo each command.
set -x

# Exit on error.
set -e

DEPS_DIR="$HOME/local"

# Install build and docs deps into a dedicated conda environment.
conda create -y -q -p "$DEPS_DIR" \
    c-compiler cxx-compiler cmake ninja \
    "libboost>=1.73" "fmt>=10" "heyoka>=7" "heyoka.py>=7" spdlog \
    "xtensor>=0.26" xtensor-blas pagmo-devel pygmo \
    pybind11 sgp4 spiceypy matplotlib scipy \
    python=3.13 \
    "sphinx<9" sphinx-book-theme sphinxcontrib-bibtex myst-nb graphviz

# Build and install kep3 + pykep.
conda run -p "$DEPS_DIR" cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$DEPS_DIR" \
    -DCMAKE_INSTALL_PREFIX="$DEPS_DIR" \
    -Dkep3_BUILD_TESTS=OFF \
    -Dkep3_BUILD_BENCHMARKS=OFF \
    -Dkep3_BUILD_PYTHON_BINDINGS=ON
conda run -p "$DEPS_DIR" cmake --build build --target install --parallel 2

# Build HTML docs and soft-fail on external link flakiness.
(
    cd doc
    conda run -p "$DEPS_DIR" make html
    if ! conda run -p "$DEPS_DIR" make linkcheck; then
        echo "Warning: Sphinx linkcheck failed (likely transient external outage); continuing."
    fi
)

set +e
set +x