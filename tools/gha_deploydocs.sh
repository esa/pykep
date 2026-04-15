#!/usr/bin/env bash

# Echo each command.
set -x

# Exit on error.
set -e

# Define the conda environment prefix for docs build.
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
# Configure the project with CMake inside the docs environment.
conda run -p "$DEPS_DIR" cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$DEPS_DIR" \
    -DCMAKE_INSTALL_PREFIX="$DEPS_DIR" \
    -Dkep3_BUILD_TESTS=OFF \
    -Dkep3_BUILD_BENCHMARKS=OFF \
    -Dkep3_BUILD_PYTHON_BINDINGS=ON
# Build and install targets inside the docs environment.
conda run -p "$DEPS_DIR" cmake --build build --target install --parallel 2

# Build HTML docs and soft-fail on external link flakiness.
(
    # Enter the docs directory.
    cd doc
    # Build HTML documentation.
    conda run -p "$DEPS_DIR" make html
    # Check external links and continue if transient failures occur.
    if ! conda run -p "$DEPS_DIR" make linkcheck; then
        # Emit a warning when link checking fails.
        echo "Warning: Sphinx linkcheck failed (likely transient external outage); continuing."
    fi
)

# Disable exit-on-error at script end.
set +e
# Disable command echoing at script end.
set +x