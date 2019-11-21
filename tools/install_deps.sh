#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

if [[ "${DCGP_BUILD}" != manylinux* ]]; then
    if [[ "${TRAVIS_OS_NAME}" == "osx" ]]; then
        wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh;
    else
        wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh;
    fi
    export deps_dir=$HOME/local
    export PATH="$HOME/miniconda/bin:$PATH"
    bash miniconda.sh -b -p $HOME/miniconda
    conda config --add channels conda-forge --force

    # obake-devel is needed as far as the conda package audi does not list it as a dependency
    conda_pkgs="cmake eigen boost boost-cpp tbb tbb-devel pagmo audi symengine obake-devel"

    if [[ "${DCGP_BUILD}" == "Python37" || "${DCGP_BUILD}" == "OSXPython37" ]]; then
        conda_pkgs="$conda_pkgs python=3.7 pyaudi pygmo"
    elif [[ "${DCGP_BUILD}" == "Python27" || "${DCGP_BUILD}" == "OSXPython27" ]]; then
        conda_pkgs="$conda_pkgs python=2.7 pyaudi pygmo"
    fi

    # We create the conda environment and activate it
    conda create -q -p $deps_dir -y
    source activate $deps_dir
    conda install $conda_pkgs -y

    #if [[ "${DCGP_BUILD}" == Python* ]]; then
    #    conda install doxygen graphviz -y
    #fi
fi

set +e
set +x
