[![codecov](https://codecov.io/gh/esa/kep3/branch/main/graph/badge.svg?token=K49OWJCRB9)](https://codecov.io/gh/esa/kep3)

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://esa.github.io/kep3/">
    <img src="doc/_static/pykep_logo.png" alt="Logo" width="280">
  </a>
  <p align="center">
    A coolbox for trajectory design. 
    <br />
    <a href="https://esa.github.io/kep3/"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/esa/kep3/issues/new/choose">Report bug</a>
    ·
    <a href="https://github.com/esa/kep3/issues/new/choose">Request feature</a>
    ·
    <a href="https://github.com/esa/kep3/discussions">Discuss</a>
  </p>
</p>

kep3 is an astrodynamics "coolbox" focused on trajectory design, with a modern C++ core and Python bindings (pykep). 

## Installation

Recommended workflow:

1. Install and activate the conda development environment.

```bash
conda env create -f kep3_devel.yml
conda activate kep3_devel
```

2. Configure, compile, and install the C++ library into the active conda environment.

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
  -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake --build build --config Release --target install
```

3. Build and install the Python bindings (pykep) in the same environment.

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
  -DCMAKE_PREFIX_PATH=$CONDA_PREFIX \
  -Dkep3_BUILD_PYTHON_BINDINGS=ON
cmake --build build --config Release --target install
```

### Notes On Configure Flags

The flag -DCMAKE_EXPORT_COMPILE_COMMANDS=1 tells CMake to generate a compile_commands.json file in the build directory. This file contains the exact compiler invocation used for each translation unit and is commonly used by IDE tooling, language servers, static analyzers, and code quality tools.

Flags used in the commands above:

- -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
  Installs headers, libraries, and package config files into the currently active conda environment.

- -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
  Tells CMake where to search for dependencies first (for example Boost, fmt, heyoka, xtensor, xtensor-blas) inside the active conda environment.

Common optional project flags:

- -Dkep3_BUILD_TESTS=ON
  Builds C++ unit tests.

- -Dkep3_BUILD_BENCHMARKS=ON
  Builds benchmark executables.

- -Dkep3_BUILD_PYTHON_BINDINGS=ON
  Builds and installs the pykep Python module.

- -DKEP3_VERBOSE_CONFIGURE=ON
  Enables more detailed configure-time diagnostics.

Common optional CMake flags:

- -DCMAKE_BUILD_TYPE=Debug
  Useful for debugging builds on single-config generators (for example Ninja and Unix Makefiles).

