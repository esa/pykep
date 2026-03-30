[![codecov](https://codecov.io/gh/esa/pykep/branch/main/graph/badge.svg?token=K49OWJCRB9)](https://codecov.io/gh/esa/pykep)

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://esa.github.io/pykep/">
    <img src="doc/_static/pykep_logo.png" alt="Logo" width="280">
  </a>
  <p align="center">
    A coolbox for trajectory design. 
    <br />
    <a href="https://esa.github.io/pykep/"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/esa/pykep/issues/new/choose">Report bug</a>
    ·
    <a href="https://github.com/esa/pykep/issues/new/choose">Request feature</a>
    ·
    <a href="https://github.com/esa/pykep/discussions">Discuss</a>
  </p>
</p>

# 🚀 Version 3 is here!

This is the official repo for **kep3** (C++ library) and its twin **pykep** (python package) version 3, the next-generation astrodynamics toolbox. Version 3 is not just an update — it's a full reimagining of what a space trajectory coolbox can be.

The code is still under development and we will only release the conda packages when we are sure all is well coordinated. In the meantime, feel free to play around discuss API and help us debug :)

If you care about orbital mechanics, trajectory optimization, or spacecraft mission design — **pykep belongs in your toolkit**. From students and researchers to mission designers and competition participants (hello, GTOC!), pykep is the Swiss Army knife of astrodynamics.

> ⚠️ The old pykep version is no longer be actively maintained. **version 3 is the future.**

---

## What is kep3?

kep3 is a C++ library with a rich Python interface (pykep) for **space flight mechanics** research. Built from the ground up with performance, usability, and extensibility in mind, it brings together everything the astrodynamics community has been asking for tailored at scientists who want to perform cutting edge research in space flight mechanics.

Whether you're computing **Lambert arcs**, propagating **Keplerian orbits**, or designing complex **multi-gravity-assist trajectories** and designing the next generation solvers for low-thrust optimization, kep3 has you covered with clean APIs, rigorous numerics, and serious speed.

---

## Why switch from version 2?

pykep served the community faithfully for years, but the new version raises the bar across the board:

- ⚡ **Faster** — rewritten for increased performance
- 🧼 **Cleaner API** — more intuitive, consistent, and Pythonic interfaces
- 🔬 **Better numerics** — improved solvers and propagators with higher accuracy
- 📦 **Modern packaging** — easy installation, better dependency management
- 📖 **Richer documentation** — with examples, tutorials, and full API reference

---

## Installation

Recommended workflow (conda):

1. Install and activate the conda development environment.

```bash
conda env create -f kep3_devel.yml
conda activate kep3_devel
```

2. Configure, compile, and install the C++ library into the active conda environment.

```bash
mkdir build && cd build
cmake -S . -G Ninja -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
ninja install
```

3. Build and install the Python bindings (pykep) in the same environment.

```bash
cmake -S . -G Ninja -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -Dkep3_BUILD_PYTHON_BINDINGS=ON
ninja install```

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

