[![C++ CI](https://github.com/esa/pykep/actions/workflows/ci-cpp.yml/badge.svg)](https://github.com/esa/pykep/actions/workflows/ci-cpp.yml)
[![Python CI](https://github.com/esa/pykep/actions/workflows/ci-python.yml/badge.svg)](https://github.com/esa/pykep/actions/workflows/ci-python.yml)
[![Manylinux CI](https://github.com/esa/pykep/actions/workflows/ci-manylinux.yml/badge.svg)](https://github.com/esa/pykep/actions/workflows/ci-manylinux.yml)
[![PyPI version](https://img.shields.io/pypi/v/pykep)](https://pypi.org/project/pykep/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/pykep?label=conda-forge)](https://anaconda.org/conda-forge/pykep)
[![codecov](https://codecov.io/gh/esa/pykep/branch/master/graph/badge.svg)](https://codecov.io/gh/esa/pykep)

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

# Pykep version 3

This is the official repo for **kep3** (C++ library) and its twin **pykep** (python package) version 3, the next-generation astrodynamics toolbox. Version 3 is not just an update: it's a full reimagining of what a space trajectory coolbox can be.

The code is still under development and we will only release the conda and PyPi packages when we are confidnt all is well coordinated. In the meantime, feel free to play around discuss API and help us debug :)

If you care about orbital mechanics, trajectory optimization, or spacecraft mission design — **pykep belongs in your toolkit**. From students and researchers to mission designers and competition participants (hello, GTOC!), pykep is the Swiss Army knife of astrodynamics.

> ⚠️ The old pykep version is no longer actively maintained. **version 3 is the future.**

---

## What is kep3?

kep3 is a C++ library with a rich Python interface (pykep) for **space flight mechanics** research. Built from the ground up with performance, usability, and extensibility in mind, it brings together everything the astrodynamics community has been asking for tailored at scientists who want to perform cutting edge research in space flight mechanics.

Whether you're computing **Lambert arcs**, propagating **Keplerian orbits**, or designing complex **multi-gravity-assist trajectories** and designing the next generation solvers for low-thrust optimization, kep3 has you covered with clean APIs, rigorous numerics, and serious speed.

---

## Versioning policy

Starting from the v3 line, this project follows [Semantic Versioning](https://semver.org/) for releases:

- `MAJOR`: incompatible API changes.
- `MINOR`: backward-compatible new features.
- `PATCH`: backward-compatible bug fixes.

For the C++ shared library (`kep`), we track ABI compatibility with `SOVERSION`:

- `VERSION` is set to the full project version (`MAJOR.MINOR.PATCH`).
- `SOVERSION` is tied to `MAJOR`.

In practice, if `MAJOR` changes, downstream binaries should expect possible ABI incompatibilities and rebuild/relink.

---

## Installation

Choose one of the following installation paths, depending on your use case.

### 1. Conda (preferred)

Conda is the recommended option for a fully managed scientific Python stack.

```bash
conda install -c conda-forge pykep
```

Current status: conda-forge currently provides the stable v1 line. v3 packages will be published once the v3 API is stabilized.

### 2. PyPI (pip)

Use pip if you prefer Python wheels from PyPI.

```bash
pip install pykep
```

### 3. Build from source (recommended for v3 development now)

Building from source is currently the recommended path for v3 development and testing.

1. Create and activate the development environment:

```bash
conda env create -f kep3_devel.yml
conda activate kep3_devel
```

2. Configure and install `kep3` and `pykep`:

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -Dkep3_BUILD_PYTHON_BINDINGS=ON
cmake --build build --target install --parallel
```

3. Validate the installation:

```bash
python -c "import pykep; print(pykep.__version__)"
```

### Configure Flags Reference

The configure step controls where `kep3` is installed, where dependencies are discovered, and which optional targets are generated.

#### Core flags used above

```cmake
-DCMAKE_EXPORT_COMPILE_COMMANDS=1
-DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"
-DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
-Dkep3_BUILD_PYTHON_BINDINGS=ON
```

1. `-DCMAKE_EXPORT_COMPILE_COMMANDS=1`
   Produces `build/compile_commands.json`, a machine-readable compilation database used by IDEs, language servers, static analyzers, and refactoring tools.

2. `-DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"`
   Installs headers, libraries, CMake package files, and Python artifacts into the active conda environment, avoiding contamination of system paths.

3. `-DCMAKE_PREFIX_PATH="$CONDA_PREFIX"`
   Prioritizes dependency resolution from the active conda environment (for example `Boost`, `fmt`, `heyoka`, `xtensor`, `xtensor-blas`). This improves reproducibility across machines.

4. `-Dkep3_BUILD_PYTHON_BINDINGS=ON`
   Enables build and installation of the `pykep` Python extension module.

#### Optional project flags

```cmake
-Dkep3_BUILD_CPP_LIBRARY
-Dkep3_BUILD_TESTS
-Dkep3_BUILD_BENCHMARKS
-DKEP3_VERBOSE_CONFIGURE
```

1. `-Dkep3_BUILD_CPP_LIBRARY=OFF`
   Controls whether the `kep3` C++ library is built from source (default). When set to `OFF`, CMake will instead locate an already-installed `kep3` via `find_package` and report a fatal error if it is not found. This is useful when you only want to build the Python bindings (or tests/benchmarks) against a `kep3` that has already been installed.

2. `-Dkep3_BUILD_TESTS=ON`
   Builds the C++ unit-test targets.

3. `-Dkep3_BUILD_BENCHMARKS=ON`
   Builds benchmark executables under `benchmark/`.

4. `-DKEP3_VERBOSE_CONFIGURE=ON`
   Emits additional configure-time diagnostics useful for dependency and toolchain troubleshooting.

#### Common CMake build-type flag

```cmake
-DCMAKE_BUILD_TYPE=Debug
```

Use this with single-config generators (for example `Ninja` and `Unix Makefiles`) when you need debug symbols and lower optimization.

#### Example: full configure command with common optional targets

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -Dkep3_BUILD_PYTHON_BINDINGS=ON \
  -Dkep3_BUILD_TESTS=ON \
  -Dkep3_BUILD_BENCHMARKS=ON
```

