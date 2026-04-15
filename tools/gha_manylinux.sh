#!/usr/bin/env bash

# Enable strict error handling and pipe checks.
set -Eeuo pipefail
# Echo each command.
set -x

# Require build type input and GitHub workspace input.
: "${KEP3_BUILD_TYPE:?KEP3_BUILD_TYPE is required}"
: "${GITHUB_WORKSPACE:?GITHUB_WORKSPACE is required}"

# Print selected build type, GitHub ref, and workspace for debugging.
echo "KEP3_BUILD_TYPE: ${KEP3_BUILD_TYPE}"
echo "GITHUB_REF: ${GITHUB_REF:-<unset>}"
echo "GITHUB_WORKSPACE: ${GITHUB_WORKSPACE}"

# Mark workspace as a trusted git directory.
git config --global --add safe.directory "${GITHUB_WORKSPACE}"

# Map build type to manylinux Python ABI folder.
case "${KEP3_BUILD_TYPE}" in
    *314t*) PYTHON_DIR="cp314-cp314t" ;;
    *314*) PYTHON_DIR="cp314-cp314" ;;
    *313*) PYTHON_DIR="cp313-cp313" ;;
    *312*) PYTHON_DIR="cp312-cp312" ;;
    *311*) PYTHON_DIR="cp311-cp311" ;;
    *)
        echo "Invalid build type '${KEP3_BUILD_TYPE}'. Supported: Python314t, Python314, Python313, Python312, Python311"
        exit 1
        ;;
esac

# Compose Python toolchain path.
PYBIN="/opt/python/${PYTHON_DIR}/bin"
# Fail if selected Python executable is missing.
if [[ ! -x "${PYBIN}/python" ]]; then
    echo "Python executable not found at ${PYBIN}/python"
    exit 1
fi

# Set default install prefix, CMake prefix path, and library path for manylinux build.
PREFIX="${PREFIX:-/usr/local}"
export CMAKE_PREFIX_PATH="${PREFIX}:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="${PREFIX}/lib64:${PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export BLA_VENDOR="${BLA_VENDOR:-OpenBLAS}"

# Keep Python heyoka.py and C++ heyoka versions intentionally aligned:
# heyoka.py 7.10.1 is built against heyoka C++ 7.10.0.
HEYOKA_PYPI_VERSION="${HEYOKA_PYPI_VERSION:-7.10.1}"
HEYOKA_CPP_GIT_REF="${HEYOKA_CPP_GIT_REF:-v7.10.0}"

# Upgrade Python packaging tools, install build dependencies, and install runtime dependencies from PyPI.
"${PYBIN}/python" -m pip install --upgrade pip setuptools wheel
"${PYBIN}/python" -m pip install cmake auditwheel build
"${PYBIN}/python" -m pip install numpy scipy matplotlib cloudpickle sgp4 spiceypy pygmo "heyoka==${HEYOKA_PYPI_VERSION}"

# Select cmake binary from Python toolchain.
cmake_bin="${PYBIN}/cmake"

# Define a helper function to run commands with time measurement when available.
run_time() {
    if command -v /usr/bin/time >/dev/null 2>&1; then
        /usr/bin/time -v "$@"
    else
        "$@"
    fi
}

# Define a helper function to clone a git repository at a specific ref, with fallback.
clone_at_ref() {
    # Capture clone source URL.
    local url="$1"
    # Capture destination directory.
    local dir="$2"
    # Capture git ref/tag/branch.
    local ref="$3"

    # Remove existing checkout if present.
    rm -rf "${dir}"
    if [[ -n "${ref}" ]]; then
        git clone --depth 1 --branch "${ref}" "${url}" "${dir}" || {
            git clone --depth 1 "${url}" "${dir}"
            (
                cd "${dir}"
                git fetch --depth 1 origin "${ref}"
                git checkout "${ref}"
            )
        }
    else
        # Clone default branch when no ref is given.
        git clone --depth 1 "${url}" "${dir}"
    fi
}

# Define a helper function to build and install a CMake-based project from source.
build_and_install_cmake_repo() {
    # Capture source directory path.
    local src_dir="$1"
    # Shift to forward remaining CMake arguments.
    shift
    # Remove previous build directory.
    rm -rf "${src_dir}/build"
    # Create fresh build directory.
    mkdir -p "${src_dir}/build"
    (
        cd "${src_dir}/build"
        # Configure project with CMake.
        "${cmake_bin}" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
            -DCMAKE_PREFIX_PATH="${PREFIX}" \
            "$@" \
            ..
        # Build and install project.
        run_time "${cmake_bin}" --build . --target install --parallel 2
    )
}

# Define install root for manually built dependencies and project build/install.
INSTALL_ROOT="/root/install"
mkdir -p "${INSTALL_ROOT}"
cd "${INSTALL_ROOT}"

# Install BLAS/LAPACK dependencies needed by xtensor-blas consumers.
yum install -y openblas-devel lapack-devel || yum install -y openblas lapack

# Header-only stack required by xtensor/xtensor-blas.
XTL_REF="${XTL_REF:-0.8.1}"
XSIMD_REF="${XSIMD_REF:-13.2.0}"
XTENSOR_REF="${XTENSOR_REF:-0.26.0}"
XTENSOR_BLAS_REF="${XTENSOR_BLAS_REF:-0.22.0}"

# Clone xtl at requested ref.
clone_at_ref https://github.com/xtensor-stack/xtl.git xtl "${XTL_REF}"
# Build and install xtl.
build_and_install_cmake_repo xtl

# Clone xsimd at requested ref.
clone_at_ref https://github.com/xtensor-stack/xsimd.git xsimd "${XSIMD_REF}"
# Build and install xsimd.
build_and_install_cmake_repo xsimd

# Clone xtensor at requested ref.
clone_at_ref https://github.com/xtensor-stack/xtensor.git xtensor "${XTENSOR_REF}"
# Build and install xtensor.
build_and_install_cmake_repo xtensor

# Clone xtensor-blas at requested ref.
clone_at_ref https://github.com/xtensor-stack/xtensor-blas.git xtensor-blas "${XTENSOR_BLAS_REF}"
# Build and install xtensor-blas.
build_and_install_cmake_repo xtensor-blas

# Build heyoka C++ package for find_package(heyoka CONFIG REQUIRED).
# Clone heyoka C++ at requested ref.
clone_at_ref https://github.com/bluescarni/heyoka.git heyoka "${HEYOKA_CPP_GIT_REF}"
# Build and install heyoka C++.
build_and_install_cmake_repo heyoka \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DHEYOKA_BUILD_TESTS=OFF \
    -DHEYOKA_BUILD_BENCHMARKS=OFF \
    -DHEYOKA_BUILD_TUTORIALS=OFF \
    -DHEYOKA_ENABLE_IPO=OFF

# Build and install kep3 + pykep.
cd "${GITHUB_WORKSPACE}"
# Remove previous build directory.
rm -rf build
# Create fresh build directory.
mkdir -p build
(
    cd build
    "${cmake_bin}" ../ \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DCMAKE_PREFIX_PATH="${PREFIX}" \
        -DBLA_VENDOR="${BLA_VENDOR}" \
        -Dkep3_BUILD_TESTS=OFF \
        -Dkep3_BUILD_BENCHMARKS=OFF \
        -Dkep3_BUILD_PYTHON_BINDINGS=ON \
        -DPython3_EXECUTABLE="${PYBIN}/python"
    run_time "${cmake_bin}" --build . --target install --parallel 2
)

# Build wheel from the installed package tree.
# Define wheel working directory.
WHEEL_DIR="${GITHUB_WORKSPACE}/build/wheel"
# Ensure CMake generated setup.py exists.
if [[ ! -f "${WHEEL_DIR}/setup.py" ]]; then
    echo "Expected wheel metadata at ${WHEEL_DIR}/setup.py was not generated by CMake"
    exit 1
fi
# Clean any previous wheel staging artifacts.
rm -rf "${WHEEL_DIR}/pykep" "${WHEEL_DIR}/build" "${WHEEL_DIR}/dist" "${WHEEL_DIR}/dist2"

# Resolve site-packages path for selected Python.
SITE_PACKAGES_DIR="$(${PYBIN}/python -c 'import site; print(site.getsitepackages()[0])')"
# Copy installed pykep package into wheel staging dir.
cp -r "${SITE_PACKAGES_DIR}/pykep" "${WHEEL_DIR}/"

(
    # Enter wheel staging directory.
    cd "${WHEEL_DIR}"
    # Build binary wheel.
    run_time "${PYBIN}/python" setup.py bdist_wheel
    # Repair wheel for manylinux compliance.
    run_time "${PYBIN}/auditwheel" repair dist/pykep*.whl -w dist2
)

# Smoke-test repaired wheel in a clean context.
# Switch to root directory for clean import context.
cd /
# Install repaired wheel.
"${PYBIN}/python" -m pip install --force-reinstall "${WHEEL_DIR}"/dist2/pykep*.whl
# Run intermodule smoke checks.
"${PYBIN}/python" - <<'PY'
import gc
import pykep as pk
import pygmo as pg
import heyoka as hy

print("pykep", pk.__version__)
print("pygmo", pg.__version__)
print("heyoka", hy.__version__)

hy.install_custom_numpy_mem_handler()
try:
    # Intermodule check: heyoka expressions (C_++ generated) 
    # are converted correctly in python.
    dyn = pk.ta.zoh_kep_dyn()
    print("intermodule checks passed")
finally:
    # performing a defensive garbage collecting to prevent CI issues
    dyn = None
    gc.collect()
    hy.remove_custom_numpy_mem_handler()
    import os
    os._exit(0)
PY

# Disable command echoing at script end.
set +x