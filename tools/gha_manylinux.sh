#!/usr/bin/env bash

set -Eeuo pipefail
set -x

: "${KEP3_BUILD_TYPE:?KEP3_BUILD_TYPE is required}"
: "${GITHUB_WORKSPACE:?GITHUB_WORKSPACE is required}"

echo "KEP3_BUILD_TYPE: ${KEP3_BUILD_TYPE}"
echo "GITHUB_REF: ${GITHUB_REF:-<unset>}"
echo "GITHUB_WORKSPACE: ${GITHUB_WORKSPACE}"

git config --global --add safe.directory "${GITHUB_WORKSPACE}"

case "${KEP3_BUILD_TYPE}" in
    *314*) PYTHON_DIR="cp314-cp314" ;;
    *313*) PYTHON_DIR="cp313-cp313" ;;
    *312*) PYTHON_DIR="cp312-cp312" ;;
    *311*) PYTHON_DIR="cp311-cp311" ;;
    *)
        echo "Invalid build type '${KEP3_BUILD_TYPE}'. Supported: Python314, Python313, Python312, Python311"
        exit 1
        ;;
esac

PYBIN="/opt/python/${PYTHON_DIR}/bin"
if [[ ! -x "${PYBIN}/python" ]]; then
    echo "Python executable not found at ${PYBIN}/python"
    exit 1
fi

PREFIX="${PREFIX:-/usr/local}"
export CMAKE_PREFIX_PATH="${PREFIX}:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="${PREFIX}/lib64:${PREFIX}/lib:${LD_LIBRARY_PATH:-}"

"${PYBIN}/python" -m pip install --upgrade pip setuptools wheel
"${PYBIN}/python" -m pip install cmake ninja auditwheel build
"${PYBIN}/python" -m pip install numpy scipy matplotlib cloudpickle sgp4 spiceypy pygmo heyoka

cmake_bin="${PYBIN}/cmake"

if [[ ! -x "${cmake_bin}" ]]; then
    echo "CMake from pip was not found at ${cmake_bin}"
    exit 1
fi

if ! command -v git >/dev/null 2>&1; then
    echo "git is required in the container image"
    exit 1
fi

run_time() {
    if command -v /usr/bin/time >/dev/null 2>&1; then
        /usr/bin/time -v "$@"
    else
        "$@"
    fi
}

clone_at_ref() {
    local url="$1"
    local dir="$2"
    local ref="$3"

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
        git clone --depth 1 "${url}" "${dir}"
    fi
}

build_and_install_cmake_repo() {
    local src_dir="$1"
    shift
    rm -rf "${src_dir}/build"
    mkdir -p "${src_dir}/build"
    (
        cd "${src_dir}/build"
        "${cmake_bin}" -G Ninja \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
            -DCMAKE_PREFIX_PATH="${PREFIX}" \
            "$@" \
            ..
        run_time "${cmake_bin}" --build . --target install --parallel 2
    )
}

INSTALL_ROOT="/root/install"
mkdir -p "${INSTALL_ROOT}"
cd "${INSTALL_ROOT}"

# Header-only stack required by xtensor/xtensor-blas.
XTL_REF="${XTL_REF:-0.8.1}"
XSIMD_REF="${XSIMD_REF:-13.2.0}"
XTENSOR_REF="${XTENSOR_REF:-0.26.0}"
XTENSOR_BLAS_REF="${XTENSOR_BLAS_REF:-0.22.0}"

clone_at_ref https://github.com/xtensor-stack/xtl.git xtl "${XTL_REF}"
build_and_install_cmake_repo xtl

clone_at_ref https://github.com/xtensor-stack/xsimd.git xsimd "${XSIMD_REF}"
build_and_install_cmake_repo xsimd

clone_at_ref https://github.com/xtensor-stack/xtensor.git xtensor "${XTENSOR_REF}"
build_and_install_cmake_repo xtensor

clone_at_ref https://github.com/xtensor-stack/xtensor-blas.git xtensor-blas "${XTENSOR_BLAS_REF}"
build_and_install_cmake_repo xtensor-blas

# Build heyoka C++ package for find_package(heyoka CONFIG REQUIRED).
HEYOKA_REF="${HEYOKA_REF:-}"
clone_at_ref https://github.com/bluescarni/heyoka.git heyoka "${HEYOKA_REF}"
build_and_install_cmake_repo heyoka \
    -DHEYOKA_BUILD_TESTS=OFF \
    -DHEYOKA_BUILD_BENCHMARKS=OFF \
    -DHEYOKA_BUILD_TUTORIALS=OFF \
    -DHEYOKA_BUILD_PYTHON_BINDINGS=OFF \
    -DHEYOKA_ENABLE_IPO=OFF

# Build and install kep3 + pykep.
cd "${GITHUB_WORKSPACE}"
rm -rf build
mkdir -p build
(
    cd build
    "${cmake_bin}" ../ \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DCMAKE_PREFIX_PATH="${PREFIX}" \
        -Dkep3_BUILD_TESTS=OFF \
        -Dkep3_BUILD_BENCHMARKS=OFF \
        -Dkep3_BUILD_PYTHON_BINDINGS=ON \
        -DPython3_EXECUTABLE="${PYBIN}/python"
    run_time "${cmake_bin}" --build . --target install --parallel 2
)

# Build wheel from the installed package tree.
WHEEL_DIR="${GITHUB_WORKSPACE}/build/wheel"
rm -rf "${WHEEL_DIR}"
mkdir -p "${WHEEL_DIR}"

SITE_PACKAGES_DIR="$(${PYBIN}/python -c 'import site; print(site.getsitepackages()[0])')"
cp -r "${SITE_PACKAGES_DIR}/pykep" "${WHEEL_DIR}/"
cp "${GITHUB_WORKSPACE}/tools/wheel_setup.py" "${WHEEL_DIR}/setup.py"
cp "${GITHUB_WORKSPACE}/tools/wheel_setup.cfg" "${WHEEL_DIR}/setup.cfg"

PYKEP_VERSION="$(${PYBIN}/python -c 'import pykep; print(pykep.__version__)')"
sed -i "s/@PYKEP_VERSION@/${PYKEP_VERSION}/g" "${WHEEL_DIR}/setup.py"

(
    cd "${WHEEL_DIR}"
    run_time "${PYBIN}/python" setup.py bdist_wheel
    run_time "${PYBIN}/auditwheel" repair dist/pykep*.whl -w dist2
)

# Smoke-test repaired wheel in a clean context.
cd /
"${PYBIN}/python" -m pip install --force-reinstall "${WHEEL_DIR}"/dist2/pykep*.whl
"${PYBIN}/python" -c "import pykep, pygmo, heyoka; print('pykep', pykep.__version__)"

set +x