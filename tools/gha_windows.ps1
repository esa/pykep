# PowerShell script
$ErrorActionPreference = "Stop"

$depsDir = $env:CONDA_PREFIX
if (-not $depsDir) {
    throw "CONDA_PREFIX is not set"
}

Write-Host "CONDA_PREFIX: $depsDir"

conda config --set always_yes yes
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install cmake ninja c-compiler cxx-compiler "libboost>=1.73" "fmt>=10" "heyoka>=7" "heyoka.py>=7" spdlog "xtensor>=0.26" xtensor-blas pagmo-devel pybind11 sgp4 spiceypy matplotlib scipy "eigen>=3,<5" "nlopt<2.10.1"

# Build and install kep3 + pykep.
if (Test-Path build) { Remove-Item -Recurse -Force build }
New-Item -ItemType Directory -Path build | Out-Null
Push-Location build
cmake -G "Visual Studio 17 2022" -A x64 `
    -DCMAKE_BUILD_TYPE=Release `
    -DCMAKE_INSTALL_PREFIX="$depsDir" `
    -DCMAKE_PREFIX_PATH="$depsDir" `
    -Dkep3_BUILD_TESTS=OFF `
    -Dkep3_BUILD_BENCHMARKS=OFF `
    -Dkep3_BUILD_PYTHON_BINDINGS=ON `
    ..
cmake --build . --config Release --target INSTALL --parallel 2
Pop-Location

# Run tests outside the source checkout so Python resolves the installed package.
Push-Location $env:TEMP
python -c "import pykep; pykep.test.run_test_suite()"
Pop-Location