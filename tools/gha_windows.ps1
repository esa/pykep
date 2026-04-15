# PowerShell script
# Stop execution on errors.
$ErrorActionPreference = "Stop"

# Read active conda environment path.
$depsDir = $env:CONDA_PREFIX
# Fail if no conda environment is active.
if (-not $depsDir) {
    throw "CONDA_PREFIX is not set"
}

# Print active conda prefix for debugging.
Write-Host "CONDA_PREFIX: $depsDir"

# Enable non-interactive conda operations.
conda config --set always_yes yes
# Ensure conda-forge channel is configured.
conda config --add channels conda-forge
# Enforce strict channel priority.
conda config --set channel_priority strict
# Install build and runtime dependencies.
conda install cmake ninja c-compiler cxx-compiler "libboost>=1.73" "fmt>=10" "heyoka>=7" "heyoka.py>=7" spdlog "xtensor>=0.26" xtensor-blas pagmo-devel pygmo pybind11 sgp4 spiceypy matplotlib scipy "eigen>=3,<5" "nlopt<2.10.1"

# Build and install kep3 + pykep.
# Remove previous build directory if it exists.
if (Test-Path build) { Remove-Item -Recurse -Force build }
# Create a fresh build directory.
New-Item -ItemType Directory -Path build | Out-Null
# Enter the build directory.
Push-Location build
# Configure the project with CMake.
cmake -G "Visual Studio 17 2022" -A x64 `
    -DCMAKE_BUILD_TYPE=Release `
    -DCMAKE_INSTALL_PREFIX="$depsDir" `
    -DCMAKE_PREFIX_PATH="$depsDir" `
    -Dkep3_BUILD_TESTS=OFF `
    -Dkep3_BUILD_BENCHMARKS=OFF `
    -Dkep3_BUILD_PYTHON_BINDINGS=ON `
    ..
# Build and install targets.
cmake --build . --config Release --target INSTALL --parallel 2
# Return to previous directory.
Pop-Location

# Run tests outside the source checkout so Python resolves the installed package.
# Switch to temp directory to avoid local import shadowing.
Push-Location $env:TEMP
# Run the installed pykep test suite.
python -c "import pykep; pykep.test.run_test_suite()"
# Return to previous directory.
Pop-Location