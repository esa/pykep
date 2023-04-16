# Powershell script
# Install conda environment
conda config --set always_yes yes
conda create --name pykep cmake boost boost-cpp python=3.11 scipy matplotlib
conda activate pykep

# Define environment variables for clang ...
# ... and make them persistent 
cmd.exe /c "call `"C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvars64.bat`" && set > %temp%\vcvars.txt"

Get-Content "$env:temp\vcvars.txt" | Foreach-Object {
  if ($_ -match "^(.*?)=(.*)$") {
    Set-Content "env:\$($matches[1])" $matches[2]
  }
}

mkdir build
cd build

cmake `
    -G "Ninja" `
    -DCMAKE_C_COMPILER=clang-cl `
    -DCMAKE_CXX_COMPILER=clang-cl `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\pykep `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\pykep `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DPYKEP_BUILD_KEP_TOOLBOX=yes `
    -DPYKEP_BUILD_PYKEP=no `
    -DPYKEP_BUILD_SPICE=yes `
    -DPYKEP_BUILD_TESTS=yes `
    ..

cmake --build . --target install --config Release
ctest -j4 -V -C Release
