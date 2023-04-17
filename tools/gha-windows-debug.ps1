# Powershell script
# Install conda environment
conda config --set always_yes yes
conda create --name pykep cmake boost boost-cpp python=3.11 scipy matplotlib
conda activate pykep

mkdir build
cd build

cmake `
    -G "Visual Studio 16 2019" -A x64 `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\pykep `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\pykep `
    -DPYKEP_BUILD_TESTS=no `
    -DPYKEP_BUILD_SPICE=yes `
    ..

cmake --build . --config Release
cmake --build . --config Release --target install


cd ..
mkdir build_python
cd build_python

cmake `
    -G "Visual Studio 16 2019" -A x64 `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\pykep `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\pykep `
    -DPYKEP_BUILD_PYKEP=yes `
    -DPYKEP_BUILD_KEP_TOOLBOX=no `
    ..

cmake --build . --config Release
cmake --build . --config Release --target install