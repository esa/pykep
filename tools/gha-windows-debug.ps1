# Powershell script
# Install conda environment
conda config --set always_yes yes
conda create --name pykep cmake boost boost-cpp python=3.11 scipy matplotlib
conda activate pykep

mkdir build
cd build

cmake `
    -G "Visual Studio 16 2019" `
    -A x64 `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\pykep `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\pykep `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DPYKEP_BUILD_KEP_TOOLBOX=yes `
    -DPYKEP_BUILD_PYKEP=no `
    -DPYKEP_BUILD_SPICE=yes `
    -DPYKEP_BUILD_TESTS=yes `
    ..

cmake --build . --config Debug --target install
ctest -VV --output-on-failure -j4 -C Debug
