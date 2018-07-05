#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

CMAKE_VERSION="3.11.1"
BOOST_VERSION="1.67.0"

if [[ ${PYKEP_BUILD} == *37 ]]; then
	PYTHON_DIR="cp37-cp37m"
elif [[ ${PYKEP_BUILD} == *36 ]]; then
	PYTHON_DIR="cp36-cp36m"
elif [[ ${PYKEP_BUILD} == *35 ]]; then
	PYTHON_DIR="cp35-cp35m"
elif [[ ${PYKEP_BUILD} == *34 ]]; then
	PYTHON_DIR="cp34-cp34m"
elif [[ ${PYKEP_BUILD} == *27 ]]; then
	PYTHON_DIR="cp27-cp27mu"
else
	echo "Invalid build type: ${PYKEP_BUILD}"
	exit 1
fi

# HACK: for python 3.x, the include directory
# is called 'python3.xm' rather than just 'python3.x'.
# This confuses the build system of Boost.Python, thus
# we create a symlink to 'python3.x'.
cd /opt/python/${PYTHON_DIR}/include
PY_INCLUDE_DIR_NAME=`ls`
# If the include dir ends with 'm', create a symlink
# without the 'm'.
if [[ $PY_INCLUDE_DIR_NAME == *m ]]; then
	ln -s $PY_INCLUDE_DIR_NAME `echo $PY_INCLUDE_DIR_NAME|sed 's/.$//'`
fi

cd
mkdir install
cd install

# Install Boost
curl -L http://dl.bintray.com/boostorg/release/${BOOST_VERSION}/source/boost_`echo ${BOOST_VERSION}|tr "." "_"`.tar.bz2 > boost_`echo ${BOOST_VERSION}|tr "." "_"`.tar.bz2
tar xjf boost_`echo ${BOOST_VERSION}|tr "." "_"`.tar.bz2
cd boost_`echo ${BOOST_VERSION}|tr "." "_"`
sh bootstrap.sh --with-python=/opt/python/${PYTHON_DIR}/bin/python > /dev/null
./bjam --toolset=gcc link=shared threading=multi cxxflags="-std=c++11" variant=release --with-python --with-serialization --with-date_time --with-test --with-system -j2 install > /dev/null
cd ..

# Install CMake
curl -L https://github.com/Kitware/CMake/archive/v${CMAKE_VERSION}.tar.gz > v${CMAKE_VERSION}
tar xzf v${CMAKE_VERSION} > /dev/null 2>&1
cd CMake-${CMAKE_VERSION}/
./configure > /dev/null
gmake -j2 > /dev/null
gmake install > /dev/null
cd ..

# Install and compile pykep
cd /pykep
cd build
cmake -DBUILD_MAIN=no -DBUILD_PYKEP=yes -DBUILD_SPICE=yes -DBUILD_TESTS=no -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=/opt/python/${PYTHON_DIR}/bin/python ../;
make -j2 install

# Compile wheels
#cd /pykep/build/wheel
#mv ${PATH_TO_PYTHON}/lib/python${PYTHON_VERSION}/site-packages/pykep ./
# we create a copy of the libkeplerian_toolbox.so in /usr/lib so that the system linker can find it
# cp ${PATH_TO_PYTHON}/lib/libkeplerian_toolbox.so /usr/local/lib


# Compile wheels
cd wheel
# Copy the installed pyaudi files, wherever they might be in /usr/local,
# into the current dir.
cp -a `find /usr/local/lib -type d -iname 'pykep'` ./

/opt/python/${PYTHON_DIR}/bin/pip install numpy wheel

# Create the wheel and repair it.
/opt/python/${PYTHON_DIR}/bin/python setup.py bdist_wheel
auditwheel repair dist/pykep* -w ./dist2
# Try to install it and run the tests.
cd /
/opt/python/${PYTHON_DIR}/bin/pip install /pykep/build/wheel/dist2/pykep*
/opt/python/${PYTHON_DIR}/bin/python -c "import pykep; print(pykep.epoch(0))"

# Upload in PyPi
# This variable will contain something if this is a tagged build (vx.y.z), otherwise it will be empty.
export PYKEP_RELEASE_VERSION=`echo "${TRAVIS_TAG}"|grep -E 'v[0-9]+\.[0-9]+.*'|cut -c 2-`
if [[ "${PYKEP_RELEASE_VERSION}" != "" ]]; then
    echo "Release build detected, uploading to PyPi."
    /opt/python/${PYTHON_DIR}/bin/pip install twine
    /opt/python/${PYTHON_DIR}/bin/twine upload -u darioizzo /audi/build/wheel/dist2/pykep*
fi


