#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

CMAKE_VERSION="3.11.1"
BOOST_VERSION="1.67.0"

if [[ ${AUDI_BUILD} == *36 ]]; then
	PYTHON_DIR="cp36-cp36m"
elif [[ ${AUDI_BUILD} == *35 ]]; then
	PYTHON_DIR="cp35-cp35m"
elif [[ ${AUDI_BUILD} == *34 ]]; then
	PYTHON_DIR="cp34-cp34m"
elif [[ ${AUDI_BUILD} == *27 ]]; then
	PYTHON_DIR="cp27-cp27mu"
else
	echo "Invalid build type: ${AUDI_BUILD}"
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
mkdir build
cd build
# The include directory for py3 is X.Xm, while for py2 is X.X.
# In pykep Boost_PYTHON3_LIBRARY_RELEASE and Boost_PYTHON_LIBRARY_RELEASE are used
if [[ "${PYTHON_VERSION}" != "2.7" ]]; then
    cmake -DBUILD_MAIN=no -DBUILD_PYKEP=yes -DBUILD_SPICE=yes -DBUILD_TESTS=no -DCMAKE_BUILD_TYPE=Release -DBoost_PYTHON3_LIBRARY_RELEASE=/usr/local/lib/${BOOST_PYTHON_LIB_NAME} -DPYTHON_INCLUDE_DIR=${PATH_TO_PYTHON}/include/python${PYTHON_VERSION}m/ -DPYTHON_EXECUTABLE=${PATH_TO_PYTHON}/bin/python -DCMAKE_INSTALL_PREFIX=${PATH_TO_PYTHON} ../
else
    cmake -DBUILD_MAIN=no -DBUILD_PYKEP=yes -DBUILD_SPICE=yes -DBUILD_TESTS=no -DCMAKE_BUILD_TYPE=Release -DBoost_PYTHON_LIBRARY_RELEASE=/usr/local/lib/${BOOST_PYTHON_LIB_NAME} -DPYTHON_INCLUDE_DIR=${PATH_TO_PYTHON}/include/python${PYTHON_VERSION}/ -DPYTHON_EXECUTABLE=${PATH_TO_PYTHON}/bin/python -DCMAKE_INSTALL_PREFIX=${PATH_TO_PYTHON} ../
fi
make
make install

# Compile wheels
cd /pykep/build/wheel
mv ${PATH_TO_PYTHON}/lib/python${PYTHON_VERSION}/site-packages/pykep ./
# we create a copy of the libkeplerian_toolbox.so in /usl/lib so that the system linker can find it
cp ${PATH_TO_PYTHON}/lib/libkeplerian_toolbox.so /usr/local/lib

# The following line is needed as a workaround to the auditwheel problem KeyError = .lib
# Using and compiling a null extension module (see manylinux_wheel_setup.py)
# fixes the issue (TODO: probably better ways?)
touch dummy.cpp

# We install required dependncies (do it here, do not let pip install do it)
${PATH_TO_PYTHON}/bin/pip install numpy
${PATH_TO_PYTHON}/bin/pip wheel ./ -w wheelhouse/
# Bundle external shared libraries into the wheels (only py35 has auditwheel)
auditwheel repair wheelhouse/pykep*.whl -w ./wheelhouse2/
# Install packages (not sure what --no-index -f does, should also work without, but just in case)
${PATH_TO_PYTHON}/bin/pip install pykep --no-index -f wheelhouse2
# Test
${PATH_TO_PYTHON}/bin/python -c "import pykep; print(pykep.epoch(0))"

# Upload in PyPi
# This variable will contain something if this is a tagged build (vx.y.z), otherwise it will be empty.
export PYKEP_RELEASE_VERSION=`echo "${TRAVIS_TAG}"|grep -E 'v[0-9]+\.[0-9]+.*'|cut -c 2-`
if [[ "${PYKEP_RELEASE_VERSION}" != "" ]]; then
    echo "Release build detected, uploading to PyPi."
    ${PATH_TO_PYTHON}/bin/pip install twine
    ${PATH_TO_PYTHON}/bin/twine upload -u darioizzo wheelhouse2/pykep*.whl
fi
