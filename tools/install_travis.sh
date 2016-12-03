#!/bin/bash
set -e -x

cd /pykep
echo "Environment variables passed to docker from travis VM:"
echo ${BUILD_TYPE}
echo ${PATH_TO_PYTHON}
echo ${PYTHON_VERSION}
echo ${TRAVIS_TAG}

# Compile and install boost
wget --no-check-certificate https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.bz2 > /dev/null 2>&1
tar --bzip2 -xf /pykep/boost_1_62_0.tar.bz2 > /dev/null 2>&1
cd boost_1_62_0
./bootstrap.sh > /dev/null 2>&1
# removing the wrongly detected python 2.4 (deletes 5 lines after the comment "# Python configuration" )
sed -i.bak -e '/# Python configuration/,+5d' ./project-config.jam
# Defining the correct locations for python and boost_python
if [[ "${PYTHON_VERSION}" != "2.7" ]]; then #python3
    export BOOST_PYTHON_LIB_NAME=libboost_python3.so
    echo "using python" >> project-config.jam
    echo "     : ${PYTHON_VERSION}" >> project-config.jam
    echo "     : ${PATH_TO_PYTHON}/bin/python"  >> project-config.jam
    # note the m is not there !!
    echo "     : ${PATH_TO_PYTHON}/include/python${PYTHON_VERSION}m"  >> project-config.jam
    echo "     : ${PATH_TO_PYTHON}/lib"  >> project-config.jam
    echo "     ;" >> project-config.jam
else #python2
    export BOOST_PYTHON_LIB_NAME=libboost_python.so
    echo "using python" >> project-config.jam
    echo "     : ${PYTHON_VERSION}" >> project-config.jam
    echo "     : ${PATH_TO_PYTHON}/bin/python"  >> project-config.jam
    echo "     : ${PATH_TO_PYTHON}/include/python${PYTHON_VERSION}"  >> project-config.jam
    echo "     : ${PATH_TO_PYTHON}/lib"  >> project-config.jam
    echo "     ;" >> project-config.jam
fi

# Add here the boost libraries that are needed
./b2 install cxxflags="-std=c++11" --with-python --with-serialization --with-date_time --with-test --with-system > /dev/null 2>&1

# Install cmake
cd /pykep
wget --no-check-certificate https://cmake.org/files/v3.7/cmake-3.7.0.tar.gz > /dev/null 2>&1
tar xvf /pykep/cmake-3.7.0.tar.gz > /dev/null 2>&1
cd cmake-3.7.0
./bootstrap > /dev/null 2>&1
make > /dev/null 2>&1
make install > /dev/null 2>&1

# Install and compile pykep
cd /pykep
mkdir build
cd build
# The include directory for py3 is X.Xm, while for py2 is X.X.
# In pykep Boost_PYTHON3_LIBRARY_RELEASE and Boost_PYTHON_LIBRARY_RELEASE are used
if [[ "${PYTHON_VERSION}" != "2.7" ]]; then
    cmake -DBUILD_MAIN=no -DBUILD_PYKEP=yes -DBUILD_SPICE=yes -DBUILD_TESTS=no -DCMAKE_INSTALL_PREFIX=/pykep/local -DCMAKE_BUILD_TYPE=Release -DBoost_PYTHON3_LIBRARY_RELEASE=/usr/local/lib/${BOOST_PYTHON_LIB_NAME} -DPYTHON_INCLUDE_DIR=${PATH_TO_PYTHON}/include/python${PYTHON_VERSION}m/ -DPYTHON_EXECUTABLE=${PATH_TO_PYTHON}/bin/python  ../
else
    cmake -DBUILD_MAIN=no -DBUILD_PYKEP=yes -DBUILD_SPICE=yes -DBUILD_TESTS=no -DCMAKE_INSTALL_PREFIX=/pykep/local -DCMAKE_BUILD_TYPE=Release -DBoost_PYTHON_LIBRARY_RELEASE=/usr/local/lib/${BOOST_PYTHON_LIB_NAME} -DPYTHON_INCLUDE_DIR=${PATH_TO_PYTHON}/include/python${PYTHON_VERSION}/ -DPYTHON_EXECUTABLE=${PATH_TO_PYTHON}/bin/python  ../
fi
make
make install

# Compile wheels
cd /pykep/build/wheel
cp -R /pykep/local/lib/python${PYTHON_VERSION}/site-packages/pykep ./
# The following line is needed as a workaround to the auditwheel problem KeyError = .lib
# Using and compiling a null extension module (see manylinux_wheel_setup.py)
# fixes the issue (TODO: probably better ways?)
touch dummy.cpp

# We install required dependncies (do it here, do not let pip install do it)
${PATH_TO_PYTHON}/bin/pip wheel ./ -w wheelhouse/
# Bundle external shared libraries into the wheels (only py35 has auditwheel)
/opt/python/cp35-cp35m/bin/auditwheel repair wheelhouse/pykep*.whl -w ./wheelhouse2/
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
