#!/bin/bash

# build src
mkdir build
cd build
cmake ..
sudo make install
cd ..

# build doc
cd doc/sphinx
make clean && make html
