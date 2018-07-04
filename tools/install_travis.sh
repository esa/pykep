#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

if [[ "${PYKEP_BUILD}" == manylinux* ]]; then
    cd ..;
    docker pull ${DOCKER_IMAGE};
    docker run --rm -e TWINE_PASSWORD -e PYKEP_BUILD -e TRAVIS_TAG -v `pwd`:/pykep $DOCKER_IMAGE bash /pykep/tools/install_docker.sh
fi

set +e
set +x
