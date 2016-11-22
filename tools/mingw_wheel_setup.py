from setuptools import setup
from setuptools.dist import Distribution
from distutils import util
import sys

NAME = 'pykep'
VERSION = '@pykep_VERSION@'
DESCRIPTION = 'Basic space flight mechanics computations mostly based on perturbed Keplerian dynamics'
LONG_DESCRIPTION = 'PyKEP is a scientific library providing basic space flight mechanics computations mostly based on perturbed Keplerian dynamics.'
URL = 'https://github.com/esa/pykep'
AUTHOR = 'Dario Izzo'
AUTHOR_EMAIL = 'dario.izzo@gmail.com'
LICENSE = 'GPLv3+/LGPL3+'
CLASSIFIERS = [
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 4 - Beta',

    'Operating System :: OS Independent',

    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Astronomy',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics',

    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',

    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3'
]
KEYWORDS = 'space keplerian math physics interplanetary'
INSTALL_REQUIRES = ['numpy']
PLATFORMS = ['Unix','Windows','OSX']

class BinaryDistribution(Distribution):
    def has_ext_modules(foo):
        return True

# Setup the list of external dlls.
import os.path
mingw_wheel_libs = 'mingw_wheel_libs_python{}.txt'.format(sys.version_info[0])
l = open(mingw_wheel_libs,'r').readlines()
DLL_LIST = [os.path.basename(_[:-1]) for _ in l]

setup(name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    url=URL,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    keywords=KEYWORDS,
    platforms=PLATFORMS,
    install_requires=INSTALL_REQUIRES,
    packages=['pykep', 'pykep/core'],
    package_dir = {
            'pykep': 'pykep',
            'pykep/core': 'pykep/core',
            'pykep/examples': 'pykep/examples',
            'pykep/orbit_plots': 'pykep/orbit_plots',
            'pykep/phasing': 'pykep/phasing',
            'pykep/planet': 'pykep/planet',
            'pykep/sims_flanagan': 'pykep/sims_flanagan',
            'pykep/trajopt': 'pykep/trajopt',
            'pykep/util': 'pykep/util'
            },
    # Include pre-compiled extension
    package_data={'pykep': ['core/_core.pyd', 'planet/_planet.pyd', 'sims_flanagan/_sims_flanagan.pyd', 'util/_util.pyd'] + DLL_LIST},
    distclass=BinaryDistribution)
