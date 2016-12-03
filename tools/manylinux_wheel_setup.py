from distutils.core import setup, Extension
NAME = 'PyKEP'
VERSION = '@PyKEP_VERSION@'
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

extension_module = Extension(
    'dummy',
     sources=['dummy.cpp']
)

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
    ext_modules = [extension_module],
    packages=['PyKEP', 'PyKEP.core', 'PyKEP.examples', 'PyKEP.orbit_plots', 'PyKEP.phasing', 'PyKEP.planet', 'PyKEP.sims_flanagan', 'PyKEP.trajopt', 'PyKEP.util'],
    # Include pre-compiled extension
    package_data={
                'PyKEP.core': ['_core.pyd'] + DLL_LIST,
                'PyKEP.planet': ['_planet.pyd'],
                'PyKEP.sims_flanagan': ['_sims_flanagan.pyd'],
                'PyKEP.util': ['_util.pyd']
                },
)
