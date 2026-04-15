from setuptools import setup
from setuptools import find_packages
from setuptools.dist import Distribution

NAME = "pykep"
VERSION = "@PYKEP_VERSION@"
DESCRIPTION = "A coolbox for trajectory design"
LONG_DESCRIPTION = "The pykep Python package with binary extension and trajectory optimization utilities."
URL = "https://github.com/esa/kep3"
AUTHOR = "The pykep development team"
AUTHOR_EMAIL = "dario.izzo@gmail.com"
LICENSE = "MPL-2.0"
CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Mathematics",
    "Topic :: Scientific/Engineering :: Physics",
    "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
    "Programming Language :: Python :: 3",
]
KEYWORDS = "astrodynamics trajectory-optimization orbital-mechanics"
INSTALL_REQUIRES = [
    "numpy",
    "scipy",
    "matplotlib",
    "sgp4",
    "spiceypy",
    "pygmo",
    "heyoka",
]


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    url=URL,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    keywords=KEYWORDS,
    install_requires=INSTALL_REQUIRES,
    packages=find_packages(include=["pykep", "pykep.*"]),
    include_package_data=True,
    package_data={"pykep": ["*.so", "*.so.*"]},
    distclass=BinaryDistribution,
)