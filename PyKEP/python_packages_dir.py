import distutils.sysconfig
import os

print os.path.split(distutils.sysconfig.get_python_lib())[-1]
