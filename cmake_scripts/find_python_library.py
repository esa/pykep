from distutils.sysconfig import get_config_vars, get_python_version
import os

py_ver = get_python_version()
py_lib_path = get_config_vars()["LIBDIR"]

available_libs = []

for file in os.listdir(py_lib_path):
        if file.find("libpython") != -1:
            if file.find(py_ver) != -1:
                available_libs.append(file)

print(py_lib_path + "/" + available_libs[0])
