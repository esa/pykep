from distutils.sysconfig import get_config_vars, get_python_version
import os

py_ver = get_python_version()
py_lib_path = get_config_vars()["LIBDIR"]

available_libs = []

for file in os.listdir(py_lib_path):
        if file.find("libpython") != -1:
            if file.find(py_ver) != -1:
                available_libs.append(file)

dynamic_libs = [d for d in available_libs if d.find('dylib') != -1 or d.find('.so') != -1 or d.find('.dll') != -1]
static_libs = [d for d in available_libs if d.find('a')]

if dynamic_libs != []:
    print(py_lib_path + "/" + dynamic_libs[0])
else:
    print(py_lib_path + "/" + static_libs[0])
