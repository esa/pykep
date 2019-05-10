import os
import re
import sys


def wget(url, out):
    import urllib.request
    print('Downloading "' + url + '" as "' + out + '"')
    urllib.request.urlretrieve(url, out)


def rm_fr(path):
    import shutil
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def run_command(raw_command, directory=None, verbose=True):
    # Helper function to run a command and display optionally its output
    # unbuffered.
    import shlex
    from subprocess import Popen, PIPE, STDOUT
    print(raw_command)
    proc = Popen(shlex.split(raw_command), cwd=directory,
                 stdout=PIPE, stderr=STDOUT)
    if verbose:
        output = ''
        while True:
            line = proc.stdout.readline()
            if not line:
                break
            line = str(line, 'utf-8')
            # Don't print the newline character.
            print(line[:-1])
            sys.stdout.flush()
            output += line
        proc.communicate()
    else:
        output = str(proc.communicate()[0], 'utf-8')
    if proc.returncode:
        raise RuntimeError(output)
    return output


# Build type setup.
BUILD_TYPE = os.environ['BUILD_TYPE']
is_release_build = (os.environ['APPVEYOR_REPO_TAG'] == 'true') and bool(
    re.match(r'v[0-9]+\.[0-9]+.*', os.environ['APPVEYOR_REPO_TAG_NAME']))
if is_release_build:
    print("Release build detected, tag is '" +
          os.environ['APPVEYOR_REPO_TAG_NAME'] + "'")
is_python_build = 'Python' in BUILD_TYPE

# Get mingw and set the path.
wget(r'https://github.com/bluescarni/binary_deps/raw/master/x86_64-8.1.0-release-posix-seh-rt_v6-rev0.7z', 'mw64.7z')
run_command(r'7z x -oC:\\ mw64.7z', verbose=False)
ORIGINAL_PATH = os.environ['PATH']
os.environ['PATH'] = r'C:\\mingw64\\bin;' + os.environ['PATH']

# Download boost (this includes also all the boost_python libraries)
wget(r'https://github.com/bluescarni/binary_deps/raw/master/boost_mgw81-mt-x64-1_70.7z', 'boost.7z')
# Extract them.
run_command(r'7z x -aoa -oC:\\ boost.7z', verbose=False)

# Setup of the dependencies for a Python build.
if is_python_build:
    if 'Python37' in BUILD_TYPE:
        python_version = '37'
    elif 'Python36' in BUILD_TYPE:
        python_version = '36'
    elif 'Python27' in BUILD_TYPE:
        python_version = '27'
    else:
        raise RuntimeError('Unsupported Python build: ' + BUILD_TYPE)

    python_package = r'python' + python_version + r'_mingw_64.7z'
    # Remove any existing Python installation named PythonXX in c:\
    rm_fr(r'c:\\Python' + python_version)
    # Set paths.
    pinterp = r'c:\\Python' + python_version + r'\\python.exe'
    pip = r'c:\\Python' + python_version + r'\\scripts\\pip'
    twine = r'c:\\Python' + python_version + r'\\scripts\\twine'
    pykep_install_path = r'C:\\Python' + \
        python_version + r'\\Lib\\site-packages\\pykep'
    # Get Python.
    wget(r'https://github.com/bluescarni/binary_deps/raw/master/' +
         python_package, 'python.7z')
    run_command(r'7z x -aoa -oC:\\ python.7z', verbose=False)
      # Install pip and deps.
    wget(r'https://bootstrap.pypa.io/get-pip.py', 'get-pip.py')
    run_command(pinterp + ' get-pip.py --force-reinstall')
    if is_release_build:
        run_command(pip + ' install twine')

# Set the path so that the precompiled libs can be found.
os.environ['PATH'] = os.environ['PATH'] + r';c:\\local\\lib'

# Proceed to the build.
os.makedirs('build')
os.chdir('build')
common_cmake_opts = r'-DCMAKE_PREFIX_PATH=c:\\local -DCMAKE_INSTALL_PREFIX=c:\\local -DBUILD_SPICE=yes'

# Configuration step.
if is_python_build:
    # WARNING: Cannot break the following in multiple lines using \ 
    run_command(r'cmake -G "MinGW Makefiles" .. -DCMAKE_BUILD_TYPE=Release ' + common_cmake_opts + r' -DBUILD_PYKEP=yes -DBUILD_TESTS=no -DBoost_INCLUDE_DIR=c:\\local\\include -DBoost_SERIALIZATION_LIBRARY_RELEASE=c:\\local\\lib\\libboost_serialization-mgw81-mt-x64-1_70.dll -DBoost_DATE_TIME_LIBRARY_RELEASE=c:\\local\\lib\\libboost_date_time-mgw81-mt-x64-1_70.dll -DBoost_PYTHON' + python_version + r'_LIBRARY_RELEASE=c:\\local\\lib\\libboost_python' + python_version + r'-mgw81-mt-x64-1_70.dll -DPYTHON_EXECUTABLE=C:\\Python' + python_version + r'\\python.exe -DPYTHON_LIBRARY=C:\\Python' + python_version + r'\\python' + python_version + r'.dll')
elif BUILD_TYPE in ['Release', 'Debug']:
    cmake_opts = r'-DCMAKE_BUILD_TYPE=' + BUILD_TYPE + \
        r' -DBUILD_TESTS=yes ' + common_cmake_opts
    run_command(r'cmake -G "MinGW Makefiles" .. ' + cmake_opts)
else:
    raise RuntimeError('Unsupported build type: ' + BUILD_TYPE)

# Build+install step.
#run_command(r'cmake --build . --target install')
run_command(r'mingw32-make install VERBOSE=1')

# Move the shared library in the correct python DESTINATION
#run_command(r'move "C:\Program Files (x86)\Keplerian_Toolbox\lib\libkeplerian_toolbox.dll" "c:\local\lib"')

# Testing, packaging.
if is_python_build:
    # Run the Python tests.
    run_command(
        pinterp + r' -c "import pykep; print(pykep.epoch(0))"')
    # Build the wheel.
    import shutil
    os.chdir('wheel')
    shutil.move(pykep_install_path, r'.')
    wheel_libs = 'mingw_wheel_libs_python{}.txt'.format(python_version)
    DLL_LIST = [_[:-1] for _ in open(wheel_libs, 'r').readlines()]
    for _ in DLL_LIST:
        shutil.copy(_, 'pykep/core')
    run_command(pinterp + r' setup.py bdist_wheel')
    os.environ['PATH'] = ORIGINAL_PATH
    run_command(pip + r' install dist\\' + os.listdir('dist')[0])

    os.chdir('/')
    run_command(
        pinterp + r' -c "import pykep; print(pykep.epoch(0))"')
    if is_release_build:
        os.chdir('C:/projects/pykep/build/wheel')
        run_command(twine + r' upload -u darioizzo dist\\' +
                    os.listdir('dist')[0])
elif BUILD_TYPE == 'Release':
    run_command(r'ctest -VV')
elif BUILD_TYPE == 'Debug':
    run_command(r'ctest -VV')
else:
    raise RuntimeError('Unsupported build type: ' + BUILD_TYPE)
