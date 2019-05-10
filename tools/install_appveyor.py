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

# ----------------------------------SCRIPT START-----------------------------------------#
# Build type setup.
BUILD_TYPE = os.environ['BUILD_TYPE']
is_release_build = (os.environ['APPVEYOR_REPO_TAG'] == 'true') and bool(
    re.match(r'v[0-9]+\.[0-9]+.*', os.environ['APPVEYOR_REPO_TAG_NAME']))
if is_release_build:
    print("Release build detected, tag is '" +
          os.environ['APPVEYOR_REPO_TAG_NAME'] + "'")
is_python_build = 'Python' in BUILD_TYPE

# Check here for a list of installed software in the appveyor VMs: https://www.appveyor.com/docs/windows-images-software/
# USING: mingw64 8.1.0
ORIGINAL_PATH = os.environ['PATH']
os.environ['PATH'] = r'C:\\mingw-w64\\x86_64-8.1.0-posix-seh-rt_v6-rev0\\mingw64\\bin;' + os.environ['PATH']
# Set the path so that the precompiled boost libs can be found.
os.environ['PATH'] = os.environ['PATH'] + r';c:\\local\\lib'
# Download boost (this includes also all the boost_python libraries, will all be installed in \\local)
# USING: boost 1.70 from binary deps
wget(r'https://github.com/bluescarni/binary_deps/raw/master/boost_mgw81-mt-x64-1_70.7z', 'boost.7z')
# Extract them.
run_command(r'7z x -aoa -oC:\\ boost.7z', verbose=False)

# Setup of the Python build variables (version based)
if is_python_build:
    if 'Python37-x64' in BUILD_TYPE:
        python_version = r'37'
        python_folder = r'Python37-x64'
        python_library = r'C:\\' + python_folder + r'\\python37.dll '
    elif 'Python36-x64' in BUILD_TYPE:
        python_version = '36'
        python_folder = r'Python36-x64'
        python_library = r'C:\\' + python_folder + r'\\python36.dll '
    elif 'Python27-x64' in BUILD_TYPE:
        python_version = r'27'
        python_folder = r'Python27-x64'
        python_library = r'C:\\' + python_folder + r'\\libs\\python27.dll '
        # Fot py27 I could not get it to work with the appveyor python (I was close but got tired).
        # Since this is anyway going to disappear (py27 really!!!), I am handling it as a one time workaround using the old py27 patched by bluescarni
        rm_fr(r'c:\\Python27-x64')
        wget(r'https://github.com/bluescarni/binary_deps/raw/master/python27_mingw_64.7z', 'python.7z')
        run_command(r'7z x -aoa -oC:\\ python.7z', verbose=False)
        run_command(r'mv C:\\Python27 C:\\Python27-x64', verbose=False)
    else:
        raise RuntimeError('Unsupported Python build: ' + BUILD_TYPE)

    # Set paths.
    pinterp = r"C:\\" + python_folder + r'\\python.exe'
    pip = r"C:\\" + python_folder + r'\\scripts\\pip'
    twine = r"C:\\" + python_folder + r'\\scripts\\twine'
    module_install_path = r"C:\\" + python_folder + r'\\Lib\\site-packages\\pykep'
    # Install pip and deps.
    run_command(pinterp + r' --version', verbose=True)
    wget(r'https://bootstrap.pypa.io/get-pip.py', 'get-pip.py')
    run_command(pinterp + ' get-pip.py --force-reinstall')
    if is_release_build:
        run_command(pip + ' install twine')

# Proceed to the build.
os.makedirs('build')
os.chdir('build')
common_cmake_opts = r'-DCMAKE_PREFIX_PATH=c:\\local -DCMAKE_INSTALL_PREFIX=c:\\local -DBUILD_SPICE=no '

# Configuration step.
if is_python_build:
    run_command(r'cmake -G "MinGW Makefiles" .. -DCMAKE_BUILD_TYPE=Release ' + \
        common_cmake_opts + \
        r'-DBUILD_PYKEP=yes ' + \
        r'-DBUILD_MAIN=no ' + \
        r'-DBUILD_TESTS=no ' + \
        r'-DBoost_INCLUDE_DIR=c:\\local\\include ' + \
        r'-DBoost_SERIALIZATION_LIBRARY_RELEASE=c:\\local\\lib\\libboost_serialization-mgw81-mt-x64-1_70.dll ' + \
        r'-DBoost_DATE_TIME_LIBRARY_RELEASE=c:\\local\\lib\\libboost_date_time-mgw81-mt-x64-1_70.dll ' + \
        r'-DBoost_PYTHON' + python_version + r'_LIBRARY_RELEASE=c:\\local\\lib\\libboost_python' + python_version + r'-mgw81-mt-x64-1_70.dll ' + \
        r'-DPYTHON_INCLUDE_DIR=C:\\' + python_folder + r'\\include ' + \
        r'-DPYTHON_EXECUTABLE=C:\\' + python_folder + r'\\python.exe ' + \
        r'-DPYTHON_LIBRARY=' + python_library + r' ' + \
        r'-DCMAKE_CXX_FLAGS="-D_hypot=hypot"')
elif BUILD_TYPE in ['Release', 'Debug']:
    cmake_opts = r'-DCMAKE_BUILD_TYPE=' + BUILD_TYPE + \
        r' -DBUILD_TESTS=yes ' + common_cmake_opts
    run_command(r'cmake -G "MinGW Makefiles" .. ' + cmake_opts)
else:
    raise RuntimeError('Unsupported build type: ' + BUILD_TYPE)

# Build+install step.
run_command(r'mingw32-make install VERBOSE=1')

# Testing, packaging.
if is_python_build:
    # Run the Python tests.
    run_command(
        pinterp + r' -c "from pykep import test; test.run_test_suite();"')
    # Build the wheel.
    import shutil
    os.chdir('wheel')
    shutil.move(module_install_path, r'.')
    wheel_libs = 'mingw_wheel_libs_python{}.txt'.format(python_version)
    DLL_LIST = [_[:-1] for _ in open(wheel_libs, 'r').readlines()]
    for _ in DLL_LIST:
        shutil.copy(_, 'pykep/core')
    run_command(pinterp + r' setup.py bdist_wheel')
    os.environ['PATH'] = ORIGINAL_PATH
    run_command(pip + r' install dist\\' + os.listdir('dist')[0])

    os.chdir('/')
    run_command(
        pinterp + r' -c "from pykep import test; test.run_test_suite();"')
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
