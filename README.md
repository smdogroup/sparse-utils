# SparseUtils - Sparse matrix utilities for numerical discretization applications

[SparseUtils](https://github.com/smdogroup/sparse-utils) is a header-only
library that contains sparse matrix data structures and sparse direct/iterator
solvers.

## Dependencies
- LAPACK
- metis

## Installation

### CMake

CMake is preferred to install SparseUtils. For basic installation, use the
following command:
```
mkdir build && cd build && cmake .. && make install
```
This installs SparseUtils (headers and CMake files) into
```${HOME}/installs/sparse-utils```.


### In-tree
Not yet supported.

### Raw Makefile
Not yet supported.

## Use

Using SparseUtils is simple with CMake. For the application project, add
```
find_package(SparseUtils REQUIRED PATHS <path-to-sparse-utils-installation>)
```
to ```CMakeLists.txt```, and use
```
target_link_libraries(<app-target> SparseUtils::SparseUtils)
```
in the ```CMakeLists.txt``` for the application executables. See
[examples/CMakeLists.txt](examples/CMakeLists.txt) for example.



## Test
To build and execute unit tests, using the following command:
```
mkdir build && cd build && cmake .. -DSPARSE_UTILS_BUILD_TESTS=ON && make -j && ctest .
```

## CMake variables

Below is the complete table of CMake variables that SparseUtils uses to
control the compilation.

_Recall that to give the variable VARIABLE, its value VAL can be set using the
following syntax in the command line:_
```
cmake ... -DVARIABLE=VAL ...
```

| Variable | Description | Default | Choices |
|----------|-------------|---------|---------|
|SPARSE_UTILS_METIS_DIR|path to metis installation|```${HOME}/installs/metis```|a path|
|CMAKE_INSTALL_PREFIX|the path to install SparseUtils to|```${HOME}/installs/sparse-utils```|a path|
|SPARSE_UTILS_BUILD_TESTS|build unit tests or not|OFF|ON, OFF|