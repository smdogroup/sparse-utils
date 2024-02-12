# SparseUtils - Sparse matrix utilities for numerical discretization applications

sparse-utils is a header-only library that contains sparse matrix data
structures and sparse direct/iterator solvers.

## Usage: CMake


## Dependencies
- BLAS/Lapack
- metis

## CMake variables

Below is the complete table of CMake variables that sparse-utils accepts to
control the compilation.

_Recall that to give the variable VARIABLE, its value VAL can be set using the
following syntax in the command line:_
```
cmake ... -DVARIABLE=VAL ...
```

| Variable | Description | Default | Choices |
|----------|-------------|---------|---------|
|SPARSE_UTILS_METIS_DIR|path to metis installation| No default |a path|
