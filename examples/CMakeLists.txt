cmake_minimum_required(VERSION 3.12)

project(SPARSE_UTILS_EXAMPLES)

if(NOT DEFINED SPARSE_UTILS_INSTALL_DIR)
    set(SPARSE_UTILS_INSTALL_DIR $ENV{HOME}/installs/sparse-utils)
endif()

if(NOT EXISTS ${SPARSE_UTILS_INSTALL_DIR})
    message(FATAL_ERROR "${SPARSE_UTILS_INSTALL_DIR} does not exist, specify the\
    path to the SparseUtils installation by -DSPARSE_UTILS_INSTALL_DIR=..."
    )
endif()

find_package(SparseUtils REQUIRED PATHS ${SPARSE_UTILS_INSTALL_DIR})

add_subdirectory(bsr)
add_subdirectory(cholesky)
