cmake_minimum_required(VERSION 3.12)

project(SPARSE_UTILS_EXAMPLES)

if(NOT SPARSE_UTILS_INSTALL_DIR OR NOT EXISTS ${SPARSE_UTILS_INSTALL_DIR})
    message(FATAL_ERROR "specify the path to the SparseUtils installation by "
        "-DSPARSE_UTILS_INSTALL_DIR=..."
    )
endif()

find_package(SparseUtils REQUIRED PATHS ${SPARSE_UTILS_INSTALL_DIR})

add_subdirectory(bsr)
add_subdirectory(cholesky)
