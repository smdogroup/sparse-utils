# include sparse headers
include_directories(${SPARSE_UTILS_ROOT_DIR}/include)

# include metis
include_directories(${SPARSE_UTILS_METIS_DIR}/include)
link_directories(${SPARSE_UTILS_METIS_DIR}/lib)

# Add targets
add_executable(cholesky cholesky.cpp)

# Link libraries
target_link_libraries(cholesky LAPACK::LAPACK metis)

