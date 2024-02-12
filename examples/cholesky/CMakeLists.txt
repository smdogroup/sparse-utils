add_executable(cholesky cholesky.cpp)
target_link_libraries(cholesky PRIVATE SparseUtils::SparseUtils)
