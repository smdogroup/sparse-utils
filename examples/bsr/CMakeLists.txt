add_executable(bsr_to_dense bsr_to_dense.cpp)
target_link_libraries(bsr_to_dense PRIVATE SparseUtils::SparseUtils)
