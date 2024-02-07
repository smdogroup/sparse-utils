include_directories(${SPARSE_UTILS_ROOT_DIR}/include)
include_directories(${SPARSE_UTILS_ROOT_DIR}/tests) # test_commons.h

# Add targets
add_executable(test_bsr_to_csr_csc test_bsr_to_csr_csc.cpp)

# Link to the default main from Google Test
target_link_libraries(test_bsr_to_csr_csc gtest_main)

# Make tests auto-testable with CMake ctest
include(GoogleTest)
gtest_discover_tests(test_bsr_to_csr_csc)
