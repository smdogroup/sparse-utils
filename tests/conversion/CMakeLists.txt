# Add targets
add_executable(test_bsr_to_csr_csc test_bsr_to_csr_csc.cpp)

target_include_directories(
    test_bsr_to_csr_csc
    PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests/include)

# Link to the default main from Google Test
target_link_libraries(test_bsr_to_csr_csc gtest_main)

# Make tests auto-testable with CMake ctest
include(GoogleTest)
gtest_discover_tests(test_bsr_to_csr_csc)
