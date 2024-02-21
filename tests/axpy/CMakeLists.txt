# Add targets
add_executable(test_axpy test_axpy.cpp)

target_include_directories(
    test_axpy
    PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests/include)

# Link to the default main from Google Test
target_link_libraries(test_axpy gtest_main)

# Make tests auto-testable with CMake ctest
include(GoogleTest)
gtest_discover_tests(test_axpy)
