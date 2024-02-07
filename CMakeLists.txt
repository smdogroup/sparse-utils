# Specify cmake version needed
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# Set project name, language and c++ std
project(SPARSE_UTILS LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SPARSE_UTILS_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Set warning flags
if (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|GNU")
add_compile_options(
  -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare
)
endif()

# Build unit tests
# Download google test (git tag needs to be updated frequently)
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG e2239ee6043f73722e7aa812a459f54a28552929 # release-1.10.0
)
FetchContent_MakeAvailable(googletest)
enable_testing()
add_subdirectory(${SPARSE_UTILS_ROOT_DIR}/tests)

# Build examples
add_subdirectory(${SPARSE_UTILS_ROOT_DIR}/examples)
