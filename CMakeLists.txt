cmake_minimum_required(VERSION 3.12)

# Set project name, language and c++ std
project("SPARSE_UTILS"
  VERSION 1.0.0
  DESCRIPTION "A collection of sparse matrix data types and sparse solvers."
  HOMEPAGE_URL "https://github.com/smdogroup/sparse-utils"
  LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
if(NOT DEFINED SPARSE_UTILS_METIS_DIR OR NOT EXISTS ${SPARSE_UTILS_METIS_DIR})
  message(FATAL_ERROR "please specify the path to a metis installation by using
  -DSPARSE_UTILS_METIS_DIR=...")
endif()

# Some path variables
include(GNUInstallDirs)

option(SPARSE_UTILS_BUILD_TESTS "Build unit tests" ${PROJECT_IS_TOP_LEVEL})
option(SPARSE_UTILS_INSTALL_LIBRARY "Enable installation" ${PROJECT_IS_TOP_LEVEL})

find_package(LAPACK REQUIRED)

# Interface makes the target header-only that does not need to be compiled
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Dependencies
target_include_directories(${PROJECT_NAME} INTERFACE
  ${SPARSE_UTILS_METIS_DIR}/include)
target_link_directories(
  ${PROJECT_NAME} INTERFACE ${SPARSE_UTILS_METIS_DIR}/lib)

target_link_libraries(${PROJECT_NAME} INTERFACE metis LAPACK::LAPACK)

# Pass along c++17 requirements to projects using this library
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

# # Set warning flags
# # TODO: will this work for downstream projects? No
# if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|GNU")
# add_compile_options(
# -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare
# )
# endif()
if(SPARSE_UTILS_INSTALL_LIBRARY)
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

  include(CMakePackageConfigHelpers)
  write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

  configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION
    ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  # install cmake files
  install(
    EXPORT ${PROJECT_NAME}_Targets
    FILE "${PROJECT_NAME}Targets.cmake"
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
  )
  install(
    FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
  )

  # Install headers
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)
endif()

# # Build unit tests
# # Download google test (git tag needs to be updated frequently)
# if(SPARSE_UTILS_BUILD_TESTS)
# add_subdirectory(tests)
# enable_testing()
# endif()
