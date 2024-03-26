cmake_minimum_required(VERSION 3.12)

# Set project name, language and c++ std
project("SparseUtils"
  VERSION 1.0.0
  DESCRIPTION "A collection of sparse matrix data types and sparse solvers."
  HOMEPAGE_URL "https://github.com/smdogroup/sparse-utils"
  LANGUAGES CXX)

include(cmake/utils.cmake)

# Dependencies
if(NOT DEFINED SPARSE_UTILS_METIS_DIR)
  set(SPARSE_UTILS_METIS_DIR $ENV{HOME}/installs/metis)
endif()

if(NOT DEFINED SPARSE_UTILS_LAPACKE_DIR)
  set(SPARSE_UTILS_LAPACKE_DIR /opt/homebrew/opt/lapack)
endif()

check_path_and_inform(metis ${SPARSE_UTILS_METIS_DIR} SPARSE_UTILS_METIS_DIR)
check_path_and_inform(lapacke ${SPARSE_UTILS_LAPACKE_DIR} SPARSE_UTILS_LAPACKE_DIR)

find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED PATHS ${SPARSE_UTILS_LAPACKE_DIR})

# Change the default value of CMAKE_INSTALL_PREFIX
set(SPARSE_UTILS_INSTALL_DIR $ENV{HOME}/installs/sparse-utils)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${SPARSE_UTILS_INSTALL_DIR} CACHE PATH "The default\
  path to install SparseUitls to" FORCE)
  message("The directory to install SparseUtils is set to default:\
   ${SPARSE_UTILS_INSTALL_DIR}")
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# Utilities
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(SPARSE_UTILS_BUILD_TESTS "Build unit tests" OFF)
option(SPARSE_UTILS_INSTALL_LIBRARY "Enable installation" ${PROJECT_IS_TOP_LEVEL})

# Interface makes the target header-only that does not need to be compiled
add_library(${PROJECT_NAME} INTERFACE)

# Use std::complex<double> instead of the default C99 complex type for LAPACKE
target_compile_definitions(${PROJECT_NAME} INTERFACE -Dlapack_complex_double=std::complex<double>)

# Pass along c++17 requirements to projects using this library
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

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
target_link_libraries(${PROJECT_NAME} INTERFACE metis LAPACK::LAPACK lapacke)

# Set warning flags
# TODO: specify warning flags for other compilers
if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|GNU")
  target_compile_options(${PROJECT_NAME} INTERFACE
    -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare
  )
endif()

if(SPARSE_UTILS_INSTALL_LIBRARY)
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

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

if(SPARSE_UTILS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()