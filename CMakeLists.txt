cmake_minimum_required(VERSION 3.20)
project(ASC_bla)

# Set default CMAKE_BUILD_TYPE to RelWithDebInfo if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
endif()


# C++17 required for constexpr if
set (CMAKE_CXX_STANDARD 17)

# Set default CMAKE\_BUILD\_TYPE to RelWithDebInfo if not specified
if(NOT CMAKE\_BUILD\_TYPE)
  set(CMAKE\_BUILD\_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
endif()


include_directories(src)

find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)

# BLAS/LAPACK detection
find_package(OpenBLAS)
if(NOT OpenBLAS_FOUND)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()


add_subdirectory (demos)

if(WIN32)
   file (DOWNLOAD https://sourceforge.net/projects/openblas/files/v0.3.24/OpenBLAS-0.3.24-x64.zip/download "${CMAKE_BINARY_DIR}/openblas.zip" EXPECTED_MD5 fc08fe8c0dc7364da115d0e09b5a134f)
   file (ARCHIVE_EXTRACT INPUT ${CMAKE_BINARY_DIR}/openblas.zip DESTINATION ${CMAKE_BINARY_DIR}/openblas)
   set(BLAS_LIBRARIES ${CMAKE_BINARY_DIR}/openblas/lib/libopenblas.dll.a CACHE PATH "BLAS" FORCE)
   set(LAPACK_LIBRARIES ${CMAKE_BINARY_DIR}/openblas/lib/libopenblas.dll.a CACHE PATH "LAPACK" FORCE)
endif()



find_package(LAPACK REQUIRED)
add_executable (test_lapack demos/test_lapack.cpp)
target_link_libraries (test_lapack PUBLIC LAPACK::LAPACK)


pybind11_add_module(bla src/bind_bla.cpp)

#  fixing "undefined symbol: dgemm_"
if (TARGET OpenBLAS::OpenBLAS)
  target_link_libraries(bla PRIVATE OpenBLAS::OpenBLAS)
elseif (TARGET LAPACK::LAPACK)
  target_link_libraries(bla PRIVATE LAPACK::LAPACK)
else()
  target_link_libraries(bla PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

install (TARGETS bla DESTINATION ASCsoft)
install (FILES src/vector.hpp DESTINATION ASCsoft/include)


# Copy openblas dll file to build directory (no need to set PATH)
if(WIN32)
    add_custom_command(TARGET test_lapack POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${CMAKE_BINARY_DIR}/openblas/bin/libopenblas.dll"
      $<TARGET_FILE_DIR:test_lapack>
    )
endif()

