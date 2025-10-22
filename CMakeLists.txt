cmake_minimum_required(VERSION 3.20)
project(ASC_bla)

# C++17 required for constexpr if
set (CMAKE_CXX_STANDARD 17)

include_directories(src)


find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)


add_subdirectory (demos)


pybind11_add_module(bla src/bind_bla.cpp)

install (TARGETS bla DESTINATION ASCsoft)
install (FILES src/vector.hpp DESTINATION ASCsoft/include)


