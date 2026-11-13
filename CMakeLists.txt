cmake_minimum_required(VERSION 3.14)
project(voltspice CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-Wall -Wextra)

add_library(voltspice_core src/circuit.cpp src/parser.cpp src/solver.cpp src/simulator.cpp)
target_include_directories(voltspice_core PUBLIC src)

add_executable(voltspice src/main.cpp)
target_link_libraries(voltspice voltspice_core)

enable_testing()
add_executable(test_voltspice tests/test_voltspice.cpp)
target_link_libraries(test_voltspice voltspice_core)
add_test(NAME voltspice COMMAND test_voltspice)

add_executable(test_sparse tests/test_sparse.cpp)
target_link_libraries(test_sparse voltspice_core)
add_test(NAME sparse_vs_dense COMMAND test_sparse)
