cmake_minimum_required(VERSION 3.10)
project(Static-Matrices)

add_library(s-mat INTERFACE)

target_include_directories(s-mat INTERFACE include/)

enable_testing()
add_subdirectory(test/)
