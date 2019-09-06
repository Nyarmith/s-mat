cmake_minimum_required(VERSION 3.10)
project(Static-Matrices)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(s-mat INTERFACE)

target_include_directories(s-mat INTERFACE include/)


if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(test/)
endif()
