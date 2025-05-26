cmake_minimum_required(VERSION 3.10)
project(LinearAlgebraLibrary)

add_library(LinearAlgebraLibrary STATIC "Linear Algebra Library.cpp") 

target_include_directories(LinearAlgebraLibrary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})