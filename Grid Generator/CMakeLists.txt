cmake_minimum_required(VERSION 3.10)

# Define the project name
project(GridGenerator)

# Create a library for GridGenerator
add_library(GridGenerator STATIC "GridGenerator.cpp")

# Link to the LinearAlgebraLibrary target
target_link_libraries(GridGenerator PUBLIC LinearAlgebraLibrary)

# Include directories for the GridGenerator project and LinearAlgebraLibrary
target_include_directories(GridGenerator PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}  # Include current directory (for GridGenerator.h)
    ${CMAKE_SOURCE_DIR}/LinearAlgebraLibrary  # Include LinearAlgebraLibrary directory
)