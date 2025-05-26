cmake_minimum_required(VERSION 3.10)

# Define the project name
project(2DRealGasFVSLibrary) 

# Create a static library for 2DFVSLibrary
add_library(2DRealGasFVSLibrary STATIC "2D Real Gas FVS Library.cpp") 

# Link to GridGenerator and LinearAlgebraLibrary libraries
target_link_libraries(2DRealGasFVSLibrary PUBLIC GridGenerator LinearAlgebraLibrary)

# Include directories for 2DFVSLibrary, GridGenerator, and LinearAlgebraLibrary
target_include_directories(2DRealGasFVSLibrary PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}  # Include current directory (for 2DFVSLibrary.h)
    ${CMAKE_SOURCE_DIR}/GridGenerator  # Include GridGenerator directory
    ${CMAKE_SOURCE_DIR}/LinearAlgebraLibrary  # Include LinearAlgebraLibrary directory
)
