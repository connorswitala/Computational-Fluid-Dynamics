cmake_minimum_required(VERSION 3.10)
project(RealGasSolver)

# Original source file for the executable
add_executable(RealGasSolver "Real Gas Solver.cpp")

# Rename output executable
set_target_properties(RealGasSolver PROPERTIES
    OUTPUT_NAME "Real Gas Solver"
)

# Add the source files
target_sources(RealGasSolver PRIVATE
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\2D FVS Library\\2D FVS Library.cpp"
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\Linear Algebra Library\\Linear Algebra Library.cpp"
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\GridGenerator\\GridGenerator.cpp"
)

# Include the header directories
target_include_directories(RealGasSolver PUBLIC
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\2D FVS Library"
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\Linear Algebra Library"
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\GridGenerator"
)
