cmake_minimum_required(VERSION 3.10)
project(DPLRFiniteVolume)

add_executable(DPLRFiniteVolume "DPLR Finite Volume.cpp")

# Add the directories containing the source files (cpp) directly
target_sources(DPLRFiniteVolume PRIVATE
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\2D FVS Library\\2D FVS Library.cpp"
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\Linear Algebra Library\\Linear Algebra Library.cpp"
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\GridGenerator\\GridGenerator.cpp"
)

# Include the directories containing header files (h)
target_include_directories(DPLRFiniteVolume PUBLIC
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\2D FVS Library"
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\Linear Algebra Library"
    "C:\\Users\\Connor\\source\\repos\\Directed Study Spring 2025\\GridGenerator"
)
