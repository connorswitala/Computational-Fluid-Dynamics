cmake_minimum_required(VERSION 3.10)
project(11-SpeciesGibbsMinimizer)

# Create executable
add_executable(11-SpeciesGibbsFreeEnergyMinimizer 
    "11-SpeciesGibbsMinimizer.cpp"
    "C:/Users/Connor/source/repos/Computational Fluid Dynamics/2D Perfect Gas FVS Library/2D Perfect Gas FVS Library.cpp"
    "C:/Users/Connor/source/repos/Computational Fluid Dynamics/Linear Algebra Library/Linear Algebra Library.cpp"
    "C:/Users/Connor/source/repos/Computational Fluid Dynamics/Grid Generator/GridGenerator.cpp"
)

# Add include directories
target_include_directories(11-SpeciesGibbsFreeEnergyMinimizer PUBLIC
    "C:/Users/Connor/source/repos/Computational Fluid Dynamics/2D Perfect Gas FVS Library"
    "C:/Users/Connor/source/repos/Computational Fluid Dynamics/Linear Algebra Library"
    "C:/Users/Connor/source/repos/Computational Fluid Dynamics/Grid Generator"
)


