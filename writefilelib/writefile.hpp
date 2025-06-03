#pragma once

#include <fstream> 
#include <iomanip>
#include <iostream>
#include "../gridlib/grid.hpp"


Vector convert(const Vector& U);

void cfd_contour(const int Nx, const int Ny, Tensor& U, Grid& grid, string filename);

void cfd_centerline(const int Nx, const int Ny, Tensor& U, Grid& grid, string filename);