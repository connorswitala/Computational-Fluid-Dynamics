#pragma once

#include <fstream> 
#include <iomanip>
#include <iostream>
#include "../gridlib/grid.hpp"
#include "../linalglib/linalg.hpp"

struct ThermoEntry {
	double rho, e, p, T, R, cv, gamma, dpdrho, dpde;
};


Vector convert(const Vector& U);
Vector convert_real(const Vector& U, ThermoEntry thermo);

void cfd_contour(const int Nx, const int Ny, Tensor& U, Grid& grid, string filename);
void cfd_real_contour(const int Nx, const int Ny, Tensor& U, Grid& grid, vector<vector<ThermoEntry>>& cell_thermo, string filename);   


void cfd_centerline(const int Nx, const int Ny, Tensor& U, Grid& grid, string filename);

void output_grid(Grid& grid, const int& Nx, const int& Ny, string& filename);
