#pragma once

#include "../writefilelib/writefile.hpp"
#include <cstdlib>
#include <fstream> 
#include <chrono> 
#include <iomanip>
#include <sstream>
#include <algorithm>

#define TIME chrono::high_resolution_clock::now(); 
#define DURATION chrono::duration<double> duration; 

using namespace std;

// Global constants for fluid dynamics

constexpr double Ru = 8314;
constexpr double Pr = 0.71;
constexpr double perfgam = 1.4;
constexpr double perfR = 287.0;

ThermoEntry operator*(const double& s, const ThermoEntry& A);

ThermoEntry operator+(const ThermoEntry& A, const ThermoEntry& B);


inline double computeInternalEnergy(const Vector& U) {
	return U[3] / U[0] - 1.0 / (2.0 * U[0] * U[0]) * (U[1] * U[1] + U[2] * U[2]);
}
inline double computePressure(const Vector& U, double& gam) {
	double e = computeInternalEnergy(U);
	return (gam - 1.0) * U[0] * e;
}
inline double computeTemperature(const Vector& U, ThermoEntry& Thermo) {
	double e = computeInternalEnergy(U);
	return e / Thermo.cv;
}
inline double computeSoundSpeed(ThermoEntry& Thermo) {

	return sqrt(Thermo.gamma * Thermo.R * Thermo.T);
}

// This struct contains the states for inviscid Jacobian computation
struct Inviscid_State {
	double rho, u, v, p, a, k, uprime, pp, pe, h0;
};

// This function computes the states for the inviscid Jacobians
inline Inviscid_State compute_inviscid_state(const Vector& U, ThermoEntry& thermo, double nx, double ny) {
	Inviscid_State S;

	S.rho = U[0];
	S.u = U[1] / U[0];
	S.v = U[2] / U[0];
	S.p = thermo.p;
	S.a = computeSoundSpeed(thermo);
	S.k = 1 / (S.a * S.a);
	S.uprime = S.u * nx + S.v * ny;
	S.pp = thermo.dpdrho - thermo.e / U[0] * thermo.dpde + 1 / (2 * U[0]) * (S.u * S.u + S.v * S.v) * thermo.dpde;
	S.pe = 1 / U[0] * thermo.dpde;
	S.h0 = thermo.e + 0.5 * (S.u * S.u + S.v * S.v);
	return S;
}

// This set boundary conditions based on text from the UI
vector<vector<ThermoEntry>> load_csv(const string& filename);
int find_index(double target, double min, double max, int n);
ThermoEntry bilinear_interpolate(
	const std::vector<std::vector<ThermoEntry>>& table,
	const std::vector<double>& rho_vec,
	double rho, double e);
vector<double> custom_rho_spacing();

Vector primtoCons(const Vector& V, double gam);
Vector constoPrim(const Vector& U, double gam);

unique_ptr<Grid> find_grid(string& gridname, int Nx, int Ny) ; 
tuple<int, int, Tensor, unique_ptr<Grid>, string, Vector, BCMap> restart_solution(string& filename);

class Solver {

private:

	string gridtype, restart_name;

	bool chemical_eqbm_enabled, restart; 
	const int Nx, Ny, progress_update;
	double CFL, Tw, dt, inner_residual, t_tot;

	Vector V_inlet, U_inlet, Global_Residual, t, iteration, rho_vec;
	vector<vector<ThermoEntry>> chem_lookup_table, cell_thermo;
	Tensor U, dU_new, dU_old, i_Fluxes, j_Fluxes, i_rho_fluxes, j_rho_fluxes;
	Tesseract i_plus_inviscid_Jacobians, i_minus_inviscid_Jacobians, i_viscous_Jacobians, j_plus_inviscid_Jacobians, j_minus_inviscid_Jacobians, j_viscous_Jacobians, i_rho_A, j_rho_A;

	Grid& grid;
	BCMap BCs;  
	inlet_conditions INLET;



public:

	double outer_residual;

	Solver(const int Nx, const int Ny, Vector U_inlet, Grid& grid, string& gridtype, BCMap BCs,  
		double CFL, double Tw, int& progress_update, bool& restart, string& restart_name, bool& chemical_eqbm_enabled, Tensor& U);


	Matrix inviscid_boundary_2D_E(BCType type, const Vector& U, const Point& normals);
	Vector inviscid_boundary_2D_U(BCType type, const Vector& U, const Point& normals);


	void get_chemistry();
	void get_perf_chemistry();
	void get_ghost_cells();
	void compute_dt();

	void solve_inviscid();
	void solve_inviscid_timestep();

	void restart_solution(string& filename); 

	void compute_inviscid_jacobians();
	void solve_left_line_inviscid();
	void solve_middle_line_inviscid(const int i);
	void solve_right_line_inviscid();

	void compute_inner_residual();
	void compute_outer_residual();

	void write_residual_csv();


	void time(void (Solver::* func)()) {
		auto start = TIME;

		(this->*func)();

		auto end = TIME;
		DURATION duration = end - start;
		cout << "Time taken: " << duration.count() << endl;
	}
	Vector minmod(Vector& a, Vector& b);

	//Vector constoViscPrim(const Vector& U, ThermoEntry& Thermo);
	//Matrix viscous_boundary_2D_E(BoundaryCondition type, const Vector& U, const Point& normals);
	//Vector viscous_boundary_2D_U(BoundaryCondition type, const Vector& U, const Point& normals);
	//void viscous_calculations();
	//void solve_viscous(); 
	//void solve_viscous_timestep(); 
	//void compute_viscous_jacobians(); 

	//void solve_left_line_viscous();
	//void solve_middle_line_viscous(const int i);
	//void solve_right_line_viscous();
};

//
//
//struct Viscous_State {
//	double rho, u, v, p, a, k, uprime, pp, h0, T;
//};
//inline Viscous_State compute_viscous_state(const Vector& U, double nx, double ny) {
//	Viscous_State S;
//
//	S.rho = U[0];
//	S.u = U[1] / S.rho;
//	S.v = U[2] / S.rho;
//	S.p = computePressure(U);
//	S.T = S.p / (S.rho * R);
//	S.a = sqrt(gam * S.p / S.rho);
//	S.k = 1 / (S.a * S.a);
//	S.uprime = S.u * nx + S.v * ny;
//	S.pp = 0.5 * (gam - 1) * (S.u * S.u + S.v * S.v);
//	S.h0 = (U[3] + S.p) / S.rho;
//	return S;
//}
//
//
//
// 