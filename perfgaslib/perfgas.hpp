#pragma once

#include "../writefilelib/writefile.hpp" 
#include <chrono> 
#include <iomanip>
#include <tuple>

#define TIME chrono::high_resolution_clock::now(); 
#define DURATION chrono::duration<double> duration; 

using namespace std;

// Global constants for fluid dynamics
constexpr double gam = 1.4;
constexpr double R = 287;
constexpr double Ru = 8314;
constexpr double cv = R / (gam - 1);
constexpr double cp = cv + R;
constexpr double Pr = 0.71;

inline double computeInternalEnergy(const Vector& U) {
	double rho = U[0];
	double KE = 0.5 * U[0] * ((U[1] * U[1]) / (rho * rho) + (U[2] * U[2]/(rho * rho)));
	return (U[3] - KE)/rho;
}

// Inline function that computes pressure from state vector
inline double computePressure(const Vector& U) {
	double e = computeInternalEnergy(U); 
	return (gam - 1) * U[0] * e;  
}

// Inline function that compuites Temperature from state vector
inline double computeTemperature(const Vector& U) {
	double e = computeInternalEnergy(U);
	return e / cv;
}

// This struct contains the states for inviscid Jacobian computation
struct Inviscid_State {
	double rho, u, v, p, a, k, uprime, pp, h0;
};

// This struct contains the states for viscous Jacobians computation
struct Viscous_State {
	double rho, u, v, p, a, k, uprime, pp, h0, T;
};

// This function computes the states for the inviscid Jacobians
inline Inviscid_State compute_inviscid_state(const Vector& U, double nx, double ny) {
	Inviscid_State S;
	S.rho = U[0];
	S.u = U[1] / S.rho;
	S.v = U[2] / S.rho;
	S.p = computePressure(U);
	S.a = sqrt(gam * S.p / S.rho);
	S.k = 1 / (S.a * S.a);
	S.uprime = S.u * nx + S.v * ny;
	S.pp = 0.5 * (gam - 1) * (S.u * S.u + S.v * S.v);
	S.h0 = U[3] / U[0];
	return S;
}

// This function computes the states for the viscous Jacobians
inline Viscous_State compute_viscous_state(const Vector& U, double nx, double ny) {
	Viscous_State S;

	S.rho = U[0];
	S.u = U[1] / S.rho;
	S.v = U[2] / S.rho;
	S.p = computePressure(U);
	S.T = S.p / (S.rho * R);
	S.a = sqrt(gam * S.p / S.rho);
	S.k = 1 / (S.a * S.a);
	S.uprime = S.u * nx + S.v * ny;
	S.pp = 0.5 * (S.u * S.u + S.v * S.v) * (gam - 1);
	S.h0 = U[3] / U[0];
	return S;
}

Vector primtoCons(const Vector& V);
Vector constoPrim(const Vector& U);

tuple<int, int, Tensor, unique_ptr<Grid>, string, Vector, BCMap>  restart_solution(string& filename);  
unique_ptr<Grid> find_grid(string& gridname, int Nx, int Ny);  

class Solver {

private:

	string gridtype, restart_name;

	bool restart;
	int Nx, Ny, progress_update;
	double CFL, Tw, dt, inner_residual, t_tot, outer_residual;   

	Vector U_inlet, Global_Residual, t, iteration;
	Tensor U, dU_new, dU_old, i_Fluxes, j_Fluxes;
	Tesseract i_plus_inviscid_Jacobians, i_minus_inviscid_Jacobians, i_viscous_Jacobians, j_plus_inviscid_Jacobians, j_minus_inviscid_Jacobians, j_viscous_Jacobians;

	Grid& grid; 
	BCMap BCs; 
	inlet_conditions INLET;


public:

	Solver(int Nx, int Ny, Vector U_inlet, Grid& grid, string& gridtype, BCMap BCs, double CFL, double Tw,  
		int& progress_update, bool& restart, string& restart_name, Tensor& U);

	Vector constoViscPrim(const Vector& U);

	Matrix inviscid_boundary_2D_E(BCType type, const Vector& U, const Point& normals);
	Vector inviscid_boundary_2D_U(BCType type, const Vector& U, const Point& normals);

	Matrix viscous_boundary_2D_E(BCType type, const Vector& U, const Point& normals);
	Vector viscous_boundary_2D_U(BCType type, const Vector& U, const Point& normals);

	void solve_inviscid();
	void solve_viscous();
	void solve_inviscid_timestep();
	void solve_viscous_timestep();
	void compute_dt();
	void compute_ghost_cells(); 
	void compute_inviscid_jacobians();
	void compute_viscous_jacobians();

	void solve_left_line_inviscid();
	void solve_middle_line_inviscid(const int i);
	void solve_right_line_inviscid();

	void solve_left_line_viscous();
	void solve_middle_line_viscous(const int i);
	void solve_right_line_viscous();


	void compute_inner_residual();
	void compute_outer_residual();

	void viscous_calculations();

	void write_2d_csv(const string& filename);
	void write_1d_csv(const string& filename);
	void write_residual_csv();

	Vector minmod(Vector& Ui, Vector& Uii);
};