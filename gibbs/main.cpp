
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include "../linalglib/linalg.hpp"
#include "../realgaslib/realgas.hpp"  

using namespace std;
constexpr int n_species = 10;

struct molar_concentration {
	double N2, O2, NO, N, O, Ar, Arp, Np, Op, em, total;
};

class Chemistry {
private:

	double rho, e, gam, T, p, R_mix, cv_mix;
	Vector CP_0, H_0, S_0, mu_k0, U, mk, Xk, Yk, Ck, R_k, MW, theta_v, V, h_f, theta_f, q;
	Matrix T_coeff, Int_const, hi_t_coeff, hi_t_int_const, lo_t_coeff, lo_t_int_const, middle_t_coeff, middle_t_int_const, a;

	molar_concentration initial_moles;

public:

	Chemistry() : rho(rho), e(e), T(T), p(p), R_mix(R_mix), cv_mix(cv_mix), gam(gam), q(n_species), a(3, Vector(n_species)),
		CP_0(n_species), H_0(n_species), S_0(n_species), mu_k0(n_species), U(n_species), mk(n_species), Xk(n_species),
		Yk(n_species), R_k(n_species), MW(n_species), theta_v(3), V(4), h_f(n_species), theta_f(n_species), Ck(n_species),
		T_coeff(n_species, Vector(7)), Int_const(n_species, Vector(2)), hi_t_coeff(n_species, Vector(7)), hi_t_int_const(n_species, Vector(2)), lo_t_coeff(n_species, Vector(7)),
		lo_t_int_const(n_species, Vector(2)), middle_t_coeff(n_species, Vector(7)), middle_t_int_const(n_species, Vector(2)) {

		// As a reminder, the order of species is  		
		// =[ N2, O2, NO, N, O, Ar, Ar+, N+, O+, e- ]=		
		// (Yk) is the mass fraction and (Xk) is the molar fraction. 

		MW = { 28.0134, 31.998, 30.008, 14.0067, 15.9994, 39.948, 39.9474514, 14.0061514, 15.9988514, 0.000548579903 };	// Molecular weights

		// Calculate species gas constants
		for (int i = 0; i < n_species; ++i) {
			R_k[i] = Ru / MW[i]; // Individual gas constants 
		}

		theta_v = { 3395.0, 2239.0, 2817.0 };	// Characteristic temperatures of vibration for N2, O2, NO
		theta_f = { 0.0, 0.0, 2.996120e+6, 3.362160e+7, 1.542000e+7, 0.0, 3.82155e7, 1.34337e8, 9.80594e7, 0.0000 };	// Enthalpies of formation 

		// N, O, Ar Stoichiometric coefficients
		a = {
			{
				{2, 0, 1, 1, 0, 0, 0, 1, 0, 0},
				{0, 2, 1, 0, 1, 0, 0, 0, 1, 0},
				{0, 0, 0, 0, 0, 1, 1, 0, 0, 0}
			}
		};

		// Charge of each species
		q = { 0, 0, 0, 0, 0, 0, 1, 1, 1, -1 };

		// Three sets of NASA Polynomials for certain temperature ranges [lo = 200 - 1000], middle = [1000 - 6000], hi = [6000 - 20000])
		// N2, O2, NO, N, O, Ar, Ar+, N+, O+, NO+ e-

		lo_t_coeff = {
			{
				{2.210371497e+4, -3.818461820e+2, 6.082738360, -8.530914410e-3, 1.384646189E-05, -9.625793620e-9, 2.519705809e-12},
				{-3.425563420e+4, 4.847000970e+2, 1.119010961, 4.293889240e-3, -6.836300520e-7, -2.023372700e-9, 1.039040018e-12},
				{-1.143916503e+4, 1.536467592e+2, 3.431468730, -2.668592368e-3, 8.481399120e-6, -7.685111050e-9, 2.386797655e-12},
				{ 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0},
				{-7.953611300e+3, 1.607177787e+2, 1.966226438, 1.013670310e-3, -1.110415423e-6, 6.517507500e-10, -1.584779251e-13},
				{ 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0 },
				{-57312.09170, 793.079147, -1.717121217, 0.01044184018, -1.180207501e-05, 6.52813478e-09, -1.44755813e-12},
				{5237.07921, 2.299958315, 2.487488821, 2.737490756e-05, -3.134447576e-08, 1.850111332e-11, -4.447350984e-15},
				{0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0},
				{0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0}
			}
		};
		lo_t_int_const = {
			{
				{7.108460860e+2, -1.076003744e+1},
				{-3.391454870e+3, 1.849699470e+1},
				{9.098214410e+3, 6.728725490},
				{5.610463780e+4, 4.193905036},
				{2.840362437e+4, 8.404241820},
				{-745.375, 4.37967491},
				{179057.223, 29.4915095},
				{225628.4738, 5.076830786},
				{187935.2842, 4.39337676},
				{-745.375, -11.72081224}
			}
		};
		middle_t_coeff = {
			{
				   { 5.877124060e+5, -2.239249073e+3, 6.066949220, -6.139685500e-4, 1.491806679e-7, -1.923105485e-11, 1.061954386e-15 },
				   { -1.037939022e+6, 2.344830282e+3, 1.819732036, 1.267847582e-3, -2.188067988e-7, 2.053719572e-11, -8.193467050e-16 },
				   { 2.239018716e+5, -1.289651623e+3, 5.433936030, -3.656034900e-4, 9.880966450e-8, -1.416076856e-11, 9.380184620e-16 },
				   { 8.876501380e+4, -1.071231500e+2, 2.362188287, 2.916720081e-4, -1.729515100e-7, 4.012657880e-11, -2.677227571e-15 },
				   { 2.619020262e+5, -7.298722030e+2, 3.317177270, -4.281334360e-4, 1.036104594e-7, -9.438304330e-12, 2.725038297e-16 },
				   {20.10538475, -0.0599266107, 2.500069401, -3.99214116e-08, 1.20527214e-11, -1.819015576e-15, 1.078576636e-19},
				   {-383596.54, 816.20197, 2.301342628, -4.95298377e-06, 1.205108477e-08, -2.185050286e-12, 1.265493898e-16},
				   {290497.0374, -855.790861, 3.47738929, -5.28826719e-04, 1.352350307e-07, -1.389834122e-11, 5.046166279e-16},
				   {-216651.3208, 666.545615, 1.702064364, 4.71499281e-04, -1.427131823e-07, 2.016595903e-11, -9.107157762e-16},
				   {0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0}
			}
		};
		middle_t_int_const = {
			{
				{1.283210415e+4, -1.586640027e+1},
				{-1.689010929e+4, 1.738716506e+1},
				{1.750317656e+4, -8.501669090},
				{ 5.697351330e+4, 4.865231506 },
				{3.392428060e+4, -6.679585350e-1},
				{-744.993961, 4.37918011},
				{177181.1455, 7.94750748},
				{231080.9984, -1.994146545},
				{183719.1966, 10.05690382},
				{-745.375, -11.72081224}
			}
		};
		hi_t_coeff = {
			{
				{ 8.310139160e+8, -6.420733540e+5, 2.020264635e+2, -3.065092046e-2, 2.486903333e-6, -9.705954110e-11, 1.437538881e-15 },
				{ 4.975294300e+8, -2.866106874e+5, 6.690352250e+1, -6.169959020e-3, 3.016396027e-7, -7.421416600e-12, 7.278175770e-17 },
				{ -9.575303540e+8, 5.912434480e+5, -1.384566826e+2, 1.694339403e-2, -1.007351096e-6, 2.912584076e-11, -3.295109350e-16 },
				{ 5.475181050e+8, -3.107574980e+5, 6.916782740e+1, -6.847988130e-3, 3.827572400e-7, -1.098367709e-11, 1.277986024e-16 },
				{ 1.779004264e+8, -1.082328257e+5, 2.810778365e+1, -2.975232262e-3, 1.854997534e-7, -5.796231540e-12, 7.191720164e-17 },
				{-9.95126508e+08, 6.45888726e+05, -167.5894697, 0.02319933363, -1.721080911e-06, 6.53193846e-11, -9.740147729e-16},
				{10068848.27, -6624.36128, 4.4469082, -3.017567664e-04, 2.612882069e-08, -1.201637769e-12, 2.299206903e-17},
				{1.646092148e+07, -11131.65218, 4.97698664, -2.005393583e-04, 1.022481356e-08, -2.691430863e-13, 3.539931593e-18},
				{-2.143835383e+08, 1.469518523e+05, -36.8086454, 5.03616454e-03, -3.087873854e-07, 9.18683487e-12, -1.074163268e-16},
				{0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0}
			}
		};
		hi_t_int_const = {
			{
				{4.938707040e+6, -1.672099740e+3},
				{2.293554027e+6, -5.530621610e+2},
				{-4.677501240e+6, 1.242081216e+3},
				{ 2.550585618e+6, -5.848769753e+2 },
				{8.890942630e+5, -2.181728151e+2},
				{-5.07830034e+06, 1465.298484},
				{234950.4137, -10.32262257},
				{313628.4696, -17.06646380},
				{-961420.896, 342.619308},
				{-745.375, -11.72081224}
			}
		};

		// Set polynomials based on temperature
		if (T > 200.0 && T < 1000.0) {
			T_coeff = lo_t_coeff;
			Int_const = lo_t_int_const;
		}
		else if (T >= 1000.0 && T < 6000.0) {
			T_coeff = middle_t_coeff;
			Int_const = middle_t_int_const;
		}
		else {
			T_coeff = hi_t_coeff;
			Int_const = hi_t_int_const;
		}

	}

	void display_molar_fractions() {
		cout << "Molar Fractions: " << endl << endl;
		cout << "N2: " << Xk[0] << endl;
		cout << "O2: " << Xk[1] << endl;
		cout << "NO: " << Xk[2] << endl;
		cout << "N:  " << Xk[3] << endl;
		cout << "O:  " << Xk[4] << endl;
		cout << "Ar: " << Xk[5] << endl;
		cout << "Ar+: " << Xk[6] << endl;
		cout << "N+: " << Xk[7] << endl;
		cout << "O+:  " << Xk[8] << endl;
		cout << "e-:  " << Xk[9] << endl;
	}
	void display_mass_fraction() {
		cout << "Mass Fractions: " << endl << endl;
		cout << "N2:  " << Yk[0] << endl;
		cout << "O2:  " << Yk[1] << endl;
		cout << "NO:  " << Yk[2] << endl;
		cout << "N:   " << Yk[3] << endl;
		cout << "O:   " << Yk[4] << endl;
		cout << "Ar:  " << Yk[5] << endl;
		cout << "Ar+: " << Yk[6] << endl;
		cout << "N+:  " << Yk[7] << endl;
		cout << "O+:  " << Yk[8] << endl;
		cout << "e-:  " << Yk[9] << endl;
		cout << "Temperature: " << T << endl;
	}
	void compute_mass_fractions() {

		double MW_mix = 0.0; // Mixture Molecular Weight 

		for (int i = 0; i < n_species; ++i) {
			MW_mix += Xk[i] * MW[i];
		}

		for (int i = 0; i < n_species; ++i) {
			Yk[i] = (Xk[i] * MW[i]) / MW_mix;
		}
	}
	void compute_h0() {
		for (int j = 0; j < n_species; ++j) {
			H_0[j] = Ru * T * (-T_coeff[j][0] * 1 / (T * T) + T_coeff[j][1] * log(T) / T + T_coeff[j][2] + T_coeff[j][3] * T / 2
				+ T_coeff[j][4] * T * T / 3 + T_coeff[j][5] * T * T * T / 4 + T_coeff[j][6] * T * T * T * T / 5 + Int_const[j][0] / T);
		}
	}
	void compute_s0() {
		for (int j = 0; j < n_species; ++j) {
			S_0[j] = Ru * (-T_coeff[j][0] * 1 / (2 * T * T) - T_coeff[j][1] / T + T_coeff[j][2] * log(T) + T_coeff[j][3] * T
				+ T_coeff[j][4] * T * T / 2 + T_coeff[j][5] * T * T * T / 3 + T_coeff[j][6] * T * T * T * T / 4 + Int_const[j][1]);
		}
	}
	void compute_mu0() {

		compute_h0();
		compute_s0();

		for (int j = 0; j < n_species; ++j) {
			mu_k0[j] = H_0[j] - T * S_0[j];
		}
	}
	double norm(Vector& v1, Vector& v2) {
		double result = 0.0;
		for (int i = 0; i < 14; ++i) {
			result += fabs(v1[i] - v2[i]) * fabs(v1[i] - v2[i]);
		}
		return sqrt(result);
	}
	bool safe_compute_molar_fractions() {  // Function that calls compute_molar_fractions(). At lower temperatures it fails so it just reverts to perfect gas assumption. 
		try {
			compute_molar_fractions(); // Try to solve chemical equilibrium
			return true;
		}
		catch (const exception& ex) {
			cerr << "[Warning] Chemical equilibrium failed: " << ex.what() << std::endl;

			// Fallback: assume frozen air (N2 and O2 only)
			// N2 = 78%, O2 = 21%, Ar = 0.97%, everything else 0%
			Xk = { 0.7808, 0.2095, 0.0, 0.0, 0.0, 0.0097, 0.0, 0.0, 0.0, 0.0 };
			return false;
		}
	}

	// Gibbs free-energy minimizaton proceudure 
	void compute_molar_fractions() {

		compute_mu0(); // Recalculate reference chemical potentials
		Vector X_new(14), X(14), dx(14);
		double ph = initial_moles.total / n_species;

		// Initialize guess: equal distribution among species
		X = { ph, ph, ph, ph, ph, ph, ph, ph, ph, ph, 1e-4, 1e-4, 1e-4, 1e-4 };
		X_new = X / 2; // Slightly disturbed initial guess


		double residual = norm(X_new, X);

		int iteration = 0;
		const int max_iterations = 1000; // Cap iterations to prevent infinite loops

		while (residual >= 1e-10 && iteration < max_iterations) {

			iteration++;

			// Reset intermediate variables

			X[10] = 0; X[11] = 0; X[12] = 0; X[13] = 0;

			// Build Jacobian
			Matrix J = {
				{
				{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -a[0][0], -a[1][0], -a[2][0], -q[0]}, // Row 0: N2  
				{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -a[0][1], -a[1][1], -a[2][1], -q[1]}, // Row 1: O2  
				{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -a[0][2], -a[1][2], -a[2][2], -q[2]}, // Row 2: NO 
				{0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -a[0][3], -a[1][3], -a[2][3], -q[3]}, // Row 3: N 
				{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -a[0][4], -a[1][4], -a[2][4], -q[4]}, // Row 4: O
				{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -a[0][5], -a[1][5], -a[2][5], -q[5]}, // Row 5: Ar 
				{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -a[0][6], -a[1][6], -a[2][6], -q[6]}, // Row 6: Ar+ 
				{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -a[0][7], -a[1][7], -a[2][7], -q[7]}, // Row 7: N+ 
				{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -a[0][8], -a[1][8], -a[2][8], -q[8]}, // Row 8: O+
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -a[0][9], -a[1][9], -a[2][9], -q[9]}, // Row 9: e-
				{a[0][0] * X[0], a[0][1] * X[1], a[0][2] * X[2], a[0][3] * X[3], a[0][4] * X[4], a[0][5] * X[5], a[0][6] * X[6], a[0][7] * X[7], a[0][8] * X[8], a[0][9] * X[9], 0, 0, 0, 0},  // Row 10: Nitrogen conservation
				{a[1][0] * X[0], a[1][1] * X[1], a[1][2] * X[2], a[1][3] * X[3], a[1][4] * X[4], a[1][5] * X[5], a[1][6] * X[6], a[1][7] * X[7], a[1][8] * X[8], a[1][9] * X[9], 0, 0, 0, 0},  // Row 11: Oxygen conservation
				{a[2][0] * X[0], a[2][1] * X[1], a[2][2] * X[2], a[2][3] * X[3], a[2][4] * X[4], a[2][5] * X[5], a[2][6] * X[6], a[2][7] * X[7], a[2][8] * X[8], a[2][9] * X[9], 0, 0, 0, 0},  // Row 12: Argon conservation
				{q[0] * X[0], q[1] * X[1], q[2] * X[2], q[3] * X[3], q[4] * X[4], q[5] * X[5], q[6] * X[6], q[7] * X[7], q[8] * X[8], q[9] * X[9], 0, 0, 0, 0},		// Row 13: Charge neutrality 
				}
			};

			// Gibbs free energy vector

			Vector F(14);
			for (int i = 0; i < n_species; ++i) {
				double Xi_safe = max(X[i], 1e-10); // Protect log from zero  
				F[i] = -(mu_k0[i] + Ru * T * log(Xi_safe * p / 101325)) / (Ru * T);
			}


			F[10] = initial_moles.N2 - (2 * X[0] + X[2] + X[3] + X[7]); // Nitrogen  
			F[11] = initial_moles.O2 - (2 * X[1] + X[2] + X[4] + X[8]); // Oxygen
			F[12] = initial_moles.Ar - (X[5] + X[6]); // Argon
			F[13] = -(X[6] + X[7] + X[8] - X[9]); // Electron 

			// Solve for Newton update
			dx = F / J;

			// Update molar fractions safely
			for (int i = 0; i < 14; ++i) {

				double dx_safe = min(max(dx[i], -50.0), 50.0); // Clamp dx 
				X_new[i] = X[i] * exp(dx_safe);

				if (std::isnan(X_new[i]) || std::isinf(X_new[i])) {
					throw std::runtime_error("Nonphysical molar fraction computed (NaN or Inf)");
				}

				if (X_new[i] < 0) {
					throw std::runtime_error("Negative molar fraction computed");
				}

			}

			residual = norm(X_new, X);
			X = X_new;
		}

		if (iteration >= max_iterations) {
			throw std::runtime_error("Chemical equilibrium solver did not converge");
		}

		double sum = 0.0; // Normalize molar fractions 
		// Store final molar fractions
		for (int i = 0; i < n_species; ++i) {
			sum += X_new[i];
		}

		for (int i = 0; i < n_species; ++i) {
			Xk[i] = X_new[i] / sum; // Normalize molar fractions
		}
	}

	// Solve for chemical equilibrium temperature using energy conservation
	void compute_temperature() {

		double e_new = 0, cv_new = 0;

		while (fabs(e_new - e) >= 1) {

			if (T > 200.0 && T < 1000.0) {
				T_coeff = lo_t_coeff;
				Int_const = lo_t_int_const;
			}
			else if (T >= 1000.0 && T < 6000.0) {
				T_coeff = middle_t_coeff;
				Int_const = middle_t_int_const;
			}
			else {
				T_coeff = hi_t_coeff;
				Int_const = hi_t_int_const;
			}

			bool equilibrium_success = safe_compute_molar_fractions();
			compute_mass_fractions();

			if (equilibrium_success == false) {
				R_mix = 287.0;
				cv_mix = 717;
				T = e / cv_mix;
				p = rho * R_mix * T;
				break;
			}
			else {

				e_new = 0.0;

				for (int i = 0; i < 3; ++i) {
					e_new += Yk[i] * (2.5 * R_k[i] * T + R_k[i] * theta_v[i] / (exp(theta_v[i] / T) - 1) + theta_f[i]); //diatomic species
				}

				for (int i = 3; i < n_species; ++i) {
					e_new += Yk[i] * (1.5 * R_k[i] * T + theta_f[i]); // atomic species
				}

				cv_new = e_new / T;

				T = T - 0.1 * (e_new - e) / cv_new;


				double sum = 0.0;
				for (int i = 0; i < n_species; ++i) {
					sum += Xk[i] * MW[i];
				}

				R_mix = Ru / sum; // Mixture gas constant  
				p = rho * R_mix * T;
			}


		}

		cv_mix = cv_new;

	}

	Vector compute_equilibrium(double Rho, double E, Vector& initial_mol) {

		double sum = 0.0;
		R_mix = 0.0;
		initial_moles.total = 0.0;

		initial_moles.N2 = initial_mol[0];
		initial_moles.O2 = initial_mol[1];
		initial_moles.NO = initial_mol[2];
		initial_moles.N = initial_mol[3];
		initial_moles.O = initial_mol[4];
		initial_moles.Ar = initial_mol[5];
		initial_moles.Arp = initial_mol[6];
		initial_moles.Np = initial_mol[7];
		initial_moles.Op = initial_mol[8];
		initial_moles.em = initial_mol[9];


		for (int i = 0; i < n_species; ++i) {
			initial_moles.total += initial_mol[i];
			sum += initial_mol[i] * MW[i];
		}

		R_mix = Ru / sum; // Mixture gas constant    


		rho = Rho;
		e = E;
		T = e / 717;
		p = rho * R_mix * T;  // Initial pressure estimate	  

		compute_temperature();
		// Vector mass_fractions(11); 

		// for (int i = 0; i < n_species; ++i) {
		// 	mass_fractions[i] = Xk[i]; // Store mass fractions    
		// } 
		// mass_fractions[10] = T; // Store temperature  

		gam = 1 + R_mix / cv_mix;
		Vector thermo = { rho, e, p, T, R_mix, cv_mix, gam, 0.0, 0.0 }; 

		//////  Compute derivatives ////

		// dp/dp | e
		e = E + 1;
		T = e / 717;
		p = rho * 287.0 * T;
		compute_temperature();
		double p1 = p;


		e = E - 1;
		T = e / 717;
		p = rho * 287.0 * T;
		compute_temperature();
		double p2 = p;

		double dpde = (p1 - p2) / 2; // dp/de 


		// dp/de| rho
		rho = Rho + 0.00001;
		T = E / 717;
		p = rho * 287.0 * T;
		compute_temperature();
		p1 = p;

		rho = Rho - 0.00001;
		T = E / 717;
		p = rho * 287.0 * T;
		compute_temperature();
		p2 = p;
		double dpdrho = (p1 - p2) / 0.00002; // dp/drho

		thermo[7] = dpdrho; // dp/drho | e 
		thermo[8] = dpde; // dp/de | rho  

		return thermo;
	}
};

void display_entries(ThermoEntry& thermo) {
	cout << "rho: " << thermo.rho << endl;
	cout << "e: " << thermo.e << endl;
	cout << "p: " << thermo.p << endl;
	cout << "T: " << thermo.T << endl;
	cout << "R: " << thermo.R << endl;
	cout << "cv: " << thermo.cv << endl;
	cout << "gam: " << thermo.gamma << endl;
	cout << "dpdrho: " << thermo.dpdrho << endl;
	cout << "dpde: " << thermo.dpde << endl;
}

int main() {

	int n = 400;

	Vector e(n); 
	for (int i = 0; i < n; ++i) {
		e[i] = 717 * 600 + (3.5e7 - 717*600)/(n-1)*i; 
	}

	Vector rho = custom_rho_spacing();  

	Chemistry chem; 
	Vector mass_fractions;  
	Vector initial_moles = { 0.7808, 0.2095, 0.0, 0.0, 0.0, 0.0097, 0.0, 0.0, 0.0, 0.0 };  
	ofstream file("thermochemical_table.csv"); 
	file << "rho, e, p, T, R, cv, gam, dpdrho, dpde" << endl; // Header for the CSV file   
	Tensor thermo(n, Matrix(n, Vector(9)));  

	for (int i = 0; i < n; ++i) { 
		for (int j = 0; j < n; ++j) {  

			
			try {
				thermo[i][j] = chem.compute_equilibrium(rho[i], e[j], initial_moles); 
					// ... write line
			} catch (const exception& ex) {
				cerr << "Error at i=" << i << ", j=" << j << ": " << ex.what() << endl;
			} 			
			file << rho[i] << ", " << e[j] << ", " << thermo[i][j][2] << ", " << thermo[i][j][3] << ", " << thermo[i][j][4] 
				<< ", " << thermo[i][j][5] << ", " << thermo[i][j][6] << ", " << thermo[i][j][7] << ", " << thermo[i][j][8] << endl; // Write to CSV file  

		}
		cout << i * n << endl; 
	}


	
	file.close(); 
	return 0;
}
