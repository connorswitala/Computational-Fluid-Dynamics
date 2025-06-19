
#include "../realgaslib/realgas.hpp"

int main(int argc, char* argv[]) {

	// IMPORTANT:  Due to how this code is set up, you can change Nx (number of cells in 'x' direction) and Ny (number of cells in 'y' direction) in the 
	// 2D FVS Library.h file. It should be right at the top. 

	auto start = TIME;				// For timing the solution

	bool restart = false;  
	bool chemical_eqbm_enabled = true;
	string restart_name = "";
	string gridtype;   
	unique_ptr<Grid> gridg; 

	int Nx, Ny, progress_update = 50;
	double CFL = 1.0, Wall_Temp = 50;

	inlet_conditions INLET; 

	INLET.p = 66.2970,							// Inlet Pressure (SET)
	INLET.T = 231,							// Inlet Temperature (SET)
	INLET.M = 2.5,							// Inlet Mach speed (SET)
	INLET.a = sqrt(perfgam * perfR * INLET.T),		// Inlet Sound Speed
	INLET.u = INLET.M * INLET.a,			// Inlet u-velocity
	INLET.v = 0,							// Inlet v-velocity
	INLET.rho = 0.001;						// Inlet density  
	// INLET.rho = INLET.p/(R * INLET.T);

	Vector V_inlet = {INLET.rho, INLET.u, INLET.v, INLET.p};
	Vector U_inlet = primtoCons(V_inlet, perfgam); 

	if (argc == 1) {

		progress_update = 50;		// For displaying stats
		Wall_Temp = 300; 		// For Isothermal Wall boundary condition only (viscous solver). 
		CFL = 1.0; 				// Pretty much only works at 1.0
		Nx = 400, Ny = 200; 


		// BCMap BCs(BCType::Inlet, BCType::Outlet, BCType::Symmetry, BCType::Symmetry);      
		// gridg = make_unique<RampGrid>(Nx, Ny, 3, 0.66, 15);   
	

		BCMap BCs(BCType::Outlet, BCType::Symmetry, BCType::Symmetry, BCType::Inlet);   
		gridg = make_unique<CylinderGrid>(Nx, Ny, 0.1, 0.3, 0.45, 0.0001, pi / 2, pi);

		Tensor U = Tensor(Nx + 2, Matrix(Ny + 2, U_inlet));   

		
		if (dynamic_cast<RampGrid*>(gridg.get())) gridtype = "Ramp";
		else if (dynamic_cast<CylinderGrid*>(gridg.get())) gridtype = "Cylinder";
		else if (dynamic_cast<FlatPlateGrid*>(gridg.get())) gridtype = "Plate";
		else if (dynamic_cast<DoubleConeGrid*>(gridg.get())) gridtype = "Double";
		else if (dynamic_cast<MirroredGrid*>(gridg.get())) gridtype = "Mirrored";
		else gridtype = "Unknown"; 

		for (int i = 0; i < Nx + 2; ++i) {
			for (int j = 0; j < Ny + 2; ++j) {
				U[i][j] = U_inlet;
			}
		}

		Solver solver(Nx, Ny, U_inlet, *gridg, gridtype, BCs, CFL, Wall_Temp, progress_update, restart, restart_name, chemical_eqbm_enabled, U);
		solver.solve_inviscid();  
	}
	
	if (argc == 2) {

		restart = true;
		restart_name = "../plotfiles/CEI_" + string(argv[1]) + ".dat"; 
		auto [Nx, Ny, U, grid_ptr, gridname, U_inlet, BCs] = restart_solution(restart_name);
		Grid& grid = *grid_ptr;


		Solver solver(Nx, Ny, U_inlet, grid, gridname, BCs, CFL, Wall_Temp, progress_update, restart, restart_name, chemical_eqbm_enabled, U);
		solver.solve_inviscid();
	}
	if (argc == 3) {

		progress_update = stoi(argv[2]); 
		restart = true;
		restart_name = "../plotfiles/CEI_" + string(argv[1]) + ".dat"; 
		cout << restart_name << endl;
		auto [Nx, Ny, U, grid_ptr, gridname, U_inlet, BCs] = restart_solution(restart_name);
		Grid& grid = *grid_ptr;


		Solver solver(Nx, Ny, U_inlet, grid, gridname, BCs, CFL, Wall_Temp, progress_update, restart, restart_name, chemical_eqbm_enabled, U);
		solver.solve_inviscid();
	}

	return 0;
}


