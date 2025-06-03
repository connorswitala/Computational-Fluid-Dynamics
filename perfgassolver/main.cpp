
#include "../perfgaslib/perfgas.hpp"   

using namespace std;


int main() {

	// IMPORTANT:  Due to how this code is set up, you can change Nx (number of cells in 'x' direction) and Ny (number of cells in 'y' direction) in the 
	// 2D FVS Library.h file. It should be right at the top. 

	auto start = TIME;

	int Nx_holder, Ny_holder;

	double l1, l2, l3, l4, l5, R1, R2, R3, theta1, theta2, d_min;
	int progress_update;

	double Wall_Temp = 300; // For Isothermal Wall boundary condition only (viscous solver). 
	double CFL; // Pretty much only works at 1.0

	inlet_conditions INLET;
	unique_ptr<Grid> gridg;

	string solver_type, grid_type, leftBC, rightBC, bottomBC, topBC, preset;

	cout << "\033[36mRun the preset solver or enter your own (type 'preset' to solve defaule mirrored ramp grid, or 'set' to set your own grid, boundary conditions, and flow properties): \033[0m"; cin >> preset; cout << endl;

	while (preset != "preset" && preset != "set") {
		cout << "\033[36mUnknown entry. Type 'preset' for preset or 'set' for set your own: \033[0m"; cin >> preset; cout << endl;
	}


	if (preset == "set") {

		cout << "\033[36m'viscous' or 'inviscid' solver: \033[0m";	cin >> solver_type; cout << endl;

		while (solver_type != "viscous" && solver_type != "inviscid") {
			cout << "\033[36mUnknown solver type. Type 'viscous' for preset or 'inviscid' for set your own: \033[0m"; cin >> preset; cout << endl;
		}

		cout << "\033[36mGrid type (enter 'ramp' for ramp, 'cylinder' for cylinder, 'plate' for flat plate, 'double' for double ramp, \n 'mirrored' for double ramp that is mirrored above x-axis): \033[0m";	cin >> grid_type; cout << endl;

		while (grid_type != "ramp" && grid_type != "cylinder" && grid_type != "plate" && grid_type != "double" && grid_type != "mirrored") {
			cout << "\033[36mUnknown grid type. Please choose either 'ramp' 'cylinder' 'plate' or 'double': \033[0m"; cin >> grid_type; cout << endl;
		}

		cout << "\033[36mNumber of cells in x-direction: \033[0m"; cin >> Nx_holder; cout << "\033[36mNumber of cells in y-direction: \033[0m"; cin >> Ny_holder;

		const int Nx = Nx_holder, Ny = Ny_holder;

		cout << "\033[36mCFL (try to keep at or below 5.0): \033[0m";	cin >> CFL; cout << endl;

		if (grid_type == "ramp") {
			cout << "\033[31mInlet section length: \033[0m"; cin >> l1;  cout << "\033[31mRamp section length: \033[0m"; cin >> l2; cout << "\033[31mOutlet section length: \033[0m"; cin >> l3;
			cout << "\033[31mDomain height: \033[0m"; cin >> l4; cout << "\033[31mRamp angle (in degrees): \033[0m"; cin >> theta1;
			gridg = make_unique<RampGrid>(Nx, Ny, l1, l2, l3, l4, theta1);
		}

		if (grid_type == "cylinder") {
			cout << "\033[31mCylinder radius: \033[0m"; cin >> R1;  cout << "\033[31mMinimum grid radius: \033[0m"; cin >> R2; cout << "\033[31mMaximum grid radius: \033[0m"; cin >> R3; cout << "\033[31mMinumum cell thickness: \033[0m"; cin >> d_min;
			gridg = make_unique<CylinderGrid>(Nx, Ny, R1, R2, R3, d_min, pi / 2, 3 * pi / 2);
		}

		if (grid_type == "plate") {
			cout << "\033[31mx-length: \033[0m"; cin >> l1;  cout << "\033[31my-length: \033[0m"; cin >> l2; cout << "\033[31mMinimum cell thickness: \033[0m"; cin >> d_min;
			gridg = make_unique<FlatPlateGrid>(Nx, Ny, l1, l2, d_min);
		}

		if (grid_type == "double") {
			cout << "\033[31mInlet section length: \033[0m"; cin >> l1;  cout << "\033[31mFirst ramp section length: \033[0m"; cin >> l2; cout << "\033[31mSecond ramp section length: \033[0m"; cin >> l3; cout << "\033[31mOutlet section length: \033[0m"; cin >> l4;
			cout << "\033[31mInlet section height: \033[0m"; cin >> l5; cout << "\033[31mFirst ramp angle (in degrees): \033[0m"; cin >> theta1; cout << "\033[31mSecond ramp angle (in degrees): \033[0m"; cin >> theta2;
			gridg = make_unique<DoubleConeGrid>(Nx, Ny, l1, l2, l3, l4, theta1, theta2, l5);
		}

		if (grid_type == "mirrored") {
			cout << "\033[31mInlet section length: \033[0m"; cin >> l1;  cout << "\033[31mFirst ramp section length: \033[0m"; cin >> l2; cout << "\033[31mSecond ramp section length: \033[0m"; cin >> l3; cout << "\033[31mOutlet section length: \033[0m"; cin >> l4;
			cout << "\033[31mInlet section height: \033[0m"; cin >> l5; cout << "\033[31mFirst ramp angle (in degrees): \033[0m"; cin >> theta1; cout << "\033[31mSecond ramp angle (in degrees): \033[0m"; cin >> theta2;
			gridg = make_unique<MirroredGrid>(Nx, Ny, l1, l2, l3, l4, theta1, theta2, l5);
		}

		cout << endl << "\033[32mAvailable boundary conditions are: 'inlet', 'outlet', 'symmetry', 'adiabatic wall', 'isothermal wall' \nNOTE: isothermal and adiabatic are only necessary for the viscous solver" << endl << endl << "Left boundary condition: \033[0m";	cin >> leftBC;
		cout << "\033[32mRight boundary condition: \033[0m";	cin >> rightBC;
		cout << "\033[32mBottom boundary condition: \033[0m";	cin >> bottomBC;
		cout << "\033[32mTop boundary condition: \033[0m"; cin >> topBC;


		cout << endl << "\033[35mInput flow properties : \033[0m" << endl << "\033[35mPressure: \033[0m"; cin >> INLET.p;
		cout << "\033[35mTemperature: \033[0m"; cin >> INLET.T;	cout << "\033[35mMach: \033[0m"; cin >> INLET.M;


		cout << "\033[35mWall Temperature (only for isothermal wall boundary condition, enter any number if not necessary): \033[0m"; cin >> Wall_Temp; cout << endl;
		cout << "\033[36mAfter how many iterations do you want to see a status update: \033[0m"; cin >> progress_update; cout << endl;

		INLET.a = sqrt(gam * R * INLET.T),	// Inlet Sound Speed
			INLET.u = INLET.M * INLET.a,			// Inlet u-velocity
			INLET.v = 0,							// Inlet v-velocity
			INLET.rho = INLET.p / (R * INLET.T);	// Inlet density

		BoundaryConditions BCs(getBoundaryCondition(leftBC), getBoundaryCondition(rightBC), getBoundaryCondition(bottomBC), getBoundaryCondition(topBC));

		// Sets up the solver class with inlet and boundary conditions, grid type, CFL, wall temperature, and the iteration number for progress     
		Solver solver(Nx, Ny, INLET, *gridg, BCs, CFL, Wall_Temp, progress_update);

		// Calls solver type ( use -> solver.solve_viscous() if you want viscous solver)
		if (solver_type == "inviscid") solver.solve_inviscid();         
		else solver.solve_viscous();
	}

	else if (preset == "preset") {
	

		INLET.p = 10000.0,						// Inlet Pressure (SET)
			INLET.T = 300,						// Inlet Temperature (SET)
			INLET.M = 2.5,							// Inlet Mach speed (SET)
			INLET.a = sqrt(gam * R * INLET.T),	// Inlet Sound Speed
			INLET.u = INLET.M * INLET.a,			// Inlet u-velocity
			INLET.v = 0,							// Inlet v-velocity
			INLET.rho = INLET.p/(R * INLET.T);	// Inlet density

		int progress_update = 50;  // This number prints a status update after the number of iterations declared here. 
		const int Nx = 50, Ny = 25;
		CFL = 1.0; 
		BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);     
		RampGrid grid(Nx, Ny, 1, 1, 1, 0.65, 15);  

		// BoundaryConditions BCs(BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry, BoundaryCondition::Inlet);
		// CylinderGrid grid(Nx, Ny, 0.1, 0.3, 0.45, 0.001, pi, 3 * pi / 2);

		//BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::IsothermalWall, BoundaryCondition::Symmetry);   
		//FlatPlateGrid grid(Nx, Ny, 1e-3, 1e-3, 5e-6);  

		//BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);  
		//DoubleConeGrid grid(Nx, Ny, 1, 1, 1, 1, 25, 50, 2.5);  

		//BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);   
		//MirroredGrid grid(Nx, Ny, 1, 1, 1, 2, 15, 30, 2.5);   


		// Sets up the solver class with inlet and boundary conditions, grid type, CFL, wall temperature, and the iteration number for progress     
		Solver solver(Nx, Ny, INLET, grid, BCs, CFL, Wall_Temp, progress_update);

		// Calls solver type ( use -> solver.solve_viscous() if you want viscous solver)
		solver.solve_inviscid();


	}

	else {
		cout << "Unknown solver type";
	}

	cout << endl << "\033[36mBring newly created '.csv' file into Paraview for plotting. -> Click on the file in the pipeline. \nAt top of screen click Filters -> Alphabetical -> Tables to Points." << endl
		<< "Change x, y, and z column in properties to x, y, and z - points in dropdown. \nClick on Table to Points in pipeline. \nClick Filters -> Alphabetical -> Delauney2D \nIn properties section, change from surface to points." << endl
		<< "In coloring change to whatever property you want to view and click the little broken link square next to the Edit button. You can click edit to change coloring scheme. \033[0m";



	system("pause");
	return 0;
}
